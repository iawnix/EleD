"""Real-engine cached/resumed single points and reactivity, including new processes."""
from dataclasses import replace
import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import (XTBDiskCache,ReactivityJob,XTBOptions,run_xtb_batch,
    run_reactivity_batch,XTBNativeNotFoundError,DescriptorUnavailableError)
from rdkit_eth.runtime import _find_xtb_library,_extended_libxtb_build,_native_basis_build

pytestmark=pytest.mark.integration


@pytest.fixture
def library():
    try:lib=_find_xtb_library()
    except XTBNativeNotFoundError:pytest.skip('extended libxtb unavailable')
    if not _extended_libxtb_build(lib):pytest.skip('reviewed extended libxtb required')
    return lib


def molecule():
    mol=Chem.AddHs(Chem.MolFromSmiles('O'));assert AllChem.EmbedMolecule(mol,randomSeed=1)==0
    return mol


def test_real_single_point_disk_cache_and_process_roundtrip(library,tmp_path):
    mol=molecule();path=tmp_path/'cache.sqlite';cache=XTBDiskCache(path)
    opts=XTBOptions(xtb_library=library,keep_overlap=bool(_native_basis_build(library)))
    first=run_xtb_batch([mol],options=opts,cache=cache)[0]
    with patch('rdkit_eth.batch._run_xtb',side_effect=AssertionError('recalculated')):
        again=run_xtb_batch([mol],options=opts,cache=XTBDiskCache(path))[0]
    assert first==again
    from rdkit_eth.batch import _prepare,_key
    config,metadata=_prepare(mol,opts);key=_key(config,metadata)
    source = Path(__file__).resolve().parents[1] / 'src'
    command=(f'import sys; sys.path.insert(0, {str(source)!r}); '
             'from rdkit_eth import XTBDiskCache; import sys,json; '
             'r=XTBDiskCache(sys.argv[1]).get(sys.argv[2]); '
             'print(json.dumps([r.total_energy_hartree,r.partial_charges]))')
    child=subprocess.run([sys.executable,'-I','-c',command,str(path),key],capture_output=True,text=True,check=True)
    energy,charges=json.loads(child.stdout)
    assert energy==first.total_energy_hartree and charges==list(first.partial_charges)


@pytest.mark.parametrize('interface',['native','compiled'])
def test_real_reactivity_resume_preserves_failures_and_requires_capabilities(library,tmp_path,interface):
    if interface=='compiled':pytest.importorskip('rdkit_eth.rdXTBTools')
    opts=XTBOptions(xtb_library=library)
    jobs=[ReactivityJob('water',molecule(),opts),
          ReactivityJob('forced-scf-failure',molecule(),replace(opts,max_iterations=1))]
    rows=run_reactivity_batch(jobs,interface=interface,cache=XTBDiskCache(tmp_path/'cache.sqlite'),checkpoint=tmp_path/'rows.sqlite')
    assert rows[0]['error'] is None and rows[1]['error']['type']=='ChargeStateCalculationError'
    with patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=AssertionError('recalculated')):
        resumed=run_reactivity_batch(jobs,interface=interface,checkpoint=tmp_path/'rows.sqlite')
    assert resumed[0]['values']==rows[0]['values'] and resumed[1]['error']==rows[1]['error']
    if not _native_basis_build(library):
        with pytest.raises(DescriptorUnavailableError):
            run_reactivity_batch(jobs,interface=interface,population_scheme='lowdin')
