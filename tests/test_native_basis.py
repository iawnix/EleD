"""Native same-calculation basis and population contracts; no second engine."""
import builtins
from dataclasses import replace
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit_eth import (
    run_xtb, run_xtb_batch, XTBCache, DescriptorUnavailableError,
    calc_orbital_populations, calc_frontier_descriptors, calc_xtb_bond_descriptors, annotate_mol,
)
from rdkit_eth.runtime import _find_xtb_library, _native_basis_build, _BASIS_EXPORT_PATCH

CASES = json.loads((Path(__file__).parent/'data/orbital_alignment_cases.json').read_text())['cases']


@pytest.fixture
def library():
    from rdkit_eth import XTBNativeNotFoundError
    try: path = _find_xtb_library()
    except XTBNativeNotFoundError: pytest.skip('libxtb not installed')
    if not _native_basis_build(path):
        if os.environ.get('RDKIT_ETH_EXPECT_BASIS') == '1': pytest.fail('verified native basis extension required')
        pytest.skip('native basis extension unavailable')
    return path


def molecule(smiles='O'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol,randomSeed=1)==0
    return mol


@pytest.mark.integration
@pytest.mark.parametrize('method',['GFN1-xTB','GFN2-xTB'])
@pytest.mark.parametrize('case',CASES,ids=[c['name'] for c in CASES])
def test_native_basis_on_frozen_corpus(method,case,library):
    mol=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
    conf=Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
    for i,xyz in enumerate(case['positions_angstrom']): conf.SetAtomPosition(i,xyz)
    mol.AddConformer(conf)
    result=run_xtb(mol,method=method,xtb_library=library,keep_overlap=True)
    c,s=np.array(result.orbital_coefficients),np.array(result.overlap_matrix)
    np.testing.assert_allclose(c.T@s@c,np.eye(len(c)),atol=1e-8,rtol=0)
    ao,ash,sa,l=map(np.asarray,[result.ao_atom_indices,result.ao_shell_indices,result.shell_atom_indices,result.shell_angular_momenta])
    np.testing.assert_array_equal(sa[ash],ao)
    np.testing.assert_array_equal(np.bincount(ash),2*l+1)
    assert set(ao)==set(range(mol.GetNumAtoms()))
    for scheme in ['lowdin','mulliken']:
        pop=calc_orbital_populations(result,scheme=scheme)
        assert sum(pop['electron_populations_e'])==pytest.approx(result.number_of_electrons,abs=1e-7)
        np.testing.assert_allclose(np.sum(pop['atom_populations'],axis=0),1.,atol=1e-10)
    assert result.metadata['overlap_source']=='same libxtb calculation'
    # Native density populations must recover independent engine charges.
    valence={1:1,6:4,7:5,8:6,9:7,14:4,15:5,16:6,17:7,35:7,53:7}
    if all(a.GetAtomicNum() in valence for a in mol.GetAtoms()):
        reference=[valence[a.GetAtomicNum()] for a in mol.GetAtoms()]
        np.testing.assert_allclose(np.asarray(reference)-pop['electron_populations_e'],result.partial_charges,atol=1e-7,rtol=0)


@pytest.mark.integration
@pytest.mark.parametrize('method',['gfn0','gfn1','gfn2'])
@pytest.mark.parametrize('smiles',['CCO','[CH3]','[NH4+]'])
def test_all_methods_radicals_and_compiled_basis(method,smiles,library):
    rdXTBTools = pytest.importorskip("rdkit_eth.rdXTBTools")
    m=molecule(smiles)
    r=run_xtb(m,method=method,xtb_library=library,keep_overlap=True)
    ok,compiled=rdXTBTools.RunMol(m,method=method,xtbLibrary=library,keepOverlap=True)
    assert ok and compiled.overlapRetained and compiled.basisExportVerified
    np.testing.assert_allclose(compiled.GetOverlapMatrix(),r.overlap_matrix,atol=1e-12,rtol=0)
    for getter,field in [('GetAOAtomMap','ao_atom_indices'),('GetAOShellMap','ao_shell_indices'),
                         ('GetShellAtomMap','shell_atom_indices'),('GetShellAngularMomenta','shell_angular_momenta')]:
        value=getattr(compiled,getter)(); expected=getattr(r,field)
        assert value.dtype.kind=='i' and value.flags.owndata
        np.testing.assert_array_equal(value,expected)
        value[:]=-1
        np.testing.assert_array_equal(getattr(compiled,getter)(),expected)
    for scheme in ['lowdin','mulliken']:
        np.testing.assert_allclose(calc_orbital_populations(compiled,scheme=scheme)['atom_populations'],
                                   calc_orbital_populations(r,scheme=scheme)['atom_populations'],atol=1e-8)
    f=calc_frontier_descriptors(m,compiled)
    assert f['homo'] is not None
    assert calc_xtb_bond_descriptors(m,compiled)
    assert annotate_mol(m,compiled).GetAtomWithIdx(0).HasProp('xtb_wiberg_sum')
    old=compiled.GetOverlapMatrix().copy()
    run_xtb(molecule('CC'),xtb_library=library,keep_overlap=True)
    np.testing.assert_array_equal(compiled.GetOverlapMatrix(),old)
    m.GetConformer().SetAtomPosition(0,(12.,0.,0.))
    with pytest.raises(ValueError,match='coordinates'): calc_frontier_descriptors(m,compiled)


@pytest.mark.integration
def test_opt_in_cache_identity_and_no_other_engine(library,monkeypatch):
    original=builtins.__import__
    def guarded(name,*args,**kwargs):
        if name.split('.')[0] in ('tblite','ase'): pytest.fail('external engine imported')
        return original(name,*args,**kwargs)
    monkeypatch.setattr(builtins,'__import__',guarded)
    m=molecule(); cache=XTBCache()
    a=run_xtb_batch([m],cache=cache,xtb_library=library)[0]
    b=run_xtb_batch([m],cache=cache,xtb_library=library,keep_overlap=True)[0]
    assert len(cache)==2 and not a.overlap_matrix and b.overlap_matrix
    assert a.total_energy_eh==pytest.approx(b.total_energy_eh,abs=1e-12)
    with pytest.raises(DescriptorUnavailableError): calc_orbital_populations(a)
    assert calc_frontier_descriptors(m,b)['homo']


@pytest.mark.integration
def test_benzene_degenerate_subspaces(library):
    m=Chem.AddHs(Chem.MolFromSmiles('c1ccccc1')); conf=Chem.Conformer(12); conf.Set3D(True)
    for i in range(6):
        for index,radius in ((i,1.4),(i+6,2.48)):
            conf.SetAtomPosition(index,(radius*np.cos(i*np.pi/3),radius*np.sin(i*np.pi/3),0.))
    m.AddConformer(conf)
    r=run_xtb(m,xtb_library=library,keep_overlap=True)
    baseline=calc_frontier_descriptors(m,r); c=np.array(r.orbital_coefficients)
    for label in ['homo','lumo']:
        selected=baseline[label]['orbital_indices']; assert len(selected)==2
        c[:,selected]=c[:,selected]@np.array([[.6,-.8],[.8,.6]])
    changed=calc_frontier_descriptors(m,replace(r,orbital_coefficients=tuple(map(tuple,c))))
    for label in ['homo','lumo']:
        np.testing.assert_allclose(baseline[label]['atom_populations'],changed[label]['atom_populations'],atol=1e-12)


@pytest.mark.parametrize('damage',['missing','capability','patch','hash','binary','version','malformed'])
def test_capability_requires_exact_verified_binary(tmp_path,damage):
    lib=tmp_path/'lib/libxtb.so'; lib.parent.mkdir(); lib.write_bytes(b'test-binary')
    manifest=tmp_path/'share/xtb/rdkit-eth/build.json'; manifest.parent.mkdir(parents=True)
    record={'schema_version':1,'capabilities':{'orbital_basis_api':1},
            'source':{'commit':'test','basis_extension':{'api_version':1,'patch_sha256':_BASIS_EXPORT_PATCH}},
            'files':{'lib/libxtb.so':hashlib.sha256(lib.read_bytes()).hexdigest()}}
    manifest.write_text(json.dumps(record)); assert _native_basis_build(str(lib))
    if damage=='missing': manifest.unlink()
    elif damage=='binary': lib.write_bytes(b'different')
    elif damage=='malformed': manifest.write_text('[]')
    else:
        if damage=='capability': record['capabilities']['orbital_basis_api']=0
        if damage=='patch': record['source']['basis_extension']['patch_sha256']='other'
        if damage=='hash': record['files']['lib/libxtb.so']='other'
        if damage=='version': record['source']['basis_extension']['api_version']=2
        manifest.write_text(json.dumps(record))
    assert _native_basis_build(str(lib)) is None


@pytest.mark.integration
def test_standard_library_basic_support_and_missing_capability():
    from rdkit_eth.runtime import _find_xtb_library
    from rdkit_eth import XTBNativeNotFoundError
    try: library=_find_xtb_library()
    except XTBNativeNotFoundError: pytest.skip('libxtb not installed')
    if _native_basis_build(library): pytest.skip('requires official or older library')
    m=molecule()
    assert run_xtb(m,xtb_library=library).total_energy_eh<0
    with pytest.raises(DescriptorUnavailableError,match='verified native basis'):
        run_xtb(m,xtb_library=library,keep_overlap=True)
    rdXTBTools = pytest.importorskip("rdkit_eth.rdXTBTools")
    with pytest.raises(ValueError,match='verified native basis'):
        rdXTBTools.RunMol(m,xtbLibrary=library,keepOverlap=True)
