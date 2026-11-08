from helpers import calculate_xtb
from dataclasses import replace
from pathlib import Path
import json
from unittest.mock import patch

import pytest

Chem = pytest.importorskip('rdkit.Chem')
from rdkit.Chem import AllChem
from eled import (
    XTBOptions, XTBResult, XTBRunError, XTBCache,
    AnnotateMol, RunXTBBatch, XTBCacheKey,
)
from eled.input import _input_metadata, _xyz_from_mol


def molecule(smiles='CCO'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    return mol


def sample_result(mol):
    n = mol.GetNumAtoms()
    return XTBResult(
        'GFN2-xTB', 'test', -1.0, 2.0, 2, 2, 0,
        orbital_energies_ev=(-2., 0.), orbital_occupations=(2., 0.),
        partial_charges=tuple(i * .01 for i in range(n)),
        wiberg_bond_orders=tuple(tuple(0. if i == j else .8 for j in range(n)) for i in range(n)),
        metadata=_input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1)),
    )


@pytest.fixture
def engine(tmp_path, monkeypatch):
    path = tmp_path / 'libxtb.so'
    path.write_text('test engine A')
    from types import SimpleNamespace
    import os
    monkeypatch.setattr('eled.xtb._find_xtb_library',
                        lambda explicit: str(explicit or os.environ['XTB_LIBRARY']))
    monkeypatch.setattr('eled.xtb._NativeXTB', lambda path: SimpleNamespace(library=path))
    # This fixture exercises cache-key plumbing without a real extended engine.
    monkeypatch.setattr('eled.xtb._require_extended_libxtb', lambda path: {})
    return path


def test_options_are_serializable_and_validate():
    opts = XTBOptions(method='gfn2', xtb_library=Path('/test/xtb'))
    assert XTBOptions(**json.loads(json.dumps(opts.to_dict()))) == opts
    assert opts.backend == 'native'
    with pytest.raises(ValueError, match='backend'):
        XTBOptions(backend='unknown')


def test_mapping_preserves_input_and_bond_types():
    mol = molecule()
    before = mol.ToBinary()
    result = sample_result(mol)
    target = AnnotateMol(mol, result)
    assert mol.ToBinary() == before
    assert [b.GetBondType() for b in target.GetBonds()] == [b.GetBondType() for b in mol.GetBonds()]
    assert target.GetAtomWithIdx(2).GetDoubleProp('xtb_partial_charge') == .02
    assert target.GetBondWithIdx(0).GetDoubleProp('xtb_wbo') == .8
    assert json.loads(target.GetProp('xtb_metadata_json'))['conformer_id'] == 0


def test_mapping_rejects_wrong_order_isomer_geometry_and_legacy_result():
    mol = molecule()
    result = sample_result(mol)
    permuted = Chem.RenumberAtoms(mol, [1, 0, *range(2, mol.GetNumAtoms())])
    for other in (permuted, molecule('COC')):
        with pytest.raises(ValueError, match='structure or atom order'):
            AnnotateMol(other, result)
    mol.GetConformer().SetAtomPosition(0, (0., 0., 0.))
    with pytest.raises(ValueError, match='coordinates or atom order'):
        AnnotateMol(mol, result)
    with pytest.raises(ValueError, match='provenance'):
        AnnotateMol(mol, replace(result, metadata={}))


def test_invalid_arrays_fail_before_inplace_mutation():
    mol = molecule()
    result = sample_result(mol)
    for damaged in (replace(result, partial_charges=(.1,)),
                    replace(result, wiberg_bond_orders=((float('nan'),),))):
        with pytest.raises(ValueError):
            AnnotateMol(mol, damaged, copy=False)
        assert not mol.GetAtomWithIdx(0).HasProp('xtb_partial_charge')


def test_atom_hydrogen_counts_are_not_explicit_hydrogen_vertices():
    mol = Chem.MolFromSmiles('[NH4+]')
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    mol.AddConformer(conf)
    with pytest.raises(ValueError, match='Chem.AddHs'):
        calculate_xtb(mol, xtb_library='/missing')


def test_coordinate_and_spin_validation_happens_before_engine_lookup():
    mol = molecule()
    with pytest.raises(ValueError, match='electron count'):
        calculate_xtb(mol, uhf=1, xtb_library='/missing')
    with pytest.raises(ValueError, match='integers'):
        calculate_xtb(mol, charge=.5, xtb_library='/missing')
    mol.GetConformer().Set3D(False)
    with pytest.raises(ValueError, match='3D'):
        calculate_xtb(mol, xtb_library='/missing')
    mol.GetConformer().Set3D(True)
    mol.GetConformer().SetAtomPosition(0, (float('nan'), 0., 0.))
    with pytest.raises(ValueError, match='finite'):
        calculate_xtb(mol, xtb_library='/missing')


def test_cache_keys_include_geometry_options_library_and_parameters(engine, monkeypatch):
    mol = molecule()
    base = XTBCacheKey(mol, xtb_library=engine)
    assert base == XTBCacheKey(mol, xtb_library=str(engine), method='gfn2', conf_id=0)
    assert base != XTBCacheKey(mol, xtb_library=engine, method='gfn1')
    assert base != XTBCacheKey(mol, xtb_library=engine, charge=1, uhf=1)
    mol.GetConformer().SetAtomPosition(0, (0., 0., 0.))
    assert base != XTBCacheKey(mol, xtb_library=engine)
    mol = molecule()
    engine.write_text('test engine B')
    assert base != XTBCacheKey(mol, xtb_library=engine)
    changed_binary = XTBCacheKey(mol, xtb_library=engine)
    parameter = engine.parent / 'share/xtb/param_gfn2-xtb.txt'
    parameter.parent.mkdir(parents=True)
    parameter.write_text('test parameters A')
    monkeypatch.setenv('XTBPATH', str(parameter.parent))
    with_parameter = XTBCacheKey(mol, xtb_library=engine)
    assert changed_binary != with_parameter
    parameter.write_text('test parameters B')
    assert with_parameter != XTBCacheKey(mol, xtb_library=engine)
    monkeypatch.setenv('XTB_LIBRARY', str(engine))
    assert XTBCacheKey(mol) == XTBCacheKey(mol, xtb_library=engine)


def test_batch_reuses_results_and_isolates_metadata(engine):
    mol = molecule()
    cache = XTBCache()
    with patch('eled.batch._run_xtb', return_value=sample_result(mol)) as calculate:
        results = RunXTBBatch([mol] * 12, cache=cache, xtb_library=engine)
    assert calculate.call_count == 1
    assert len(cache) == 1
    results[0].metadata['charge'] = 100
    assert results[1].metadata['charge'] == 0
    assert cache.get(XTBCacheKey(mol, xtb_library=engine)).metadata['charge'] == 0


def test_batch_keeps_order_and_retries_failed_calculations(engine):
    a, b = molecule(), molecule('CO')
    cache = XTBCache()
    def fake(mol, **options):
        if mol.GetNumAtoms() == b.GetNumAtoms():
            raise XTBRunError('test failure')
        return sample_result(mol)
    with patch('eled.batch._run_xtb', side_effect=fake):
        results = RunXTBBatch([a, b, a], cache=cache, xtb_library=engine,
                                return_exceptions=True)
    assert isinstance(results[1], XTBRunError)
    assert results[0] == results[2]
    assert len(cache) == 1
    with patch('eled.batch._run_xtb', return_value=sample_result(b)) as calculate:
        RunXTBBatch([b], cache=cache, xtb_library=engine)
    assert calculate.call_count == 1
    with pytest.raises(ValueError, match='max_workers'):
        RunXTBBatch([a], backend='native', max_workers=0)
    with pytest.raises(TypeError):
        RunXTBBatch([a], unsupported_option=True, cache=cache)


def test_gfnff_rejected_before_calculation():
    with pytest.raises(ValueError, match='no molecular orbitals'):
        calculate_xtb(molecule(), method='GFN-FF')


@pytest.mark.parametrize('backend', ['cli', 'auto', 'tblite'])
def test_removed_backends_fail_before_calculation(backend, monkeypatch):
    import subprocess
    def forbidden(*args, **kwargs):
        pytest.fail('removed backends must not launch a process or load libxtb')
    monkeypatch.setattr(subprocess, 'Popen', forbidden)
    monkeypatch.setattr('eled.xtb._find_xtb_library', forbidden)
    with pytest.raises(ValueError, match='removed'):
        calculate_xtb(molecule(), backend=backend)
    with pytest.raises(ValueError, match='removed'):
        RunXTBBatch([molecule()], backend=backend)


@pytest.mark.parametrize('option,value', [('xtb_path', '/bin/xtb'), ('keep_files', True),
                                        ('timeout', 30), ('optimize', True), ('solvent', 'water')])
def test_removed_options_are_not_silently_ignored(option, value):
    with pytest.raises(TypeError, match=option):
        calculate_xtb(molecule(), **{option: value})
    with pytest.raises(TypeError, match=option):
        XTBOptions(**{option: value})
    if option != 'timeout':  # wall-clock timeout now belongs to isolated batch execution
        with pytest.raises(TypeError, match=option):
            RunXTBBatch([molecule()], **{option: value})
