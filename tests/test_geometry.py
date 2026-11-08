"""Actual coordinate stereo, failure boundaries and deterministic preparation."""
from helpers import calculate_xtb
import json
import os

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import geometry, prepare_conformers


FAST = dict(fallback_methods=(), attempts_per_strategy=1, retry_random_coords=False)


def test_deterministic_owned_conformer_pool_and_partial_coverage():
    source = Chem.MolFromSmiles('CCO'); before = source.ToBinary()
    a = prepare_conformers(source, num_conformers=8, **FAST)
    b = prepare_conformers(source, num_conformers=8, **FAST)
    assert source.ToBinary() == before
    assert a.success and 0 < a.molecule.GetNumConformers() <= 8
    assert a.report == b.report
    for ca, cb in zip(a.molecule.GetConformers(), b.molecule.GetConformers()):
        assert np.array_equal(ca.GetPositions(), cb.GetPositions())
    assert a.report['n_missing'] + a.report['n_conformers'] == 8
    assert a.report['status'] == ('partial' if a.report['n_missing'] else 'complete')
    assert len(a.report['attempts']) == 1
    assert all(not a.GetNumImplicitHs() and not a.GetNumExplicitHs() for a in a.molecule.GetAtoms())
    assert all(r['optimization_status'] == 0 for r in a.report['attempts'][0]['candidates'] if r['accepted'])
    assert json.loads(json.dumps(a.report)) == a.report


@pytest.mark.parametrize('smiles', ['C[C@H](O)C(=O)O', 'F/C=C/F', 'F/C=C\\F'])
def test_coordinates_preserve_specified_tetrahedral_and_double_bond_stereo(smiles):
    result = prepare_conformers(smiles, **FAST)
    assert result.success
    observed = Chem.Mol(result.molecule); Chem.RemoveStereochemistry(observed)
    Chem.AssignStereochemistryFrom3D(observed)
    assert Chem.MolToSmiles(Chem.RemoveHs(observed)) == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_reflected_coordinates_rejected_even_with_original_chiral_tags():
    r = prepare_conformers('C[C@H](O)C(=O)O', optimize=None, **FAST)
    mol = r.molecule; expected = geometry._stereo_constraints(mol)
    conf = mol.GetConformer()
    for i, xyz in enumerate(conf.GetPositions()):
        xyz[0] *= -1; conf.SetAtomPosition(i, xyz)
    assert geometry._check_geometry(mol, 0, expected)['reason'] == 'stereo_mismatch'


def test_optimizer_inversion_is_not_returned_as_a_success(monkeypatch):
    def invert(mol, **kwargs):
        c = mol.GetConformer()
        for i, xyz in enumerate(c.GetPositions()):
            xyz[0] *= -1; c.SetAtomPosition(i, xyz)
        return [(0, -1.)]
    monkeypatch.setattr(AllChem, 'MMFFOptimizeMoleculeConfs', invert)
    r = prepare_conformers('C[C@H](O)C(=O)O', **FAST)
    assert not r.success and r.molecule.GetNumConformers() == 0
    assert r.report['attempts'][0]['candidates'][0]['after_optimization']['reason'] == 'stereo_mismatch'


def test_nonconverged_force_field_is_not_returned(monkeypatch):
    monkeypatch.setattr(AllChem, 'MMFFOptimizeMoleculeConfs', lambda *a, **k: [(1, 1.)])
    result = prepare_conformers('CCO', **FAST)
    assert not result.success
    assert result.report['attempts'][0]['candidates'][0]['rejection'] == 'force_field_not_converged'


def test_missing_parameters_reported_before_embedding(monkeypatch):
    monkeypatch.setattr(AllChem, 'MMFFHasAllMoleculeParams', lambda mol: False)
    def forbidden(*a, **k): raise AssertionError('must not embed')
    monkeypatch.setattr(AllChem, 'EmbedMultipleConfs', forbidden)
    result = prepare_conformers('CCO')
    assert not result.success and not result.report['attempts']
    assert result.report['error']['reason'] == 'force_field_parameters_unavailable'


def test_retries_use_fresh_seeds_and_keep_stereo(monkeypatch):
    embed = AllChem.EmbedMultipleConfs; seen = []
    def fail_first(mol, **kwargs):
        p = kwargs['params']; seen.append((p.randomSeed, p.useRandomCoords, p.enforceChirality))
        return [] if len(seen) < 3 else embed(mol, **kwargs)
    monkeypatch.setattr(AllChem, 'EmbedMultipleConfs', fail_first)
    r = prepare_conformers('CCO', fallback_methods=())
    assert r.success and len(seen) == 3
    assert len({s[0] for s in seen}) == 3
    assert [s[1] for s in seen] == [False, False, True]
    assert all(s[2] for s in seen)


def test_unexpected_errors_propagate(monkeypatch):
    def broken(*a, **k): raise RuntimeError('unexpected backend bug')
    monkeypatch.setattr(AllChem, 'EmbedMultipleConfs', broken)
    with pytest.raises(RuntimeError, match='unexpected'): prepare_conformers('CCO')


def test_unavailable_default_preset_is_explicit_but_explicit_request_raises(monkeypatch):
    monkeypatch.delattr(AllChem, 'ETKDGv4', raising=False)
    r = prepare_conformers('CCO')
    assert r.success and r.report['retry_policy']['unavailable_default_fallbacks'] == ['ETKDGv4']
    with pytest.raises(ValueError, match='lacks'): prepare_conformers('CCO', fallback_methods=['ETKDGv4'])


def test_unspecified_stereo_is_flagged_or_rejected():
    r = prepare_conformers('CC(O)C(=O)O', **FAST)
    assert r.success and r.report['unspecified_stereo'] and r.report['warnings']
    with pytest.raises(ValueError, match='unspecified'):
        prepare_conformers('CC(O)C(=O)O', require_assigned_stereo=True)


def test_input_coordinates_are_discarded_from_a_copy():
    mol = Chem.AddHs(Chem.MolFromSmiles('O')); AllChem.EmbedMolecule(mol, randomSeed=1)
    before = mol.ToBinary(); r = prepare_conformers(mol, **FAST)
    assert r.success and r.report['discarded_input_conformers'] == 1
    assert not r.report['source_coordinates_used'] and mol.ToBinary() == before


@pytest.mark.parametrize('kwargs', [{'num_conformers': 0}, {'num_conformers': True},
    {'attempts_per_strategy': 0}, {'random_seed': -1}, {'random_seed': 2**32},
    {'max_embedding_iterations': 0}, {'optimization_iterations': 1.5}, {'optimize': 'other'},
    {'method': 'MolFromSmiles'}, {'fallback_methods': ['typo']}, {'fallback_methods': 'ETKDGv4'}])
def test_invalid_configuration_rejected_before_embedding(monkeypatch, kwargs):
    def forbidden(*a, **k): raise AssertionError('must not embed')
    monkeypatch.setattr(AllChem, 'EmbedMultipleConfs', forbidden)
    with pytest.raises(ValueError): prepare_conformers('CCO', **kwargs)


@pytest.mark.parametrize('source', ['', '*C', '[Na+].[Cl-]', 'invalid smiles'])
def test_unsupported_inputs_raise(source):
    with pytest.raises(ValueError): prepare_conformers(source)


def test_raw_failure_counters_and_exhausted_plan_are_retained():
    r = prepare_conformers('C[C@H]1[C@H]2[C@@H]3C[C@]3(C)N21', optimize=None, **FAST)
    assert not r.success
    counts = r.report['attempts'][0]['failure_counts']
    assert counts['available'] and sum(counts['raw']) == sum(counts['named'].values()) > 0
    assert r.report['error']['reason'] == 'no_accepted_conformer'


@pytest.mark.integration
def test_prepared_molecule_connects_to_existing_xtb_api():
    if os.environ.get('ELED_REQUIRE_RDXTB') != '1': pytest.skip('native integration profile')
    r = prepare_conformers('O', **FAST)
    result = calculate_xtb(r.molecule, conf_id=r.report['recommended_conformer_id'])
    assert result.metadata['geometry_sha256'] == r.report['conformers'][0]['geometry_sha256']
