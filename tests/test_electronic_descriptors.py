"""Analytical contracts for descriptors; no quantum engine is needed."""
from dataclasses import replace
import json

import numpy as np
import pytest
from rdkit import Chem

from eled import (
    XTBResults, XTBOptions, CalcXTBDescriptors, CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors, CalcXTBFragmentDescriptors,
    CalcOrbitalPopulations, CalcFrontierDescriptors, DescriptorUnavailableError,
    AnnotateMol,
)
from eled.input import _input_metadata, _xyz_from_mol


@pytest.mark.parametrize('overlap, mapping', [
    (None, (0,)), (((1.0,),), None), ((), ()), (None, None),
])
def test_missing_population_inputs_raise_descriptor_unavailable(overlap, mapping):
    result = XTBResults('GFN2-xTB', None, -1., None, 1, 1, 1,
                        overlapMatrix=overlap, aoAtomMap=mapping)
    with pytest.raises(DescriptorUnavailableError):
        CalcOrbitalPopulations(result)


def legacy_compiled_result(metadata):
    from types import SimpleNamespace
    result_type = type('XTBResults', (SimpleNamespace,), {'__module__': 'eled.rdXTBTools'})
    return result_type(
        success=True, overlapRetained=False, inputProvenanceJSON=json.dumps(metadata),
        method='GFN2-xTB', initialization='independent', maxIterations=250,
        accuracy=1., electronicTemperature=300., effectiveElectronicTemperature=300.,
        libraryPath='test', librarySHA256='test', sourceCommit='test', apiVersion=1,
        parameterFile=None, totalEnergyHartree=-1., numAtoms=2, numOrbitals=2,
        numElectrons=2, numUnpairedElectrons=0,
        GetAtomicCharges=lambda: (0.2, -0.2), GetOrbitalEnergies=lambda: (-1., 1.),
        GetOrbitalOccupancies=lambda: (2., 0.), GetDipole=lambda: (0., 0., 0.),
        GetWibergBondOrders=lambda: ((0., 1.), (1., 0.)),
    )


@pytest.mark.parametrize('numbers', [None, (), (1, 17)])
def test_legacy_compiled_cm5_requires_atomic_numbers(numbers):
    from eled.electronic_descriptors import _as_result
    metadata = {'coordinates_angstrom': ((0., 0., 0.), (1.3, 0., 0.)),
                'atomic_numbers': numbers}
    result = _as_result(legacy_compiled_result(metadata))
    assert result.atomicCharges == (0.2, -0.2)
    if numbers:
        assert len(result.cm5Charges) == 2
        assert sum(result.cm5Charges) == pytest.approx(0.)
    else:
        assert result.cm5Charges == ()


def test_legacy_compiled_cm5_rejects_inconsistent_atom_counts():
    from eled.electronic_descriptors import _as_result
    result = legacy_compiled_result({'atomic_numbers': (1, 17),
                                     'coordinates_angstrom': ((0., 0., 0.),)})
    with pytest.raises(ValueError, match='matching atom counts'):
        _as_result(result)


@pytest.fixture
def simple():
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    conf = Chem.Conformer(3)
    for i, xyz in enumerate(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.))):
        conf.SetAtomPosition(i, xyz)
    conf.Set3D(True)
    mol.AddConformer(conf)
    # Synthetic arrays have intentionally nonzero diagonal WBOs, which must
    # never enter pair/atom sums, and an off-topology H-H pair of 0.2.
    result = XTBResults(
        'GFN2-xTB', 'test', -1., 3., 4, 4, 0,
        orbitalEnergiesEV=(-10., -2., 1., 3.), orbitalOccupancies=(2., 2., 0., 0.),
        atomicCharges=(-.6, .2, .4), dipoleAtomicUnits=(3., 4., 0.),
        wibergBondOrders=((9., .8, .9), (.8, 9., .2), (.9, .2, 9.)),
        metadata=_input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1)),
    )
    return mol, result


def test_scalar_units_frontier_depth_and_missing_boundaries(simple):
    _, result = simple
    d = CalcXTBDescriptors(result, frontierDepth=3)
    assert d['xtb_gap_ev'] == 3
    assert d['xtb_homo_m1_ev'] == -10
    assert d['xtb_lumo_p1_ev'] == 3
    assert d['xtb_homo_m2_ev'] is None
    assert d['xtb_lumo_p3_ev'] is None
    assert d['xtb_dip_ebohr'] == 5
    assert d['xtb_dip_debye'] == pytest.approx(12.708732365)
    assert d['xtb_q_sum_e'] == pytest.approx(0)
    assert d['xtb_q_sd_e'] == pytest.approx(np.sqrt(.56 / 3))
    assert d['xtb_wbo_sum'] == pytest.approx(1.9)
    assert json.loads(json.dumps(d, allow_nan=False)) == d
    assert not any('minus' in k for k in CalcXTBDescriptors(result, frontierDepth=0))
    # Fermi smearing must not move the frontier: a partially filled LUMO stays the LUMO.
    for occupations in ((2., 2., 1e-8, 0.), (2., 1.99, 0.01, 0.), (2., 1.6, 0.4, 0.)):
        smeared = CalcXTBDescriptors(replace(result, orbitalOccupancies=occupations))
        assert (smeared['xtb_homo_idx'], smeared['xtb_lumo_idx']) == (1, 2)
        assert smeared['xtb_gap_ev'] == pytest.approx(3.)
        assert 'xtb_occ_cut' not in smeared
    doublet = replace(result, numElectrons=3, numUnpairedElectrons=1, orbitalOccupancies=(2., 1., 0., 0.))
    assert (CalcXTBDescriptors(doublet)['xtb_homo_idx'], CalcXTBDescriptors(doublet)['xtb_lumo_idx']) == (1, 2)
    triplet = replace(result, numElectrons=4, numUnpairedElectrons=2, orbitalOccupancies=(2., 1., 1., 0.))
    assert (CalcXTBDescriptors(triplet)['xtb_homo_idx'], CalcXTBDescriptors(triplet)['xtb_lumo_idx']) == (2, 3)
    assert replace(result, orbitalOccupancies=(2., 1.6, 0.4, 0.)).homoEV == -2.
    assert replace(result, orbitalOccupancies=(2., 1.6, 0.4, 0.)).lumoEV == 1.
    filled = replace(result, numElectrons=8, orbitalOccupancies=(2., 2., 2., 2.))
    assert CalcXTBDescriptors(filled)['xtb_lumo_ev'] is None
    assert CalcXTBDescriptors(filled)['xtb_gap_ev'] is None


def test_missing_arrays_are_none_and_not_fake_zero(simple):
    _, result = simple
    missing = replace(result, atomicCharges=(), dipoleAtomicUnits=(), wibergBondOrders=())
    d = CalcXTBDescriptors(missing)
    assert d['xtb_q_sum_e'] is None
    assert d['xtb_dip_debye'] is None
    assert d['xtb_wbo_sum'] is None
    with pytest.raises(ValueError, match='partial charges'):
        CalcXTBAtomDescriptors(missing)


def test_atom_bond_and_fragment_accounting(simple):
    mol, result = simple
    atoms = CalcXTBAtomDescriptors(result)
    assert [a['wbo_sum'] for a in atoms] == pytest.approx([1.7, 1., 1.1])
    bonds = CalcXTBBondDescriptors(mol, result)
    assert [b['wbo'] for b in bonds] == [.8, .9]
    f = CalcXTBFragmentDescriptors(result, {'oxygen': [0], 'hydrogens': [1, 2]})
    assert f['fragments']['oxygen']['q_e'] == -.6
    assert f['fragments']['hydrogens']['q_e'] == pytest.approx(.6)
    assert f['fragments']['hydrogens']['wbo_intra'] == .2
    assert f['fragments']['oxygen']['wbo_inter'] == pytest.approx(1.7)
    assert f['pairs'][0]['wbo'] == pytest.approx(1.7)
    assert not f['unassigned_atom_indices']
    partial = CalcXTBFragmentDescriptors(result, {'OH': [0, 1]})
    assert partial['unassigned_atom_indices'] == [2]
    assert partial['fragments']['OH']['wbo_inter'] == pytest.approx(1.1)
    assert json.loads(json.dumps(f)) == f
    altered = Chem.RenumberAtoms(mol, [1, 0, 2])
    with pytest.raises(ValueError, match='structure or atom order'):
        CalcXTBBondDescriptors(altered, result)


@pytest.mark.parametrize('fragments', [{}, {'x': []}, {'x': [0, 0]}, {'x': [-1]},
    {'x': [3]}, {'x': [True]}, {'x': [1.5]}, {'': [0]}, {'a': [0], 'b': [0, 1]}])
def test_fragment_indices_fail_explicitly(simple, fragments):
    with pytest.raises(ValueError):
        CalcXTBFragmentDescriptors(simple[1], fragments)


@pytest.mark.parametrize('change', [
    {'dipoleAtomicUnits': (1., 2.)}, {'atomicCharges': (0.,)},
    {'atomicCharges': (0., float('nan'), 0.)},
    {'wibergBondOrders': ((0., 1., 0.), (0., 0., 0.), (0., 0., 0.))},
    {'orbitalEnergiesEV': (-2., -10., 1., 3.)}, {'orbitalOccupancies': (2.,)},
])
def test_invalid_data_rejected_before_inplace_annotation(simple, change):
    mol, result = simple
    before = mol.ToBinary()
    with pytest.raises(ValueError):
        AnnotateMol(mol, replace(result, **change), copy=False)
    assert mol.ToBinary() == before


@pytest.fixture
def basis(simple):
    mol, result = simple
    # Nonorthogonal SPD metric with four AOs on three atoms; construct C so
    # Lowdin amplitudes equal a known orthogonal Q, independent of our code.
    s = np.array([[1., .3, 0., 0.], [.3, 1., .1, 0.],
                  [0., .1, 1., .2], [0., 0., .2, 1.]])
    eigenvalues, vectors = np.linalg.eigh(s)
    q = np.array([[1., 1., 1., 1.], [1., -1., 1., -1.],
                  [1., 1., -1., -1.], [1., -1., -1., 1.]]) / 2
    c = ((vectors / np.sqrt(eigenvalues)) @ vectors.T) @ q
    return mol, replace(result, overlapMatrix=tuple(map(tuple, s)),
                        orbitalCoefficients=tuple(map(tuple, c)), aoAtomMap=(0, 0, 1, 2)), q


def test_lowdin_and_mulliken_use_the_metric_and_conserve_electrons(basis):
    _, result, q = basis
    lowdin = CalcOrbitalPopulations(result)
    np.testing.assert_allclose(lowdin['atom_pop'], [[.5]*4, [.25]*4, [.25]*4], atol=1e-12)
    np.testing.assert_allclose(lowdin['mo_n_e'], [2., 1., 1.], atol=1e-12)
    c = np.array(result.orbitalCoefficients); s = np.array(result.overlapMatrix)
    mulliken = CalcOrbitalPopulations(result, scheme='mulliken', orbIdx=[3, 0])
    expected = c * (s @ c)
    expected = np.vstack([expected[:2].sum(axis=0), expected[2:]])[:, [3, 0]]
    np.testing.assert_allclose(mulliken['atom_pop'], expected, atol=1e-12)
    assert sum(mulliken['mo_n_e']) == pytest.approx(2.)
    # Squared raw AO coefficients are not the Lowdin answer.
    assert not np.allclose(c**2, q**2)
    changed = replace(result, orbitalCoefficients=tuple(map(tuple, c * [-1, 1, -1, 1])))
    np.testing.assert_allclose(CalcOrbitalPopulations(changed)['atom_pop'], lowdin['atom_pop'], atol=1e-12)


def test_degenerate_subspace_and_spatial_invariances(basis):
    mol, result, q = basis
    c = np.array(result.orbitalCoefficients) @ q.T
    result = replace(result, orbitalEnergiesEV=(-2., -2., 1., 1.),
                     orbitalCoefficients=tuple(map(tuple, c)))
    baseline = CalcFrontierDescriptors(mol, result, fragments={'O': [0], 'H': [1, 2]})
    assert baseline['homo']['orb_idx'] == [0, 1]
    assert baseline['lumo']['orb_idx'] == [2, 3]
    assert baseline['homo']['pr'] == pytest.approx(1.)
    assert baseline['homo']['frag_pop']['H'] == pytest.approx(0., abs=1e-12)
    assert baseline['hl_dist_a'] == pytest.approx(np.sqrt(.5))
    c = np.array(result.orbitalCoefficients)
    rotation = np.array([[.6, -.8], [.8, .6]])
    c[:, :2] = c[:, :2] @ rotation
    changed = replace(result, orbitalCoefficients=tuple(map(tuple, c)))
    after = CalcFrontierDescriptors(mol, changed)
    np.testing.assert_allclose(after['homo']['atom_pop'], baseline['homo']['atom_pop'], atol=1e-12)
    # Translation/rotation affect centroids, not distances or populations.
    moved = Chem.Mol(mol)
    spatial_rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    positions = mol.GetConformer().GetPositions() @ spatial_rotation + [2., 3., 4.]
    for i, xyz in enumerate(positions): moved.GetConformer().SetAtomPosition(i, xyz)
    with pytest.raises(ValueError, match='coordinates'):
        CalcFrontierDescriptors(moved, result)
    moved_result = replace(result, metadata=_input_metadata(moved, -1, 0, 0, _xyz_from_mol(moved, -1)))
    spatial = CalcFrontierDescriptors(moved, moved_result)
    np.testing.assert_allclose(spatial['homo']['cent_a'],
                               np.array(baseline['homo']['cent_a']) @ spatial_rotation + [2., 3., 4.])
    assert spatial['hl_dist_a'] == pytest.approx(baseline['hl_dist_a'], abs=1e-12)
    json.dumps(spatial, allow_nan=False)


@pytest.mark.parametrize('change', [
    {'aoAtomMap': (0, 0, 1, 5)}, {'aoAtomMap': (0, 0, 0, 0)},
    {'overlapMatrix': tuple(map(tuple, -np.eye(4)))},
    {'orbitalCoefficients': tuple(map(tuple, np.eye(4)))},
])
def test_bad_or_mismatched_basis_fails(basis, change):
    with pytest.raises(ValueError):
        CalcOrbitalPopulations(replace(basis[1], **change))


def test_missing_basis_has_no_implicit_fallback(simple, monkeypatch):
    def forbidden(*args, **kwargs): pytest.fail('descriptor extraction must not run a calculation')
    monkeypatch.setattr('eled.xtb.run_xtb_native', forbidden)
    with pytest.raises(DescriptorUnavailableError, match='same-calculation'):
        CalcOrbitalPopulations(simple[1])
    assert XTBOptions(keepOverlap=True).ToDict()['keepOverlap'] is True


def test_mulliken_signed_populations_are_not_clipped(simple):
    _, template = simple
    overlap = np.array([[1., .8], [.8, 1.]])
    first = np.array([1., -.2]); first /= np.sqrt(first @ overlap @ first)
    second = np.array([0., 1.]); second -= first * (first @ overlap @ second)
    second /= np.sqrt(second @ overlap @ second)
    c = np.column_stack([first, second])
    result = replace(template, numOrbitals=2,
                     numElectrons=2, orbitalEnergiesEV=(-1., 1.),
                     orbitalOccupancies=(2., 0.), atomicCharges=(0., 0.),
                     orbitalCoefficients=tuple(map(tuple, c)), overlapMatrix=tuple(map(tuple, overlap)),
                     aoAtomMap=(0, 1), metadata={'atom_count': 2})
    populations = CalcOrbitalPopulations(result, scheme='mulliken', orbIdx=[0])
    assert populations['atom_pop'][1][0] < 0
    assert np.sum(populations['atom_pop']) == pytest.approx(1.)
    assert np.min(CalcOrbitalPopulations(result)['atom_pop']) >= 0


def test_population_ao_permutation_preserves_atom_results(basis):
    _, result, _ = basis
    permutation = [2, 0, 3, 1]
    c, s, mapping = np.array(result.orbitalCoefficients), np.array(result.overlapMatrix), np.array(result.aoAtomMap)
    permuted = replace(result, orbitalCoefficients=tuple(map(tuple, c[permutation])),
                       overlapMatrix=tuple(map(tuple, s[np.ix_(permutation, permutation)])),
                       aoAtomMap=tuple(mapping[permutation]))
    for scheme in ['lowdin', 'mulliken']:
        np.testing.assert_allclose(CalcOrbitalPopulations(result, scheme=scheme)['atom_pop'],
                                   CalcOrbitalPopulations(permuted, scheme=scheme)['atom_pop'], atol=1e-12)


def test_annotation_includes_wbo_sum_and_prevents_property_collision(simple):
    mol, result = simple
    target = AnnotateMol(mol, result)
    assert target.GetAtomWithIdx(0).GetDoubleProp('xtb_wbo_sum') == pytest.approx(1.7)
    assert target.GetDoubleProp('xtb_dip_debye') == pytest.approx(12.708732365)
    assert not mol.GetAtomWithIdx(0).HasProp('xtb_wbo_sum')
    with pytest.raises(ValueError, match='must differ'):
        AnnotateMol(mol, result, copy=False, atomChargeProp='same', atomWibergProp='same')
    assert not mol.GetAtomWithIdx(0).HasProp('same')
