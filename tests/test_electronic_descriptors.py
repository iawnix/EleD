"""Analytical contracts for descriptors; no quantum engine is needed."""
from dataclasses import replace
import json

import numpy as np
import pytest
from rdkit import Chem

from rdkit_eth import (
    XTBResult, XTBOptions, calc_xtb_descriptors, calc_xtb_atom_descriptors,
    calc_xtb_bond_descriptors, calc_xtb_fragment_descriptors,
    calc_orbital_populations, calc_frontier_descriptors, DescriptorUnavailableError,
    run_xtb, annotate_mol,
)
from rdkit_eth.input import _input_metadata, _xyz_from_mol


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
    result = XTBResult(
        'GFN2-xTB', 'test', -1., 3., 4, 4, 0,
        orbital_energies_ev=(-10., -2., 1., 3.), orbital_occupations=(2., 2., 0., 0.),
        partial_charges=(-.6, .2, .4), dipole_e_bohr=(3., 4., 0.),
        wiberg_bond_orders=((9., .8, .9), (.8, 9., .2), (.9, .2, 9.)),
        metadata=_input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1)),
    )
    return mol, result


def test_scalar_units_frontier_depth_and_missing_boundaries(simple):
    _, result = simple
    d = calc_xtb_descriptors(result, frontier_depth=3)
    assert d['xtb_gap_ev'] == 3
    assert d['xtb_homo_minus_1_ev'] == -10
    assert d['xtb_lumo_plus_1_ev'] == 3
    assert d['xtb_homo_minus_2_ev'] is None
    assert d['xtb_lumo_plus_3_ev'] is None
    assert d['xtb_dipole_norm_e_bohr'] == 5
    assert d['xtb_dipole_norm_debye'] == pytest.approx(12.708732365)
    assert d['xtb_charge_sum_e'] == pytest.approx(0)
    assert d['xtb_charge_std_e'] == pytest.approx(np.sqrt(.56 / 3))
    assert d['xtb_wiberg_pair_sum'] == pytest.approx(1.9)
    assert json.loads(json.dumps(d, allow_nan=False)) == d
    assert not any('minus' in k for k in calc_xtb_descriptors(result, frontier_depth=0))
    tiny = replace(result, orbital_occupations=(2., 2., 1e-8, 0.))
    assert calc_xtb_descriptors(tiny)['xtb_homo_index'] == 1
    assert calc_xtb_descriptors(tiny, occupation_threshold=1e-9)['xtb_homo_index'] == 2
    filled = replace(result, orbital_occupations=(2., 2., 2., 2.))
    assert calc_xtb_descriptors(filled)['xtb_lumo_ev'] is None
    assert calc_xtb_descriptors(filled)['xtb_gap_ev'] is None


def test_missing_arrays_are_none_and_not_fake_zero(simple):
    _, result = simple
    missing = replace(result, partial_charges=(), dipole_e_bohr=(), wiberg_bond_orders=())
    d = calc_xtb_descriptors(missing)
    assert d['xtb_charge_sum_e'] is None
    assert d['xtb_dipole_norm_debye'] is None
    assert d['xtb_wiberg_pair_sum'] is None
    with pytest.raises(ValueError, match='partial charges'):
        calc_xtb_atom_descriptors(missing)


def test_atom_bond_and_fragment_accounting(simple):
    mol, result = simple
    atoms = calc_xtb_atom_descriptors(result)
    assert [a['wiberg_sum'] for a in atoms] == pytest.approx([1.7, 1., 1.1])
    bonds = calc_xtb_bond_descriptors(mol, result)
    assert [b['wiberg_order'] for b in bonds] == [.8, .9]
    f = calc_xtb_fragment_descriptors(result, {'oxygen': [0], 'hydrogens': [1, 2]})
    assert f['fragments']['oxygen']['charge_e'] == -.6
    assert f['fragments']['hydrogens']['charge_e'] == pytest.approx(.6)
    assert f['fragments']['hydrogens']['internal_wiberg_sum'] == .2
    assert f['fragments']['oxygen']['external_wiberg_sum'] == pytest.approx(1.7)
    assert f['pairs'][0]['wiberg_sum'] == pytest.approx(1.7)
    assert not f['unassigned_atom_indices']
    partial = calc_xtb_fragment_descriptors(result, {'OH': [0, 1]})
    assert partial['unassigned_atom_indices'] == [2]
    assert partial['fragments']['OH']['external_wiberg_sum'] == pytest.approx(1.1)
    assert json.loads(json.dumps(f)) == f
    altered = Chem.RenumberAtoms(mol, [1, 0, 2])
    with pytest.raises(ValueError, match='structure or atom order'):
        calc_xtb_bond_descriptors(altered, result)


@pytest.mark.parametrize('fragments', [{}, {'x': []}, {'x': [0, 0]}, {'x': [-1]},
    {'x': [3]}, {'x': [True]}, {'x': [1.5]}, {'': [0]}, {'a': [0], 'b': [0, 1]}])
def test_fragment_indices_fail_explicitly(simple, fragments):
    with pytest.raises(ValueError):
        calc_xtb_fragment_descriptors(simple[1], fragments)


@pytest.mark.parametrize('change', [
    {'dipole_e_bohr': (1., 2.)}, {'partial_charges': (0.,)},
    {'partial_charges': (0., float('nan'), 0.)},
    {'wiberg_bond_orders': ((0., 1., 0.), (0., 0., 0.), (0., 0., 0.))},
    {'orbital_energies_ev': (-2., -10., 1., 3.)}, {'orbital_occupations': (2.,)},
])
def test_invalid_data_rejected_before_inplace_annotation(simple, change):
    mol, result = simple
    before = mol.ToBinary()
    with pytest.raises(ValueError):
        annotate_mol(mol, replace(result, **change), copy=False)
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
    return mol, replace(result, overlap_matrix=tuple(map(tuple, s)),
                        orbital_coefficients=tuple(map(tuple, c)), ao_atom_indices=(0, 0, 1, 2)), q


def test_lowdin_and_mulliken_use_the_metric_and_conserve_electrons(basis):
    _, result, q = basis
    lowdin = calc_orbital_populations(result)
    np.testing.assert_allclose(lowdin['atom_populations'], [[.5]*4, [.25]*4, [.25]*4], atol=1e-12)
    np.testing.assert_allclose(lowdin['electron_populations_e'], [2., 1., 1.], atol=1e-12)
    c = np.array(result.orbital_coefficients); s = np.array(result.overlap_matrix)
    mulliken = calc_orbital_populations(result, scheme='mulliken', orbital_indices=[3, 0])
    expected = c * (s @ c)
    expected = np.vstack([expected[:2].sum(axis=0), expected[2:]])[:, [3, 0]]
    np.testing.assert_allclose(mulliken['atom_populations'], expected, atol=1e-12)
    assert sum(mulliken['electron_populations_e']) == pytest.approx(2.)
    # Squared raw AO coefficients are not the Lowdin answer.
    assert not np.allclose(c**2, q**2)
    changed = replace(result, orbital_coefficients=tuple(map(tuple, c * [-1, 1, -1, 1])))
    np.testing.assert_allclose(calc_orbital_populations(changed)['atom_populations'], lowdin['atom_populations'], atol=1e-12)


def test_degenerate_subspace_and_spatial_invariances(basis):
    mol, result, q = basis
    c = np.array(result.orbital_coefficients) @ q.T
    result = replace(result, orbital_energies_ev=(-2., -2., 1., 1.),
                     orbital_coefficients=tuple(map(tuple, c)))
    baseline = calc_frontier_descriptors(mol, result, fragments={'O': [0], 'H': [1, 2]})
    assert baseline['homo']['orbital_indices'] == [0, 1]
    assert baseline['lumo']['orbital_indices'] == [2, 3]
    assert baseline['homo']['participation_ratio'] == pytest.approx(1.)
    assert baseline['homo']['fragment_populations']['H'] == pytest.approx(0., abs=1e-12)
    assert baseline['homo_lumo_atom_centroid_separation_angstrom'] == pytest.approx(np.sqrt(.5))
    c = np.array(result.orbital_coefficients)
    rotation = np.array([[.6, -.8], [.8, .6]])
    c[:, :2] = c[:, :2] @ rotation
    changed = replace(result, orbital_coefficients=tuple(map(tuple, c)))
    after = calc_frontier_descriptors(mol, changed)
    np.testing.assert_allclose(after['homo']['atom_populations'], baseline['homo']['atom_populations'], atol=1e-12)
    # Translation/rotation affect centroids, not distances or populations.
    moved = Chem.Mol(mol)
    spatial_rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    positions = mol.GetConformer().GetPositions() @ spatial_rotation + [2., 3., 4.]
    for i, xyz in enumerate(positions): moved.GetConformer().SetAtomPosition(i, xyz)
    with pytest.raises(ValueError, match='coordinates'):
        calc_frontier_descriptors(moved, result)
    moved_result = replace(result, metadata=_input_metadata(moved, -1, 0, 0, _xyz_from_mol(moved, -1)))
    spatial = calc_frontier_descriptors(moved, moved_result)
    np.testing.assert_allclose(spatial['homo']['atom_centroid_angstrom'],
                               np.array(baseline['homo']['atom_centroid_angstrom']) @ spatial_rotation + [2., 3., 4.])
    assert spatial['homo_lumo_atom_centroid_separation_angstrom'] == pytest.approx(baseline['homo_lumo_atom_centroid_separation_angstrom'], abs=1e-12)
    json.dumps(spatial, allow_nan=False)


@pytest.mark.parametrize('change', [
    {'ao_atom_indices': (0, 0, 1, 5)}, {'ao_atom_indices': (0, 0, 0, 0)},
    {'overlap_matrix': tuple(map(tuple, -np.eye(4)))},
    {'orbital_coefficients': tuple(map(tuple, np.eye(4)))},
])
def test_bad_or_mismatched_basis_fails(basis, change):
    with pytest.raises(ValueError):
        calc_orbital_populations(replace(basis[1], **change))


def test_missing_basis_has_no_implicit_fallback(simple, monkeypatch):
    def forbidden(*args, **kwargs): pytest.fail('descriptor extraction must not run a calculation')
    monkeypatch.setattr('rdkit_eth.xtb.run_xtb_native', forbidden)
    with pytest.raises(DescriptorUnavailableError, match='same-calculation'):
        calc_orbital_populations(simple[1])
    assert XTBOptions(keep_overlap=True).to_dict()['keep_overlap'] is True


def test_mulliken_signed_populations_are_not_clipped(simple):
    _, template = simple
    overlap = np.array([[1., .8], [.8, 1.]])
    first = np.array([1., -.2]); first /= np.sqrt(first @ overlap @ first)
    second = np.array([0., 1.]); second -= first * (first @ overlap @ second)
    second /= np.sqrt(second @ overlap @ second)
    c = np.column_stack([first, second])
    result = replace(template, number_of_molecular_orbitals=2,
                     number_of_electrons=2, orbital_energies_ev=(-1., 1.),
                     orbital_occupations=(2., 0.), partial_charges=(0., 0.),
                     orbital_coefficients=tuple(map(tuple, c)), overlap_matrix=tuple(map(tuple, overlap)),
                     ao_atom_indices=(0, 1), metadata={'atom_count': 2})
    populations = calc_orbital_populations(result, scheme='mulliken', orbital_indices=[0])
    assert populations['atom_populations'][1][0] < 0
    assert np.sum(populations['atom_populations']) == pytest.approx(1.)
    assert np.min(calc_orbital_populations(result)['atom_populations']) >= 0


def test_population_ao_permutation_preserves_atom_results(basis):
    _, result, _ = basis
    permutation = [2, 0, 3, 1]
    c, s, mapping = np.array(result.orbital_coefficients), np.array(result.overlap_matrix), np.array(result.ao_atom_indices)
    permuted = replace(result, orbital_coefficients=tuple(map(tuple, c[permutation])),
                       overlap_matrix=tuple(map(tuple, s[np.ix_(permutation, permutation)])),
                       ao_atom_indices=tuple(mapping[permutation]))
    for scheme in ['lowdin', 'mulliken']:
        np.testing.assert_allclose(calc_orbital_populations(result, scheme=scheme)['atom_populations'],
                                   calc_orbital_populations(permuted, scheme=scheme)['atom_populations'], atol=1e-12)


def test_annotation_includes_wbo_sum_and_prevents_property_collision(simple):
    mol, result = simple
    target = annotate_mol(mol, result)
    assert target.GetAtomWithIdx(0).GetDoubleProp('xtb_wiberg_sum') == pytest.approx(1.7)
    assert target.GetDoubleProp('xtb_dipole_norm_debye') == pytest.approx(12.708732365)
    assert not mol.GetAtomWithIdx(0).HasProp('xtb_wiberg_sum')
    with pytest.raises(ValueError, match='must differ'):
        annotate_mol(mol, result, copy=False, atom_charge_prop='same', atom_wiberg_prop='same')
    assert not mol.GetAtomWithIdx(0).HasProp('same')
