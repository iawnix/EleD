from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from rdkit_eth import (
    XTBResult, calc_density_matrix, calc_shell_populations,
    calc_mayer_bond_descriptors, calc_frontier_bond_descriptors,
    calc_pi_descriptors,
)
from rdkit_eth.input import _input_metadata, _xyz_from_mol
from rdkit_eth.reports import descriptor_sensitivity
from rdkit_eth import descriptor_catalog, descriptor_columns


def result():
    # Two atoms, each with one s and one Cartesian p shell.  Identity S/C
    # keeps this a small analytical test while retaining the real shell maps.
    atom_map = (0, 0, 0, 0, 1, 1, 1, 1)
    shell_map = (0, 1, 1, 1, 1, 2, 3, 3, 3)  # corrected below for 8 AOs
    shell_map = (0, 1, 1, 1, 2, 3, 3, 3)
    return XTBResult(
        'GFN2-xTB', 'test', -1.0, 1.0, 8, 8, 0,
        orbital_energies_ev=tuple(float(i) for i in range(8)),
        orbital_occupations=(2., 2., 2., 2., 0., 0., 0., 0.),
        partial_charges=(0., 0.), dipole_e_bohr=(0., 0., 0.),
        wiberg_bond_orders=((0., 1.), (1., 0.)),
        orbital_coefficients=tuple(map(tuple, np.eye(8))),
        overlap_matrix=tuple(map(tuple, np.eye(8))),
        ao_atom_indices=atom_map, ao_shell_indices=shell_map,
        shell_atom_indices=(0, 0, 1, 1), shell_angular_momenta=(0, 1, 0, 1),
        metadata={'atom_count': 2, 'atomic_numbers': [6, 6], 'charge': 0},
    )


def molecule():
    mol = Chem.MolFromSmiles('C=C')
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(((0., 0., 0.), (1.34, 0., 0.))):
        conf.SetAtomPosition(i, xyz)
    conf.Set3D(True); mol.AddConformer(conf)
    return mol


def test_density_and_shell_conservation():
    r = result()
    d = calc_density_matrix(r)
    assert d['trace_PS_e'] == pytest.approx(8.)
    assert d['electron_conservation_error_e'] == pytest.approx(0.)
    shells = calc_shell_populations(r)
    assert sum(x['lowdin_s_electron_count_e'] + x['lowdin_p_electron_count_e']
               for x in shells['atoms']) == pytest.approx(8.)
    assert shells['atoms'][0]['lowdin_p_to_s_ratio'] == pytest.approx(3.)
    assert shells['atoms'][0]['charge_reconstruction_e'] == pytest.approx(-4.)


def test_mayer_and_frontier_bonds_and_open_shell_missing_reason():
    mol = Chem.MolFromSmiles('CC')
    mol.RemoveAllConformers()
    r = result()
    mayer = calc_mayer_bond_descriptors(mol, r)
    assert mayer['bonds'][0]['wiberg_bond_order'] == 1.
    assert mayer['bonds'][0]['mayer_bond_order'] == pytest.approx(0.)
    frontier = calc_frontier_bond_descriptors(mol, r)
    assert frontier['bonds']
    open_shell = replace(r, number_of_unpaired_electrons=1)
    assert calc_mayer_bond_descriptors(mol, open_shell)['bonds'][0]['mayer_bond_order'] is None
    assert calc_mayer_bond_descriptors(mol, open_shell)['bonds'][0]['mayer_missing_reason']


def test_gfn0_shell_population_is_not_eeq_charge():
    r = replace(result(), method='GFN0-xTB', partial_charges=(.3, -.3))
    rows = calc_shell_populations(r)['atoms']
    assert rows[0]['charge_reconstruction_e'] is None


def test_density_descriptors_are_invariant_to_ao_permutation_and_phase():
    r = result(); mol = Chem.MolFromSmiles('CC')
    permutation = np.array([2, 0, 3, 1, 6, 4, 7, 5])
    c = np.asarray(r.orbital_coefficients)[permutation]
    s = np.asarray(r.overlap_matrix)[np.ix_(permutation, permutation)]
    ao_atom = np.asarray(r.ao_atom_indices)[permutation]
    ao_shell = np.asarray(r.ao_shell_indices)[permutation]
    changed = replace(r, orbital_coefficients=tuple(map(tuple, c * np.array([1,-1,1,-1,1,-1,1,-1])[:,None])),
                      overlap_matrix=tuple(map(tuple, s)), ao_atom_indices=tuple(ao_atom),
                      ao_shell_indices=tuple(ao_shell))
    a = calc_shell_populations(r)['atoms']; b = calc_shell_populations(changed)['atoms']
    for key in ('mulliken_s_electron_count_e', 'mulliken_p_electron_count_e',
                'lowdin_s_electron_count_e', 'lowdin_p_electron_count_e'):
        np.testing.assert_allclose([x[key] for x in a], [x[key] for x in b], atol=1e-12)
    np.testing.assert_allclose(
        [x['mulliken_overlap_population'] for x in calc_frontier_bond_descriptors(mol, r)['bonds']],
        [x['mulliken_overlap_population'] for x in calc_frontier_bond_descriptors(mol, changed)['bonds']], atol=1e-12)


def test_pi_missing_is_explicit_for_result_without_matching_p_map():
    report = calc_pi_descriptors(molecule(), result())
    assert report['missing']


def test_sensitivity_keeps_missing_comparisons_explicit():
    reference = {'values': {'x': 1.0}, 'statistics': {'x': {'std': 0.2}}}
    other = {'values': {'x': 1.5}}
    d = descriptor_sensitivity(reference, gfn1=other, gfn2={'values': {'x': 1.0}})
    assert d['values']['x']['conformer_dispersion'] == .2
    assert d['values']['x']['gfn1_vs_gfn2_abs_difference'] == .5
    assert d['values']['x']['electronic_temperature_300_vs_1000_abs_difference'] is None
    assert d['values']['x']['missing']['electronic_temperature']


def test_descriptor_catalog_and_scalar_units_are_explicit():
    catalog = descriptor_catalog()
    assert {'molecular', 'atom', 'bond', 'density', 'reactivity', 'diagnostic'} <= set(catalog)
    assert catalog['bond']['mayer_bond_order']['unit'] == 'dimensionless'
    assert catalog['atom']['pi_electron_count_e']['unit'] == 'e'
    assert catalog['reactivity']['vertical_ip_ev']['unit'] == 'eV'
    assert catalog['molecular']['xtb_dipole_*_e_bohr']['unit'] == 'e*Bohr'
    assert descriptor_columns()['xtb_occupation_threshold'] == 'dimensionless'
    assert all(item['missing'] for group in catalog.values() for item in group.values())
