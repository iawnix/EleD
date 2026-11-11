from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from eled import (
    XTBResult, CalcDensityMatrix, CalcShellPopulations,
    CalcMayerBondDescriptors, CalcFrontierBondDescriptors,
    CalcPiDescriptors,
)
from eled.input import _input_metadata, _xyz_from_mol
from eled.reports import DescriptorSensitivity
from eled import DescriptorCatalog, DescriptorColumns


def result():
    # Two atoms, each with one s and one Cartesian p shell.  Identity S/C
    # keeps this a small analytical test while retaining the real shell maps.
    atom_map = (0, 0, 0, 0, 1, 1, 1, 1)
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
        metadata=_input_metadata(molecule(), -1, 0, 0, _xyz_from_mol(molecule(), -1)),
    )


def molecule():
    mol = Chem.MolFromSmiles('[C]=[C]')
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(((0., 0., 0.), (1.34, 0., 0.))):
        conf.SetAtomPosition(i, xyz)
    conf.Set3D(True); mol.AddConformer(conf)
    return mol


def test_density_and_shell_conservation():
    r = result()
    d = CalcDensityMatrix(r)
    assert d['tr_ps_e'] == pytest.approx(8.)
    assert d['n_err_e'] == pytest.approx(0.)
    shells = CalcShellPopulations(r)
    assert sum(x['low_s_e'] + x['low_p_e']
               for x in shells['atoms']) == pytest.approx(8.)
    assert shells['atoms'][0]['low_ps'] == pytest.approx(3.)
    assert shells['atoms'][0]['low_qrec_e'] == pytest.approx(-4.)


def test_mayer_and_frontier_bonds_and_open_shell_missing_reason():
    mol = molecule()
    r = result()
    mayer = CalcMayerBondDescriptors(mol, r)
    assert mayer['bonds'][0]['wbo'] == 1.
    assert mayer['bonds'][0]['mbo'] == pytest.approx(0.)
    frontier = CalcFrontierBondDescriptors(mol, r)
    assert frontier['bonds']
    open_shell = replace(r, number_of_unpaired_electrons=1)
    assert CalcMayerBondDescriptors(mol, open_shell)['bonds'][0]['mbo'] is None
    assert CalcMayerBondDescriptors(mol, open_shell)['bonds'][0]['mayer_missing_reason']


def test_gfn0_shell_population_is_not_eeq_charge():
    r = replace(result(), method='GFN0-xTB', partial_charges=(.3, -.3))
    rows = CalcShellPopulations(r)['atoms']
    assert rows[0]['low_qrec_e'] is None


def test_density_descriptors_are_invariant_to_ao_permutation_and_phase():
    r = result(); mol = molecule()
    permutation = np.array([2, 0, 3, 1, 6, 4, 7, 5])
    c = np.asarray(r.orbital_coefficients)[permutation]
    s = np.asarray(r.overlap_matrix)[np.ix_(permutation, permutation)]
    ao_atom = np.asarray(r.ao_atom_indices)[permutation]
    ao_shell = np.asarray(r.ao_shell_indices)[permutation]
    changed = replace(r, orbital_coefficients=tuple(map(tuple, c * np.array([1,-1,1,-1,1,-1,1,-1])[:,None])),
                      overlap_matrix=tuple(map(tuple, s)), ao_atom_indices=tuple(ao_atom),
                      ao_shell_indices=tuple(ao_shell))
    a = CalcShellPopulations(r)['atoms']; b = CalcShellPopulations(changed)['atoms']
    for key in ('mul_s_e', 'mul_p_e',
                'low_s_e', 'low_p_e'):
        np.testing.assert_allclose([x[key] for x in a], [x[key] for x in b], atol=1e-12)
    np.testing.assert_allclose(
        [x['mul_ovlp_e'] for x in CalcFrontierBondDescriptors(mol, r)['bonds']],
        [x['mul_ovlp_e'] for x in CalcFrontierBondDescriptors(mol, changed)['bonds']], atol=1e-12)


def test_pi_missing_is_explicit_for_result_without_matching_p_map():
    report = CalcPiDescriptors(molecule(), result())
    assert report['missing']


def test_sensitivity_keeps_missing_comparisons_explicit():
    from test_descriptor_reports import sample
    from eled import CalcDescriptorReport
    mol, r = sample()
    r = replace(r, metadata={**r.metadata, 'engine': {'sha256': 'fixture'}})
    reference = CalcDescriptorReport(mol, r, include_orbitals=False)
    other = CalcDescriptorReport(mol, replace(r, method='GFN1-xTB', total_energy_hartree=-1.5), include_orbitals=False)
    d = DescriptorSensitivity(reference, gfn1=other, gfn2=reference)
    assert d['values']['xtb_energy_ha']['conformer_dispersion'] is None
    assert d['values']['xtb_energy_ha']['gfn12_delta'] == .5
    assert d['values']['xtb_energy_ha']['et300_1000_delta'] is None
    assert d['values']['xtb_energy_ha']['missing']['electronic_temperature']


def test_descriptor_catalog_and_scalar_units_are_explicit():
    catalog = DescriptorCatalog()
    assert {'molecular', 'atom', 'bond', 'density', 'reactivity', 'diagnostic'} <= set(catalog)
    assert catalog['bond']['mbo']['unit'] == 'dimensionless'
    assert catalog['atom']['pi_n_e']['unit'] == 'e'
    assert catalog['atom']['mul_s_e']['unit'] == 'e'
    assert catalog['atom']['low_qrec_e']['unit'] == 'e'
    assert catalog['density']['n_e']['unit'] == 'e'
    assert catalog['reactivity']['ip_v_ev']['unit'] == 'eV'
    assert catalog['molecular']['xtb_dip_*_ebohr']['unit'] == 'e*Bohr'
    assert DescriptorColumns()['xtb_occ_cut'] == 'dimensionless'
    assert all(item['missing'] for group in catalog.values() for item in group.values())
