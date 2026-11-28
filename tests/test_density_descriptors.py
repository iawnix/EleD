from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from eled import (
    XTBResults, CalcDensityMatrix, CalcShellPopulations,
    CalcFrontierBondDescriptors,
    CalcPiDescriptors, CalcMCI, CalcMCIDescriptors,
    DescriptorCatalog, DescriptorColumns, CalcDescriptorReport,
)
from eled.input import _input_metadata, _xyz_from_mol


def result():
    # Two atoms, each with one s and one Cartesian p shell. Identity S/C
    # keeps this a small analytical test while retaining the real shell maps.
    atom_map = (0, 0, 0, 0, 1, 1, 1, 1)
    shell_map = (0, 1, 1, 1, 2, 3, 3, 3)
    return XTBResults(
        'GFN2-xTB', 'test', -1.0, 1.0, 8, 8, 0,
        orbitalEnergiesEV=tuple(float(i) for i in range(8)),
        orbitalOccupancies=(2., 2., 2., 2., 0., 0., 0., 0.),
        atomicCharges=(0., 0.), dipoleAtomicUnits=(0., 0., 0.),
        wibergBondOrders=((0., 1.), (1., 0.)),
        orbitalCoefficients=tuple(map(tuple, np.eye(8))),
        overlapMatrix=tuple(map(tuple, np.eye(8))),
        aoAtomMap=atom_map, aoShellMap=shell_map,
        shellAtomMap=(0, 0, 1, 1), shellAngularMomenta=(0, 1, 0, 1),
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


def test_mci_uses_giambiagi_bultinck_matrix_trace_and_invariance():
    r = result()
    mol = molecule()
    with pytest.raises(ValueError):
        CalcMCI(mol, r, [0, 1, 0])
    # Hückel-like benzene pi system (plus orthogonal sigma core): exact I_ring=32/729, MCI=48/729.
    benz = Chem.MolFromSmiles('[CH0]1=[CH0][CH0]=[CH0][CH0]=[CH0]1')
    conf = Chem.Conformer(6); conf.Set3D(True)
    for i in range(6):
        conf.SetAtomPosition(i, (1.4 * np.cos(i * np.pi / 3), 1.4 * np.sin(i * np.pi / 3), 0.))
    benz.AddConformer(conf)
    h_pi = np.zeros((6, 6))
    for i in range(6):
        h_pi[i, (i + 1) % 6] = h_pi[(i + 1) % 6, i] = -1.0
    _, c_pi = np.linalg.eigh(h_pi)
    # 2 AOs per atom: one localized s AO (occupied, zero inter-atom coupling) + one delocalized pi AO.
    # Total 12 AOs: indices 0..5 are localized s (atoms 0..5), indices 6..11 are pi (atoms 0..5).
    c_full = np.zeros((12, 12))
    c_full[:6, :6] = np.eye(6)
    c_full[6:, 6:] = c_pi
    occ_full = (2.,) * 6 + (2., 2., 2., 0., 0., 0.)
    meta = _input_metadata(benz, -1, 0, 0, _xyz_from_mol(benz, -1))
    benz_res = XTBResults(
        'GFN2-xTB', 'test', -1.0, 1.0, 12, 18, 0,
        orbitalEnergiesEV=tuple(float(i) for i in range(12)),
        orbitalOccupancies=occ_full,
        atomicCharges=(0.,) * 6, dipoleAtomicUnits=(0., 0., 0.),
        wibergBondOrders=tuple(tuple(1. if i != j else 0. for j in range(6)) for i in range(6)),
        orbitalCoefficients=tuple(map(tuple, c_full)),
        overlapMatrix=tuple(map(tuple, np.eye(12))),
        aoAtomMap=(0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5),
        aoShellMap=tuple(range(12)),
        shellAtomMap=(0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 5),
        shellAngularMomenta=(0,) * 12,
        metadata=meta,
    )
    for scheme in ('mulliken', 'lowdin'):
        report = CalcMCI(benz, benz_res, [0, 1, 2, 3, 4, 5], scheme=scheme)
        assert report['i_ring'] == pytest.approx(32.0 / 729.0, abs=1e-12)
        assert report['mci'] == pytest.approx(48.0 / 729.0, abs=1e-12)
        assert 'Giambiagi-Bultinck' in report['definition']
        assert CalcMCIDescriptors(benz, benz_res, {'ring': [0, 1, 2, 3, 4, 5]}, scheme=scheme)['ring'] == report
        # Ring start shift and reversal must preserve both mci and i_ring.
        shifted = CalcMCI(benz, benz_res, [2, 3, 4, 5, 0, 1], scheme=scheme)
        reversed_ring = CalcMCI(benz, benz_res, [0, 5, 4, 3, 2, 1], scheme=scheme)
        assert shifted['mci'] == pytest.approx(report['mci'], abs=1e-12)
        assert shifted['i_ring'] == pytest.approx(report['i_ring'], abs=1e-12)
        assert reversed_ring['mci'] == pytest.approx(report['mci'], abs=1e-12)
        assert reversed_ring['i_ring'] == pytest.approx(report['i_ring'], abs=1e-12)
        # Flipping AO signs and applying an intra-atom orthogonal rotation must leave MCI invariant
        # (this catches element-wise block summation bugs).
        signs = np.array([1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1], dtype=float)
        rot2 = np.array([[0.6, -0.8], [0.8, 0.6]])
        t_mat = np.diag(signs)
        for a in range(6):
            idx = [a, a + 6]
            t_mat[np.ix_(idx, idx)] = rot2 @ t_mat[np.ix_(idx, idx)]
        c_rot = t_mat @ c_full
        s_rot = t_mat @ np.eye(12) @ t_mat.T
        rotated_res = replace(
            benz_res,
            orbitalCoefficients=tuple(map(tuple, c_rot)),
            overlapMatrix=tuple(map(tuple, s_rot)),
        )
        rot_report = CalcMCI(benz, rotated_res, [0, 1, 2, 3, 4, 5], scheme=scheme)
        assert rot_report['mci'] == pytest.approx(report['mci'], abs=1e-12)
        assert rot_report['i_ring'] == pytest.approx(report['i_ring'], abs=1e-12)

    # Open-shell must be explicitly withheld instead of applying closed-shell formula.
    open_shell = replace(benz_res, numUnpairedElectrons=1)
    withheld = CalcMCI(benz, open_shell, [0, 1, 2, 3, 4, 5])
    assert withheld['mci'] is None and withheld['i_ring'] is None
    assert withheld['missing']['mci'] == 'spin_density_unavailable_for_open_shell'
    assert withheld['missing']['i_ring'] == 'spin_density_unavailable_for_open_shell'


def test_descriptor_report_includes_mci_when_basis_is_retained():
    r = result()
    four = Chem.MolFromSmiles('[CH0]1[CH0][CH0][CH0]1')
    conf = Chem.Conformer(4); conf.Set3D(True)
    for i, xy in enumerate(((0., 0.), (1., 0.), (1., 1.), (0., 1.))):
        conf.SetAtomPosition(i, (xy[0], xy[1], 0.))
    four.AddConformer(conf)
    expanded = replace(r, atomicCharges=(0.,) * 4,
                       metadata=_input_metadata(four, -1, 0, 0, _xyz_from_mol(four, -1)),
                       aoAtomMap=(0, 0, 1, 1, 2, 2, 3, 3),
                       aoShellMap=tuple(range(8)),
                       shellAtomMap=(0, 0, 1, 1, 2, 2, 3, 3),
                       shellAngularMomenta=(0,) * 8,
                       wibergBondOrders=tuple(tuple(1. if i != j else 0.
                                                    for j in range(4)) for i in range(4)))
    report = CalcDescriptorReport(four, expanded, mciRings={'ring': [0, 1, 2, 3]})
    assert report['mci']['ring']['mci'] == pytest.approx(0.)
    assert report['missing'].get('mci') is None
    assert report['analysis_settings']['mciScheme'] == 'mulliken'


def test_frontier_bonds_available_from_retained_basis():
    mol = molecule()
    r = result()
    frontier = CalcFrontierBondDescriptors(mol, r)
    assert frontier['bonds']


def test_gfn0_shell_population_is_not_eeq_charge():
    r = replace(result(), method='GFN0-xTB', atomicCharges=(.3, -.3))
    rows = CalcShellPopulations(r)['atoms']
    assert rows[0]['low_qrec_e'] is None


def test_density_descriptors_are_invariant_to_ao_permutation_and_phase():
    r = result(); mol = molecule()
    permutation = np.array([2, 0, 3, 1, 6, 4, 7, 5])
    c = np.asarray(r.orbitalCoefficients)[permutation]
    s = np.asarray(r.overlapMatrix)[np.ix_(permutation, permutation)]
    ao_atom = np.asarray(r.aoAtomMap)[permutation]
    ao_shell = np.asarray(r.aoShellMap)[permutation]
    changed = replace(r, orbitalCoefficients=tuple(map(tuple, c * np.array([1,-1,1,-1,1,-1,1,-1])[:,None])),
                      overlapMatrix=tuple(map(tuple, s)), aoAtomMap=tuple(ao_atom),
                      aoShellMap=tuple(ao_shell))
    a = CalcShellPopulations(r)['atoms']; b = CalcShellPopulations(changed)['atoms']
    for key in ('mul_s_e', 'mul_p_e', 'low_s_e', 'low_p_e'):
        np.testing.assert_allclose([x[key] for x in a], [x[key] for x in b], atol=1e-12)
    np.testing.assert_allclose(
        [x['mul_ovlp_e'] for x in CalcFrontierBondDescriptors(mol, r)['bonds']],
        [x['mul_ovlp_e'] for x in CalcFrontierBondDescriptors(mol, changed)['bonds']], atol=1e-12)


def test_pi_missing_is_explicit_for_result_without_matching_p_map():
    report = CalcPiDescriptors(molecule(), result())
    assert report['missing']


def test_descriptor_catalog_and_scalar_units_are_explicit():
    catalog = DescriptorCatalog()
    assert set(catalog) == {
        'molecular', 'atom', 'fragment', 'bond', 'density', 'reactivity',
        'fod', 'realspace', 'local_geometry', 'esp_surface', 'localized_orbitals',
        'properties', 'solvation', 'thermochemistry',
    }
    assert catalog['fod']['n_fod_e']['unit'] == 'e'
    assert catalog['realspace']['alie_ev / lea_ev']['unit'] == 'eV'
    assert 'mbo' not in catalog['bond']
    assert 'eht_*' not in catalog['molecular']
    assert catalog['bond']['wbo']['unit'] == 'dimensionless'
    assert catalog['atom']['pi_n_e']['unit'] == 'e'
    assert catalog['atom']['mul_s_e']['unit'] == 'e'
    assert catalog['atom']['low_qrec_e']['unit'] == 'e'
    assert catalog['density']['n_e']['unit'] == 'e'
    assert catalog['reactivity']['ip_v_ev']['unit'] == 'eV'
    assert catalog['molecular']['xtb_dip_*_ebohr']['unit'] == 'e*Bohr'
    assert 'xtb_occ_cut' not in DescriptorColumns()
    assert all(item['missing'] for group in catalog.values() for item in group.values())


@pytest.mark.parametrize('field', ['overlapMatrix', 'orbitalCoefficients', 'aoAtomMap'])
@pytest.mark.parametrize('missing', [None, ()])
def test_optional_basis_absence_across_interfaces(field, missing):
    from eled import CalcDescriptorReport, CalcOrbitalPopulations, DescriptorUnavailableError
    r = replace(result(), **{field: missing})
    for read in (CalcDensityMatrix, CalcOrbitalPopulations):
        with pytest.raises(DescriptorUnavailableError):
            read(r)
    report = CalcDescriptorReport(molecule(), r)
    assert report['missing']['populations'] == 'basis_not_retained'
    assert report['values']['xtb_energy_ha'] == r.totalEnergyHartree


@pytest.mark.parametrize('damage', [{'aoAtomMap': (0,) * 8},
                                    {'overlapMatrix': ((1.,),)}])
def test_corrupt_retained_basis_across_interfaces(damage):
    from eled import CalcDescriptorReport, CalcOrbitalPopulations
    r = replace(result(), **damage)
    for read in (CalcDensityMatrix, CalcOrbitalPopulations,
                 lambda data: CalcDescriptorReport(molecule(), data)):
        with pytest.raises(ValueError, match='mapping|square'):
            read(r)
