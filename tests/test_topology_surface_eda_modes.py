"""Tests for EleD 0.15.0: real-space topology/ELF/LOL/NCI, directed surface extrema, fragment EDA, and Cremer-Kraka local modes."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (
    CalcBondDensityTopology,
    CalcDirectedSurfaceDescriptors,
    CalcFragmentInteraction,
    CalcLocalBondForceConstants,
    CalcLocalBondForceConstantsFromHessian,
    CalcNCIDescriptors,
    CalcRealSpaceTopology,
    CalcVibrationalThermochemistry,
    OptimizeXTB,
    XTBOptions,
)
from eled.local_modes import HARTREE_BOHR2_TO_MDYN_A
from eled.runtime import _find_xtb_library, _native_realspace_derivs_build


def _make_mol(smiles: str, seed: int = 17):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _make_diatomic(z1: int, z2: int, dist_a: float, bond_type=Chem.BondType.SINGLE):
    rw = Chem.RWMol()
    for z in (z1, z2):
        a = Chem.Atom(z)
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    rw.AddBond(0, 1, bond_type)
    mol = rw.GetMol()
    mol.UpdatePropertyCache()
    conf = Chem.Conformer(2)
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (dist_a, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def _make_water_dimer():
    """Construct a clean H-bonded water dimer (atoms 0-2 donor, atoms 3-5 acceptor)."""
    rw = Chem.RWMol()
    for z in (8, 1, 1, 8, 1, 1):
        a = Chem.Atom(z)
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    rw.AddBond(0, 1, Chem.BondType.SINGLE)
    rw.AddBond(0, 2, Chem.BondType.SINGLE)
    rw.AddBond(3, 4, Chem.BondType.SINGLE)
    rw.AddBond(3, 5, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache()
    conf = Chem.Conformer(6)
    conf.Set3D(True)
    # Donor water at origin with H1 pointing along +x toward acceptor O3 at x=2.9 A.
    coords = [
        (0.000, 0.000, 0.000),
        (0.960, 0.000, 0.000),
        (-0.240, 0.929, 0.000),
        (2.900, 0.000, 0.000),
        (3.480, 0.000, 0.760),
        (3.480, 0.000, -0.760),
    ]
    for idx, xyz in enumerate(coords):
        conf.SetAtomPosition(idx, xyz)
    mol.AddConformer(conf, assignId=True)
    return mol


def _rotate_translate(mol, angle_rad: float = 0.61, shift=(1.1, -0.9, 0.7)):
    c, s = math.cos(angle_rad), math.sin(angle_rad)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    shift_arr = np.asarray(shift, dtype=float)
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        pos = rot @ np.array([p.x, p.y, p.z], dtype=float) + shift_arr
        conf.SetAtomPosition(i, pos.tolist())
    return out


# ---------------------------------------------------------------------------
# Pure-Python analytical & guard unit tests
# ---------------------------------------------------------------------------


def test_local_bond_force_constants_analytical_diatomic_and_rotation_invariance():
    # Analytical harmonic spring along x with force constant k_au = 0.5 Ha/Bohr^2:
    # V = 0.5 * k * (x1 - x0 - r0)^2 -> H_xx = [[k, -k], [-k, k]].
    # Exact compliance C_00 = 1 / k, exact local force constant k^a = k.
    mol = _make_diatomic(1, 1, 0.74)
    k_au = 0.5
    hess = np.zeros((6, 6), dtype=float)
    hess[0, 0] = k_au
    hess[3, 3] = k_au
    hess[0, 3] = -k_au
    hess[3, 0] = -k_au

    out = CalcLocalBondForceConstantsFromHessian(mol, hess)
    assert out['is_stationary_point'] is True
    assert out['n_bonds'] == 1
    b0 = out['bonds'][0]
    assert b0['k_local_ha_bohr2'] == pytest.approx(k_au, rel=1e-12)
    assert b0['k_local_mdyn_a'] == pytest.approx(k_au * HARTREE_BOHR2_TO_MDYN_A, rel=1e-12)
    assert b0['compliance_bohr2_ha'] == pytest.approx(1.0 / k_au, rel=1e-12)
    assert b0['omega_local_cm1'] > 4000.0
    json.dumps(out, allow_nan=False)

    # Non-positive force constant (inverted saddle Hessian) -> omega_local_cm1 is None with explicit missing!
    out_neg = CalcLocalBondForceConstantsFromHessian(mol, -hess)
    assert out_neg['n_imaginary_modes'] == 1
    assert out_neg['is_stationary_point'] is False
    assert out_neg['bonds'][0]['k_local_ha_bohr2'] == pytest.approx(-k_au, rel=1e-12)
    assert out_neg['bonds'][0]['omega_local_cm1'] is None
    assert out_neg['bonds'][0]['missing']['omega_local_cm1'] == 'nonpositive_local_force_constant'


def test_fragment_eda_and_local_modes_input_guards():
    ethane = _make_mol('CC', seed=11)
    # Covalent bond crosses between C0 and C1 -> rejected unless allowCovalentCut=True.
    with pytest.raises(ValueError, match='covalent bond'):
        CalcFragmentInteraction(ethane, [0, 2, 3, 4], [1, 5, 6, 7])
    with pytest.raises(ValueError, match='disjoint'):
        CalcFragmentInteraction(ethane, [0, 1], [1, 2])
    with pytest.raises(ValueError, match='interface'):
        CalcFragmentInteraction(ethane, [0], [1], interface='cli')

    dimer = _make_water_dimer()
    # Charge mismatch on full partition rejected before running xTB.
    with pytest.raises(ValueError, match='sum to complex charge'):
        CalcFragmentInteraction(dimer, [0, 1, 2], [3, 4, 5], chargeA=1, unpairedA=1, chargeB=0, unpairedB=0)

    # Extra pair validation & NCI grid guards.
    with pytest.raises(ValueError, match='extraPairs'):
        CalcLocalBondForceConstantsFromHessian(dimer, np.zeros((18, 18)), extraPairs=[(0, 99)])
    with pytest.raises(ValueError, match='both fragmentA and fragmentB'):
        CalcNCIDescriptors(dimer, fragmentA=[0, 1, 2])
    with pytest.raises(ValueError, match='maxGridPoints'):
        CalcNCIDescriptors(dimer, gridSpacingAngstrom=0.05, maxGridPoints=100)


# ---------------------------------------------------------------------------
# Integration tests (require verified realspace ABI v2 libxtb)
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def library():
    path = _find_xtb_library(None)
    if not _native_realspace_derivs_build(path):
        pytest.skip('verified realspace ABI v2 libxtb build required')
    return path


@pytest.mark.integration
def test_bond_density_topology_ellipticity_and_elf_on_ethane_ethene_ethyne(library):
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)
    ethane = _make_mol('CC', seed=17)
    ethene = _make_mol('C=C', seed=17)
    ethyne = _make_mol('C#C', seed=17)

    t_single = CalcBondDensityTopology(ethane, options=opts)
    t_double = CalcBondDensityTopology(ethene, options=opts)
    t_triple = CalcBondDensityTopology(ethyne, options=opts)

    # Bond 0 is C-C / C=C / C#C in all three RDKit SMILES.
    cc_single = t_single['bonds'][0]
    cc_double = t_double['bonds'][0]
    cc_triple = t_triple['bonds'][0]

    # Density increases strictly with bond order: single < double < triple.
    assert cc_single['rho_au'] < cc_double['rho_au'] < cc_triple['rho_au']
    # Ellipticity distinguishes pi-asymmetric double bond (eps > 0.15) from cylindrical single/triple (~0).
    assert cc_single['ellipticity'] is not None and cc_single['ellipticity'] < 0.05
    assert cc_double['ellipticity'] is not None and cc_double['ellipticity'] > 0.15
    assert cc_triple['ellipticity'] is not None and cc_triple['ellipticity'] < 0.05
    # Covalent C-C bonds have high ELF (> 0.7).
    assert cc_single['elf'] > 0.7 and cc_double['elf'] > 0.7

    # Rigid rotation invariance on ethene.
    ethene_rot = _rotate_translate(ethene)
    t_double_rot = CalcBondDensityTopology(ethene_rot, options=opts)
    assert t_double_rot['bonds'][0]['rho_au'] == pytest.approx(cc_double['rho_au'], abs=1e-7)
    assert t_double_rot['bonds'][0]['ellipticity'] == pytest.approx(cc_double['ellipticity'], abs=1e-6)
    assert t_double_rot['bonds'][0]['elf'] == pytest.approx(cc_double['elf'], abs=1e-7)

    # Point-level CalcRealSpaceTopology works cleanly on arbitrary points.
    pt_out = CalcRealSpaceTopology(ethene, [(0.0, 0.0, 0.0), (100.0, 0.0, 0.0)], options=opts)
    assert pt_out['points'][0]['elf'] is not None
    assert pt_out['points'][1]['elf'] is None
    assert pt_out['points'][1]['missing']['elf'] == 'density_below_cutoff'
    json.dumps(t_double, allow_nan=False)


@pytest.mark.integration
def test_nci_descriptors_detect_water_dimer_hbond_attraction(library):
    dimer = _make_water_dimer()
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)
    nci = CalcNCIDescriptors(
        dimer,
        gridSpacingAngstrom=0.30,
        paddingAngstrom=0.8,
        fragmentA=[0, 1, 2],
        fragmentB=[3, 4, 5],
        options=opts,
    )
    assert nci['n_nci_points'] > 0
    assert nci['total_nci']['volume_a3'] > 0.0
    # Water dimer H-bond produces attractive negative sign(lambda2)*rho region.
    assert nci['total_nci']['sign_lambda2_rho_min_au'] < -0.005
    assert nci['attractive']['volume_a3'] + nci['vdw']['volume_a3'] > 0.0
    json.dumps(nci, allow_nan=False)


@pytest.mark.integration
def test_directed_surface_extrema_sigma_hole_hbd_hba_and_pi_hole(library):
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)

    # Bromobenzene & iodobenzene: sigma-hole on Br/I and aromatic benzene pi-face, no HBD/HBA!
    brph = _make_mol('Brc1ccccc1', seed=17)
    iph = _make_mol('Ic1ccccc1', seed=17)
    d_brph = CalcDirectedSurfaceDescriptors(brph, options=opts, hbaPointsPerAtom=36)
    d_iph = CalcDirectedSurfaceDescriptors(iph, options=opts, includeHydrogenBonds=False, includePiHoles=False)
    assert len(d_brph['sigma_holes']) == 1
    assert d_brph['v_sigma_max_au'] is not None
    assert d_brph['sigma_holes'][0]['radius_a'] > 1.8
    # Heavier halogen (I) exhibits a stronger (more positive) sigma-hole than Br.
    assert d_iph['v_sigma_max_au'] > d_brph['v_sigma_max_au']
    assert len(d_brph['pi_holes']) == 1
    assert d_brph['v_pi_max_au'] is not None
    assert d_brph['v_hbd_max_au'] is None and d_brph['missing']['v_hbd_max_au'] == 'no_hbd_sites'
    assert d_brph['v_hba_min_au'] is None and d_brph['missing']['v_hba_min_au'] == 'no_hba_sites'

    # Electron-withdrawing CF3Br > PhBr > CH3Br sigma-hole ordering on Br.
    cf3br = _make_mol('FC(F)(F)Br', seed=17)
    ch3br = _make_mol('CBr', seed=17)
    d_cf3br = CalcDirectedSurfaceDescriptors(cf3br, options=opts, includeHydrogenBonds=False, includePiHoles=False)
    d_ch3br = CalcDirectedSurfaceDescriptors(ch3br, options=opts, includeHydrogenBonds=False, includePiHoles=False)
    assert d_cf3br['v_sigma_max_au'] > d_brph['v_sigma_max_au'] > d_ch3br['v_sigma_max_au']

    # Water: positive HBD on O-H axis, negative HBA on O lone pairs, no sigma/pi holes.
    water = _make_mol('O', seed=17)
    d_h2o = CalcDirectedSurfaceDescriptors(water, options=opts, hbaPointsPerAtom=48)
    assert len(d_h2o['hbd_sites']) == 2
    assert d_h2o['v_hbd_max_au'] > 0.02
    assert d_h2o['v_hba_min_au'] < -0.02
    assert d_h2o['v_sigma_max_au'] is None and d_h2o['missing']['v_sigma_max_au'] == 'no_sigma_hole_sites'
    assert d_h2o['v_pi_max_au'] is None and d_h2o['missing']['v_pi_max_au'] == 'no_pi_hole_rings'

    # Rigid rotation invariance of exact axial HBD ray.
    water_rot = _rotate_translate(water)
    d_h2o_rot = CalcDirectedSurfaceDescriptors(water_rot, options=opts, hbaPointsPerAtom=48)
    assert d_h2o_rot['v_hbd_max_au'] == pytest.approx(d_h2o['v_hbd_max_au'], abs=1e-4)
    json.dumps(d_brph, allow_nan=False)


@pytest.mark.integration
def test_fragment_eda_water_dimer_native_compiled_parity_and_charge_transfer(library):
    dimer = _make_water_dimer()
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)

    eda_native = CalcFragmentInteraction(dimer, [0, 1, 2], [3, 4, 5], options=opts, interface='native')
    eda_compiled = CalcFragmentInteraction(dimer, [0, 1, 2], [3, 4, 5], options=opts, interface='compiled')

    # H-bonded water dimer is stabilizing (delta_e_total < 0) with dispersion < 0.
    assert eda_native['delta_e_total_kcalmol'] < -1.0
    assert eda_native['delta_e_disp_kcalmol'] < 0.0
    assert eda_native['wbo_inter'] > 0.01
    # Acceptor (B) donates lone-pair electron density to donor (A) -> ct_a_to_b < 0.
    assert eda_native['ct_lowdin_a_to_b_e'] < 0.0
    assert eda_native['chi_inter_e_per_ev'] is not None and eda_native['chi_inter_e_per_ev'] > 0.0
    assert eda_native['delta_e_xb_ha'] is None and eda_native['missing']['delta_e_xb_ha'] == 'not_applicable_to_gfn2'

    # Native vs compiled exact parity.
    assert eda_native['delta_e_total_ha'] == pytest.approx(eda_compiled['delta_e_total_ha'], abs=1e-10)
    assert eda_native['delta_e_disp_ha'] == pytest.approx(eda_compiled['delta_e_disp_ha'], abs=1e-10)
    assert eda_native['ct_scc_a_to_b_e'] == pytest.approx(eda_compiled['ct_scc_a_to_b_e'], abs=1e-10)
    assert eda_native['coop_inter_e'] == pytest.approx(eda_compiled['coop_inter_e'], abs=1e-10)

    # Swapping A and B flips charge-transfer sign and preserves total energy.
    eda_swap = CalcFragmentInteraction(dimer, [3, 4, 5], [0, 1, 2], options=opts, interface='native')
    assert eda_swap['delta_e_total_ha'] == pytest.approx(eda_native['delta_e_total_ha'], abs=1e-10)
    assert eda_swap['ct_scc_a_to_b_e'] == pytest.approx(-eda_native['ct_scc_a_to_b_e'], abs=1e-10)
    assert eda_swap['ct_lowdin_a_to_b_e'] == pytest.approx(-eda_native['ct_lowdin_a_to_b_e'], abs=1e-10)
    json.dumps(eda_native, allow_nan=False)


@pytest.mark.integration
def test_cremer_kraka_local_force_constants_water_and_thermochemistry_hessian_reuse(library):
    water = _make_mol('O', seed=17)
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)
    opt_water = OptimizeXTB(water, method='GFN2-xTB', xtbLibrary=library).GetMolecule()

    # Reuse Hessian directly from CalcVibrationalThermochemistry(..., includeHessian=True).
    thermo = CalcVibrationalThermochemistry(opt_water, options=opts, includeHessian=True)
    assert 'hessian_ha_bohr2' in thermo
    from_hess = CalcLocalBondForceConstantsFromHessian(opt_water, thermo['hessian_ha_bohr2'])
    direct = CalcLocalBondForceConstants(opt_water, options=opts, extraPairs=[(1, 2)])

    assert direct['is_stationary_point'] is True
    assert direct['n_bonds'] == 2
    assert direct['n_extra_pairs'] == 1
    for b_d, b_h in zip(direct['bonds'], from_hess['bonds']):
        assert b_d['k_local_mdyn_a'] == pytest.approx(b_h['k_local_mdyn_a'], rel=1e-10)
        assert 5.0 < b_d['k_local_mdyn_a'] < 12.0  # physical O-H stretching force constant (~7-9 mDyn/A)
        assert 3000.0 < b_d['omega_local_cm1'] < 4500.0
    # Both equivalent O-H bonds have equal local force constants.
    assert direct['bonds'][0]['k_local_mdyn_a'] == pytest.approx(direct['bonds'][1]['k_local_mdyn_a'], rel=1e-3)
    # Nonbonded H...H extra pair is much softer than covalent O-H.
    assert direct['extra_pairs'][0]['k_local_mdyn_a'] < direct['bonds'][0]['k_local_mdyn_a'] * 0.5

    # Rigid rotation invariance.
    water_rot = _rotate_translate(opt_water)
    rot_res = CalcLocalBondForceConstants(water_rot, options=opts)
    assert rot_res['bonds'][0]['k_local_mdyn_a'] == pytest.approx(direct['bonds'][0]['k_local_mdyn_a'], rel=1e-4)
    json.dumps(direct, allow_nan=False)
