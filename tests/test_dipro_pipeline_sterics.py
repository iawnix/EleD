"""Tests for 0.16.0: DIPRO electronic couplings, pipeline integration, and steric complements."""
from __future__ import annotations

import json
import math
import os

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (
    AnalyzeLocalGeometry,
    CalcBuriedVolume,
    CalcConeAngle,
    CalcDescriptorReport,
    CalcElectronicCoupling,
    CalcElectronicCouplingFromResults,
    CalcPyramidalization,
    ConformerDescriptorColumns,
    DESCRIPTOR_VERSION,
    DescriptorCatalog,
    DescriptorColumns,
    DescriptorStatus,
    DrawDescriptorMap,
    RunConformerAnalysis,
    XTBCalculator,
    XTBNativeNotFoundError,
    XTBOptions,
    XTBResults,
)
from eled.diagnostics import CheckBackend
from eled.input import _input_metadata, _xyz_from_mol
from eled.options import _validate_calculation_settings
from eled.runtime import _find_xtb_library, _native_basis_build


def _mol_from_coords(symbols: list[str], coords: list[tuple[float, float, float]], bonds: list[tuple[int, int]]):
    rw = Chem.RWMol()
    for sym in symbols:
        atom = Chem.Atom(sym)
        atom.SetNoImplicit(True)
        rw.AddAtom(atom)
    for i, j in bonds:
        rw.AddBond(i, j, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache(strict=False)
    conf = Chem.Conformer(len(symbols))
    conf.Set3D(True)
    for idx, xyz in enumerate(coords):
        conf.SetAtomPosition(idx, tuple(float(v) for v in xyz))
    mol.AddConformer(conf, assignId=True)
    return mol


def _rigid_transform(mol, rotation: np.ndarray, shift: np.ndarray):
    out = Chem.Mol(mol)
    conf = out.GetConformer(0)
    pos = np.asarray(conf.GetPositions(), dtype=float)
    new_pos = pos @ rotation.T + shift[None, :]
    for i, xyz in enumerate(new_pos):
        conf.SetAtomPosition(i, tuple(float(v) for v in xyz))
    return out


def _rotation_matrix() -> np.ndarray:
    a, b, c = 0.63, -0.41, 1.12
    rx = np.array([[1, 0, 0], [0, math.cos(a), -math.sin(a)], [0, math.sin(a), math.cos(a)]])
    ry = np.array([[math.cos(b), 0, math.sin(b)], [0, 1, 0], [-math.sin(b), 0, math.cos(b)]])
    rz = np.array([[math.cos(c), -math.sin(c), 0], [math.sin(c), math.cos(c), 0], [0, 0, 1]])
    return rz @ ry @ rx


def _benzene_dimer(distance_a: float = 3.5):
    mol_a = Chem.AddHs(Chem.MolFromSmiles('c1ccccc1'))
    assert AllChem.EmbedMolecule(mol_a, randomSeed=11) == 0
    AllChem.MMFFOptimizeMolecule(mol_a)
    # Align ring plane cleanly onto xy plane
    conf_a = mol_a.GetConformer(0)
    pos_a = np.asarray(conf_a.GetPositions(), dtype=float)
    pos_a -= pos_a.mean(axis=0)
    _, _, vt = np.linalg.svd(pos_a)
    pos_a = pos_a @ vt.T
    for i, xyz in enumerate(pos_a):
        conf_a.SetAtomPosition(i, tuple(float(v) for v in xyz))

    mol_b = Chem.Mol(mol_a)
    conf_b = mol_b.GetConformer(0)
    for i, xyz in enumerate(pos_a + np.array([0.0, 0.0, distance_a])):
        conf_b.SetAtomPosition(i, tuple(float(v) for v in xyz))

    dimer = Chem.CombineMols(mol_a, mol_b)
    n_a = mol_a.GetNumAtoms()
    frag_a = list(range(n_a))
    frag_b = list(range(n_a, dimer.GetNumAtoms()))
    return dimer, frag_a, frag_b


# ---------------------------------------------------------------------------
# Pure-Python unit tests (run in both .runtime-env and .eled/runtime)
# ---------------------------------------------------------------------------


def test_dipro_exact_two_state_valeev_formula():
    """Exact 2-state dimer matches Valeev (2006) J_eff = (H_AB - 0.5*(eA+eB)*S) / (1 - S^2)."""
    mol = _mol_from_coords(['H', 'H'], [(0.0, 0.0, 0.0), (0.0, 0.0, 3.0)], [])

    # 2 AOs per monomer (1 occupied HOMO + 1 virtual LUMO) -> 4 AOs in dimer
    s_mono = np.eye(2)
    c_a = np.eye(2)
    c_b = np.eye(2)
    e_a = np.array([-6.0, 1.5])
    e_b = np.array([-6.3, 1.8])
    occ_mono = np.array([2.0, 0.0])

    # Build dimer overlap S_AB and dimer Hamiltonian H_AB with known HOMO-HOMO & LUMO-LUMO couplings
    s_hh = 0.14
    s_ll = -0.09
    h_hh = -0.42
    h_ll = 0.31
    h_hl = 0.08
    h_lh = -0.05
    e_ha, e_hb = -6.05, -6.25
    e_la, e_lb = 1.45, 1.85

    s_ab = np.array([
        [1.0, 0.0, s_hh, 0.0],
        [0.0, 1.0, 0.0, s_ll],
        [s_hh, 0.0, 1.0, 0.0],
        [0.0, s_ll, 0.0, 1.0],
    ])
    h_ab = np.array([
        [e_ha, 0.0, h_hh, h_hl],
        [0.0, e_la, h_lh, h_ll],
        [h_hh, h_lh, e_hb, 0.0],
        [h_hl, h_ll, 0.0, e_lb],
    ])
    # Solve generalized eigenproblem H_ab C_ab = S_ab C_ab diag(eps_ab)
    evals_s, evecs_s = np.linalg.eigh(s_ab)
    s_inv_half = (evecs_s * (1.0 / np.sqrt(evals_s))[None, :]) @ evecs_s.T
    h_orth = s_inv_half @ h_ab @ s_inv_half
    e_ab, u_orth = np.linalg.eigh(h_orth)
    c_ab = s_inv_half @ u_orth
    occ_ab = np.array([2.0, 2.0, 0.0, 0.0])

    res_a = XTBResults(
        'GFN2-xTB', 'native', -1.0, 0.0, 2, 2, 0,
        orbitalEnergiesEV=tuple(e_a.tolist()),
        orbitalOccupancies=tuple(occ_mono.tolist()),
        atomicCharges=(0.0,),
        overlapMatrix=tuple(map(tuple, s_mono.tolist())),
        orbitalCoefficients=tuple(map(tuple, c_a.tolist())),
        aoAtomMap=(0, 0),
    )
    res_b = XTBResults(
        'GFN2-xTB', 'native', -1.0, 0.0, 2, 2, 0,
        orbitalEnergiesEV=tuple(e_b.tolist()),
        orbitalOccupancies=tuple(occ_mono.tolist()),
        atomicCharges=(0.0,),
        overlapMatrix=tuple(map(tuple, s_mono.tolist())),
        orbitalCoefficients=tuple(map(tuple, c_b.tolist())),
        aoAtomMap=(0, 0),
    )
    res_ab = XTBResults(
        'GFN2-xTB', 'native', -2.0, 0.0, 4, 4, 0,
        orbitalEnergiesEV=tuple(e_ab.tolist()),
        orbitalOccupancies=tuple(occ_ab.tolist()),
        atomicCharges=(0.0, 0.0),
        overlapMatrix=tuple(map(tuple, s_ab.tolist())),
        orbitalCoefficients=tuple(map(tuple, c_ab.tolist())),
        aoAtomMap=(0, 0, 1, 1),
    )

    out = CalcElectronicCouplingFromResults(mol, res_ab, res_a, res_b, [0], [1])
    expected_j_hole = abs((h_hh - 0.5 * (e_ha + e_hb) * s_hh) / (1.0 - s_hh ** 2))
    expected_j_elec = abs((h_ll - 0.5 * (e_la + e_lb) * s_ll) / (1.0 - s_ll ** 2))
    assert out['j_hole_ev'] == pytest.approx(expected_j_hole, abs=1e-11)
    assert out['j_hole_mev'] == pytest.approx(expected_j_hole * 1000.0, abs=1e-8)
    assert out['j_electron_ev'] == pytest.approx(expected_j_elec, abs=1e-11)
    assert out['j_ct_ab_ev'] == pytest.approx(abs(h_hl), abs=1e-11)
    assert out['j_ct_ba_ev'] == pytest.approx(abs(h_lh), abs=1e-11)
    assert out['hole']['overlap_rms'] == pytest.approx(abs(s_hh), abs=1e-12)

    # MO phase flip invariance
    res_a_flip = XTBResults(
        'GFN2-xTB', 'native', -1.0, 0.0, 2, 2, 0,
        orbitalEnergiesEV=tuple(e_a.tolist()),
        orbitalOccupancies=tuple(occ_mono.tolist()),
        atomicCharges=(0.0,),
        overlapMatrix=tuple(map(tuple, s_mono.tolist())),
        orbitalCoefficients=tuple(map(tuple, (-c_a).tolist())),
        aoAtomMap=(0, 0),
    )
    out_flip = CalcElectronicCouplingFromResults(mol, res_ab, res_a_flip, res_b, [0], [1])
    assert out_flip['j_hole_ev'] == pytest.approx(out['j_hole_ev'], abs=1e-12)
    assert out_flip['j_electron_ev'] == pytest.approx(out['j_electron_ev'], abs=1e-12)
    assert out_flip['hole']['j_signed_ev'] == pytest.approx(-out['hole']['j_signed_ev'], abs=1e-12)

    # Fragment swap invariance
    out_swap = CalcElectronicCouplingFromResults(mol, res_ab, res_b, res_a, [1], [0])
    assert out_swap['j_hole_ev'] == pytest.approx(out['j_hole_ev'], abs=1e-12)
    assert out_swap['j_electron_ev'] == pytest.approx(out['j_electron_ev'], abs=1e-12)
    assert out_swap['j_ct_ab_ev'] == pytest.approx(out['j_ct_ba_ev'], abs=1e-12)
    assert out_swap['j_ct_ba_ev'] == pytest.approx(out['j_ct_ab_ev'], abs=1e-12)

    # Misaligned overlap error guard
    bad_a = XTBResults(
        'GFN2-xTB', 'native', -1.0, 0.0, 2, 2, 0,
        orbitalEnergiesEV=tuple(e_a.tolist()),
        orbitalOccupancies=tuple(occ_mono.tolist()),
        atomicCharges=(0.0,),
        overlapMatrix=((1.0, 0.3), (0.3, 1.0)),
        orbitalCoefficients=tuple(map(tuple, c_a.tolist())),
        aoAtomMap=(0, 0),
    )
    with pytest.raises(ValueError, match='orthonormal|overlap'):
        CalcElectronicCouplingFromResults(mol, res_ab, bad_a, res_b, [0], [1])


def test_pyramidalization_planar_tetrahedral_and_degenerate():
    # Exact trigonal planar BH3
    r = 1.19
    planar = _mol_from_coords(
        ['B', 'H', 'H', 'H'],
        [
            (0.0, 0.0, 0.0),
            (r, 0.0, 0.0),
            (-0.5 * r, 0.5 * math.sqrt(3.0) * r, 0.0),
            (-0.5 * r, -0.5 * math.sqrt(3.0) * r, 0.0),
        ],
        [(0, 1), (0, 2), (0, 3)],
    )
    p_planar = CalcPyramidalization(planar, 0)
    assert p_planar['values']['pyramidalization_p'] == pytest.approx(0.0, abs=1e-12)
    assert p_planar['values']['pyramidalization_angle_deg'] == pytest.approx(0.0, abs=1e-12)
    assert p_planar['values']['poav_angle_deg'] == pytest.approx(90.0, abs=1e-12)
    assert p_planar['values']['angle_sum_deg'] == pytest.approx(360.0, abs=1e-10)
    assert p_planar['values']['pyramid_height_a'] == pytest.approx(0.0, abs=1e-12)
    assert p_planar['missing'] == {}

    # Exact ideal tetrahedral pyramid (bond length 1.0, z = -1/3)
    rho = math.sqrt(8.0 / 9.0)
    z = -1.0 / 3.0
    tet = _mol_from_coords(
        ['N', 'H', 'H', 'H'],
        [
            (0.0, 0.0, 0.0),
            (rho, 0.0, z),
            (-0.5 * rho, 0.5 * math.sqrt(3.0) * rho, z),
            (-0.5 * rho, -0.5 * math.sqrt(3.0) * rho, z),
        ],
        [(0, 1), (0, 2), (0, 3)],
    )
    p_tet = CalcPyramidalization(tet, 0)
    expected_theta = math.degrees(math.asin(1.0 / 3.0))
    assert p_tet['values']['pyramidalization_p'] == pytest.approx(1.0 / 3.0, abs=1e-12)
    assert p_tet['values']['pyramidalization_angle_deg'] == pytest.approx(expected_theta, abs=1e-9)
    assert p_tet['values']['poav_angle_deg'] == pytest.approx(90.0 + expected_theta, abs=1e-9)
    assert p_tet['values']['angle_sum_deg'] == pytest.approx(3.0 * (90.0 + expected_theta), abs=1e-8)
    assert p_tet['values']['pyramid_height_a'] == pytest.approx(1.0 / 3.0, abs=1e-12)

    # Rigid motion and neighbor permutation invariance
    tet_rot = _rigid_transform(tet, _rotation_matrix(), np.array([4.2, -1.7, 2.9]))
    p_rot = CalcPyramidalization(tet_rot, 0, neighborIndices=[3, 1, 2])
    assert p_rot['values'] == pytest.approx(p_tet['values'], abs=1e-10)

    # Degenerate collinear unit tips -> explicit None + missing
    collinear = _mol_from_coords(
        ['N', 'H', 'H', 'H'],
        [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        [(0, 1), (0, 2), (0, 3)],
    )
    p_bad = CalcPyramidalization(collinear, 0)
    assert p_bad['values']['pyramidalization_p'] is None
    assert p_bad['missing']['pyramidalization_p'] == 'degenerate_unit_bond_geometry'
    assert p_bad['values']['pyramid_height_a'] is None
    assert p_bad['missing']['pyramid_height_a'] == 'collinear_neighbors'
    assert json.loads(json.dumps(p_bad, allow_nan=False)) == p_bad


def test_cone_angle_exact_1_2_3_cap_fixed_axis_and_steric_ordering():
    # 1-cap exact: d = 3.0, r = 1.5 -> alpha = 30 deg -> cone = 60 deg
    m1 = _mol_from_coords(['Ni', 'C'], [(0.0, 0.0, 0.0), (0.0, 0.0, 3.0)], [(0, 1)])
    c1 = CalcConeAngle(m1, [1], center=0, radiiAngstrom=[1.0, 1.5])
    assert c1['values']['half_angle_deg'] == pytest.approx(30.0, abs=1e-10)
    assert c1['values']['cone_angle_deg'] == pytest.approx(60.0, abs=1e-10)

    # 2-cap exact: two atoms at +-20 deg from +z, d = 3.0, alpha = 10 deg -> half = 30, cone = 60
    ang20 = math.radians(20.0)
    r10 = 3.0 * math.sin(math.radians(10.0))
    m2 = _mol_from_coords(
        ['Ni', 'H', 'H'],
        [(0.0, 0.0, 0.0), (3.0 * math.sin(ang20), 0.0, 3.0 * math.cos(ang20)),
         (-3.0 * math.sin(ang20), 0.0, 3.0 * math.cos(ang20))],
        [(0, 1), (0, 2)],
    )
    c2 = CalcConeAngle(m2, [1, 2], center=0, radiiAngstrom=[1.0, r10, r10])
    assert c2['values']['half_angle_deg'] == pytest.approx(30.0, abs=1e-9)
    assert c2['values']['cone_angle_deg'] == pytest.approx(60.0, abs=1e-9)

    # 3-cap exact: 3 symmetric atoms on 30-deg cone around +z, alpha = 15 deg -> half = 45, cone = 90
    ang30 = math.radians(30.0)
    r15 = 3.0 * math.sin(math.radians(15.0))
    coords3 = [(0.0, 0.0, 0.0)]
    for k in range(3):
        phi = 2.0 * math.pi * k / 3.0
        coords3.append((
            3.0 * math.sin(ang30) * math.cos(phi),
            3.0 * math.sin(ang30) * math.sin(phi),
            3.0 * math.cos(ang30),
        ))
    m3 = _mol_from_coords(['Ni', 'H', 'H', 'H'], coords3, [(0, 1), (0, 2), (0, 3)])
    c3 = CalcConeAngle(m3, [1, 2, 3], center=0, radiiAngstrom=[1.0, r15, r15, r15])
    assert c3['values']['half_angle_deg'] == pytest.approx(45.0, abs=1e-9)
    assert c3['values']['cone_angle_deg'] == pytest.approx(90.0, abs=1e-9)
    assert np.allclose(c3['cone_axis'], [0.0, 0.0, 1.0], atol=1e-9)
    assert sorted(c3['tangent_atom_indices']) == [1, 2, 3]

    # Rigid rotation + translation invariance
    m3_rot = _rigid_transform(m3, _rotation_matrix(), np.array([-3.0, 5.5, 1.2]))
    c3_rot = CalcConeAngle(m3_rot, [3, 1, 2], center=0, radiiAngstrom=[1.0, r15, r15, r15])
    assert c3_rot['values']['cone_angle_deg'] == pytest.approx(90.0, abs=1e-9)

    # Fixed axis >= free minimum enclosing cone
    c3_fixed = CalcConeAngle(m3, [1, 2, 3], center=0, axisAtom=1, radiiAngstrom=[1.0, r15, r15, r15])
    assert c3_fixed['values']['cone_angle_deg'] > c3['values']['cone_angle_deg']

    # Apex inside sphere -> None + missing
    overlap_mol = _mol_from_coords(['Ni', 'C'], [(0.0, 0.0, 0.0), (0.0, 0.0, 1.0)], [(0, 1)])
    c_ov = CalcConeAngle(overlap_mol, [1], center=0, radiiAngstrom=[1.0, 1.7])
    assert c_ov['values']['cone_angle_deg'] is None
    assert c_ov['missing']['cone_angle_deg'] == 'apex_inside_atom_sphere'
    assert c_ov['overlapping_atom_indices'] == [1]

    # Steric ordering: Ni-PH3 < Ni-PMe3 at 2.28 A M-P bond
    def _build_metal_phosphine(smiles: str):
        lig = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(lig, randomSeed=7) == 0
        AllChem.UFFOptimizeMolecule(lig)
        p_idx = [a.GetIdx() for a in lig.GetAtoms() if a.GetAtomicNum() == 15][0]
        pos = np.asarray(lig.GetConformer(0).GetPositions(), dtype=float)
        p_pos = pos[p_idx]
        nbrs = [n.GetIdx() for n in lig.GetAtomWithIdx(p_idx).GetNeighbors()]
        u = pos[nbrs] - p_pos[None, :]
        u /= np.linalg.norm(u, axis=1, keepdims=True)
        lone_dir = -u.mean(axis=0)
        lone_dir /= np.linalg.norm(lone_dir)
        m_pos = p_pos + 2.28 * lone_dir
        rw = Chem.RWMol(lig)
        ni = Chem.Atom('Ni')
        ni.SetNoImplicit(True)
        m_idx = rw.AddAtom(ni)
        rw.GetAtomWithIdx(p_idx).SetNoImplicit(True)
        rw.AddBond(p_idx, m_idx, Chem.BondType.DATIVE)
        out = rw.GetMol()
        out.RemoveAllConformers()
        out.UpdatePropertyCache(strict=False)
        conf = Chem.Conformer(out.GetNumAtoms())
        conf.Set3D(True)
        for i, xyz in enumerate(pos):
            conf.SetAtomPosition(i, tuple(float(v) for v in xyz))
        conf.SetAtomPosition(m_idx, tuple(float(v) for v in m_pos))
        out.AddConformer(conf, assignId=True)
        lig_atoms = [i for i in range(out.GetNumAtoms()) if i not in (m_idx, p_idx)]
        return out, m_idx, p_idx, lig_atoms

    ni_ph3, m_ph3, p_ph3, a_ph3 = _build_metal_phosphine('P')
    ni_pme3, m_pme3, p_pme3, a_pme3 = _build_metal_phosphine('CP(C)C')
    ang_ph3 = CalcConeAngle(ni_ph3, a_ph3, center=m_ph3, axisAtom=p_ph3)['values']['cone_angle_deg']
    ang_pme3 = CalcConeAngle(ni_pme3, a_pme3, center=m_pme3, axisAtom=p_pme3)['values']['cone_angle_deg']
    assert 60.0 < ang_ph3 < ang_pme3 < 170.0
    free_pme3 = CalcConeAngle(ni_pme3, a_pme3, center=m_pme3)['values']['cone_angle_deg']
    assert free_pme3 <= ang_pme3 + 1e-6


def test_quadrant_buried_volume_and_analyze_local_geometry():
    # Place bulky Cl in Q1 (+x, +y, +z) relative to center=0, zAxisAtom=1, xzPlaneAtom=2
    mol = _mol_from_coords(
        ['Ni', 'H', 'H', 'Cl'],
        [(0.0, 0.0, 0.0), (0.0, 0.0, 5.0), (5.0, 0.0, 0.0), (1.3, 1.3, 1.0)],
        [(0, 1), (0, 2), (0, 3)],
    )
    bv = CalcBuriedVolume(
        mol, [3], center=0, zAxisAtom=1, xzPlaneAtom=2,
        sphereRadiusAngstrom=3.5, spacingAngstrom=0.2,
    )
    vals = bv['values']
    assert vals['vbur_q1_percent'] > vals['vbur_q2_percent']
    assert vals['vbur_q1_percent'] > vals['vbur_q3_percent']
    assert vals['vbur_q1_percent'] > vals['vbur_q4_percent']
    assert vals['vbur_z_pos_percent'] > vals['vbur_z_neg_percent']
    q_mean = 0.25 * sum(vals[f'vbur_q{k}_percent'] for k in range(1, 5))
    assert q_mean == pytest.approx(vals['vbur_percent'], abs=1e-9)
    z_mean = 0.5 * (vals['vbur_z_pos_percent'] + vals['vbur_z_neg_percent'])
    assert z_mean == pytest.approx(vals['vbur_percent'], abs=1e-9)

    # AnalyzeLocalGeometry supports pyramidalization and coneAngle
    combo = AnalyzeLocalGeometry(
        mol,
        selections={'atoms': {'lig': [3]}},
        pyramidalization={'metal': {'atomIdx': 0}},
        coneAngle={'lig': {'center': 0}},
    )
    assert 'pyramidalization:metal' in combo['reports']
    assert 'cone_angle:lig' in combo['reports']
    assert combo['errors'] == {}


def test_catalog_columns_and_visualization_new_fields():
    assert DESCRIPTOR_VERSION == '1.6'
    cols = DescriptorColumns()
    for key in (
        'e_elec_ha', 'e_rep_ha', 'e_es_ha', 'e_disp_ha',
        'mol_pol_bohr3', 'mol_c6_ha_bohr6', 'mol_c8_ha_bohr8',
        'epn_min_au', 'epn_max_au', 'epn_mean_au',
    ):
        assert key in cols
        assert key in ConformerDescriptorColumns()
        assert key in ConformerDescriptorColumns(charge=1)

    cat = DescriptorCatalog()
    assert 'j_hole_ev / j_electron_ev / j_ct_ab_ev / j_ct_ba_ev' in cat['fragment']
    assert 'cone_angle_deg / half_angle_deg' in cat['local_geometry']
    assert 'pyramidalization_p / pyramidalization_angle_deg / poav_angle_deg / angle_sum_deg / pyramid_height_a' in cat['local_geometry']
    status = DescriptorStatus()
    assert 'j_hole_ev / j_electron_ev / j_ct_ab_ev / j_ct_ba_ev' in status['recommended']

    # CalcDescriptorReport with includeExtended=False round-trips analysis_settings
    mol = _mol_from_coords(['O', 'H', 'H'], [(0, 0, 0), (1, 0, 0), (0, 1, 0)], [(0, 1), (0, 2)])
    meta = _input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))
    meta.update(calculation_settings=_validate_calculation_settings(), calculation_status='success')
    r = XTBResults(
        'GFN2-xTB', None, -1.0, 3.0, 3, 4, 0,
        orbitalEnergiesEV=(-10.0, -2.0, 1.0),
        orbitalOccupancies=(2.0, 2.0, 0.0),
        atomicCharges=(-0.6, 0.3, 0.3),
        dipoleAtomicUnits=(0.0, 0.0, 1.0),
        wibergBondOrders=((0.0, 1.0, 1.0), (1.0, 0.0, 0.0), (1.0, 0.0, 0.0)),
        metadata=meta,
    )
    rep_no_ext = CalcDescriptorReport(mol, r, includeExtended=False)
    assert rep_no_ext['analysis_settings']['includeExtended'] is False
    assert rep_no_ext['missing']['energy_decomposition'] == 'not_requested'
    assert rep_no_ext['missing']['mol_pol_bohr3'] == 'not_requested'
    assert CalcDescriptorReport(mol, r, **rep_no_ext['analysis_settings']) == rep_no_ext

    # DrawDescriptorMap handles new atom & bond fields cleanly
    fake_atom_report = {
        'nuclear_esp': {'atoms': [{'atom_index': i, 'epn_au': -1.0 - 0.1 * i} for i in range(3)], 'units': {'epn_au': 'Hartree/e'}},
        'polarizability': {'atoms': [{'atom_index': i, 'alpha_bohr3': 2.0 + i} for i in range(3)], 'units': {'alpha_bohr3': 'Bohr^3'}},
        'bond_overlap': {'bonds': [{'bond_index': i, 'coop_e': 0.4 + 0.1 * i} for i in range(2)], 'units': {'coop_e': 'e'}},
    }
    svg_epn = DrawDescriptorMap(mol, fake_atom_report, atomField='epn_au')
    svg_pol = DrawDescriptorMap(mol, fake_atom_report, atomField='alpha_bohr3')
    svg_coop = DrawDescriptorMap(mol, fake_atom_report, bondField='coop_e')
    for svg in (svg_epn, svg_pol, svg_coop):
        assert '<svg' in svg


# ---------------------------------------------------------------------------
# Integration tests (real libxtb native & compiled backends)
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def library():
    try:
        lib = _find_xtb_library()
    except XTBNativeNotFoundError:
        pytest.skip('libxtb unavailable')
    if not _native_basis_build(lib):
        pytest.skip('extended libxtb basis required')
    return lib


@pytest.mark.integration
def test_dipro_benzene_dimer_degenerate_invariance_and_compiled_parity(library):
    pytest.importorskip('eled.rdXTBTools')
    dimer, frag_a, frag_b = _benzene_dimer(distance_a=3.5)
    opts = XTBOptions(method='gfn2', xtbLibrary=library)

    res_native = CalcElectronicCoupling(dimer, frag_a, frag_b, options=opts, interface='native')
    res_compiled = CalcElectronicCoupling(dimer, frag_a, frag_b, options=opts, interface='compiled')

    assert res_native['hole']['available'] is True
    assert res_native['electron']['available'] is True
    # Benzene has doubly degenerate HOMO (e1g) and LUMO (e2u)
    assert res_native['hole']['n_orb_a'] == 2
    assert res_native['hole']['n_orb_b'] == 2
    assert res_native['electron']['n_orb_a'] == 2
    assert res_native['electron']['n_orb_b'] == 2
    assert 0.05 < res_native['j_hole_ev'] < 2.0
    assert 0.05 < res_native['j_electron_ev'] < 2.0
    assert res_native['hole']['j_signed_ev'] is None
    assert res_native['hole']['missing']['j_signed_ev'] == 'degenerate_subspace_sign_arbitrary'

    # Native vs compiled parity
    for key in ('j_hole_ev', 'j_electron_ev', 'j_ct_ab_ev', 'j_ct_ba_ev'):
        assert res_native[key] == pytest.approx(res_compiled[key], abs=1e-7)

    # Rigid 3D rotation + translation invariance
    dimer_rot = _rigid_transform(dimer, _rotation_matrix(), np.array([2.5, -4.1, 3.3]))
    res_rot = CalcElectronicCoupling(dimer_rot, frag_a, frag_b, options=opts, interface='native')
    for key in ('j_hole_ev', 'j_electron_ev', 'j_ct_ab_ev', 'j_ct_ba_ev'):
        assert res_rot[key] == pytest.approx(res_native[key], abs=1e-5)

    # Monomer swap invariance
    res_swap = CalcElectronicCoupling(dimer, frag_b, frag_a, options=opts, interface='native')
    assert res_swap['j_hole_ev'] == pytest.approx(res_native['j_hole_ev'], abs=1e-8)
    assert res_swap['j_electron_ev'] == pytest.approx(res_native['j_electron_ev'], abs=1e-8)
    assert res_swap['j_ct_ab_ev'] == pytest.approx(res_native['j_ct_ba_ev'], abs=1e-8)
    assert res_swap['j_ct_ba_ev'] == pytest.approx(res_native['j_ct_ab_ev'], abs=1e-8)


@pytest.mark.integration
def test_integrated_report_conformers_viz_and_doctor(library):
    mod = pytest.importorskip('eled.rdXTBTools')
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    AllChem.MMFFOptimizeMolecule(mol)

    native_res = XTBCalculator(method='gfn2', xtbLibrary=library, keepOverlap=True).Calculate(mol)
    ok, comp_res = mod.RunMol(mol, method='gfn2', xtbLibrary=library, keepOverlap=True)
    assert ok

    rep_n = CalcDescriptorReport(mol, native_res, recordId='eth_n')
    rep_c = CalcDescriptorReport(mol, comp_res, recordId='eth_c')

    for sec in (
        'energy_decomposition',
        'polarizability',
        'multipoles',
        'nuclear_esp',
        'orbital_weighted_fukui',
        'atomic_orbital_energies',
        'bond_overlap',
        'localized_orbitals',
        'linear_response',
    ):
        assert rep_n[sec] is not None, f'missing native section {sec}'
        assert rep_c[sec] is not None, f'missing compiled section {sec}'

    for key in (
        'e_elec_ha', 'e_rep_ha', 'e_es_ha', 'e_disp_ha',
        'mol_pol_bohr3', 'mol_c6_ha_bohr6', 'mol_c8_ha_bohr8',
        'epn_min_au', 'epn_max_au', 'epn_mean_au',
    ):
        assert rep_n['values'][key] is not None
        assert rep_n['values'][key] == pytest.approx(rep_c['values'][key], abs=1e-9)

    # Visualization directly from CalcDescriptorReport
    for field in ('epn_au', 'alpha_bohr3', 'f_w_plus', 'occ_e_mean_ev', 'neg_chi_diag_e_per_ev'):
        svg = DrawDescriptorMap(mol, rep_n, atomField=field)
        assert '<svg' in svg
    for field in ('coop_e', 'ecohp_ev'):
        svg = DrawDescriptorMap(mol, rep_n, bondField=field)
        assert '<svg' in svg

    # Conformer analysis picks up extended scalar + local sections
    conf2 = Chem.Conformer(mol.GetConformer(0))
    conf2.SetId(1)
    mol.AddConformer(conf2, assignId=False)
    conf_rep = RunConformerAnalysis(mol, options=XTBOptions(method='gfn2', xtbLibrary=library, keepOverlap=True))
    assert conf_rep['status'] == 'complete'
    assert conf_rep['statistics']['mol_pol_bohr3']['n_valid'] == 2
    assert conf_rep['statistics']['epn_mean_au']['n_valid'] == 2
    assert 'nuclear_esp.epn_au' in conf_rep['atom_statistics']['0']
    assert 'bond_overlap.coop_e' in conf_rep['bond_statistics']['0']

    # Diagnostics reports extended ABI capabilities symmetrically on native & compiled
    for backend in ('native', 'compiled'):
        diag = CheckBackend(backend, 'GFN2-xTB', library=library)
        eng = diag['engine']
        assert eng['properties_export_verified'] is True
        assert eng['esp_api_verified'] is True
        assert eng['realspace_export_verified'] is True
        assert eng['realspace_api_version'] == 2
