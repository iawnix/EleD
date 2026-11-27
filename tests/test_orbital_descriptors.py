"""Analytical, invariance, edge-case, and real-xTB tests for orbital_descriptors."""
from __future__ import annotations

from dataclasses import replace
import json
import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import DescriptorUnavailableError, XTBCalculator, XTBResults
from eled.input import _input_metadata, _xyz_from_mol
from eled.orbital_descriptors import (
    CalcAtomicOrbitalEnergies,
    CalcBondOverlapDescriptors,
    CalcOrbitalWeightedFukui,
)
from eled.runtime import _find_xtb_library, _native_basis_build, _native_gfn0_build


def _diatomic_mol(smiles="[O][C]", distance=1.2):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (distance, 0.0, 0.0))
    mol.AddConformer(conf)
    return mol


def _nonorthogonal_diatomic_result(
    *,
    overlap_ab=0.5,
    energies=(-12.0, -8.0, 1.0, 4.0),
    occupations=(2.0, 2.0, 0.0, 0.0),
    charges=(-0.4, 0.4),
    cm5=(-0.3, 0.3),
    method="GFN2-xTB",
):
    """Build a 2-atom, 4-AO result with exact S-orthonormal MOs and S != I."""
    mol = _diatomic_mol("[O][C]", 1.2)
    s = np.array(
        [
            [1.0, 0.0, overlap_ab, 0.0],
            [0.0, 1.0, 0.0, 0.25 * overlap_ab],
            [overlap_ab, 0.0, 1.0, 0.0],
            [0.0, 0.25 * overlap_ab, 0.0, 1.0],
        ],
        dtype=float,
    )
    eig, vec = np.linalg.eigh(s)
    s_inv_half = (vec / np.sqrt(eig)) @ vec.T
    # Orthogonal matrix in the Lowdin basis so C = S^(-1/2) U satisfies C^T S C = I.
    u = np.array(
        [
            [math.sqrt(0.75), 0.0, math.sqrt(0.25), 0.0],
            [0.0, math.sqrt(0.9), 0.0, math.sqrt(0.1)],
            [math.sqrt(0.25), 0.0, -math.sqrt(0.75), 0.0],
            [0.0, math.sqrt(0.1), 0.0, -math.sqrt(0.9)],
        ],
        dtype=float,
    )
    c = s_inv_half @ u
    meta = _input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))
    result = XTBResults(
        method=method,
        version="test",
        totalEnergyHartree=-5.0,
        homoLumoGapEV=energies[2] - energies[1],
        numOrbitals=4,
        numElectrons=int(round(sum(occupations))),
        numUnpairedElectrons=0,
        orbitalEnergiesEV=tuple(float(v) for v in energies),
        orbitalOccupancies=tuple(float(v) for v in occupations),
        atomicCharges=tuple(float(v) for v in charges),
        cm5Charges=tuple(float(v) for v in cm5) if cm5 is not None else (),
        dipoleAtomicUnits=(0.0, 0.0, 0.0),
        wibergBondOrders=((0.0, 1.0), (1.0, 0.0)),
        orbitalCoefficients=tuple(map(tuple, c)),
        overlapMatrix=tuple(map(tuple, s)),
        aoAtomMap=(0, 0, 1, 1),
        aoShellMap=(0, 1, 2, 3),
        shellAtomMap=(0, 0, 1, 1),
        shellAngularMomenta=(0, 0, 0, 0),
        metadata=meta,
    )
    return mol, result, c, s, u


def _rotate_translate_mol(mol, angle_z=0.73, angle_x=-0.41, shift=(1.5, -2.2, 3.7)):
    rotated = Chem.Mol(mol)
    conf = rotated.GetConformer()
    cz, sz = math.cos(angle_z), math.sin(angle_z)
    cx, sx = math.cos(angle_x), math.sin(angle_x)
    rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    rx = np.array([[1.0, 0.0, 0.0], [0.0, cx, -sx], [0.0, sx, cx]])
    rot = rz @ rx
    shift_vec = np.asarray(shift, dtype=float)
    pos = np.asarray(conf.GetPositions(), dtype=float)
    transformed = pos @ rot.T + shift_vec[None, :]
    for i, xyz in enumerate(transformed):
        conf.SetAtomPosition(i, tuple(map(float, xyz)))
    return rotated


def _embed_mol(smiles: str, seed: int = 1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


@pytest.fixture(scope="module")
def library():
    from eled import XTBNativeNotFoundError

    try:
        path = _find_xtb_library()
    except XTBNativeNotFoundError:
        pytest.skip("libxtb not installed")
    if not _native_basis_build(path):
        pytest.skip("verified native basis extension unavailable")
    return path


def test_orbital_weighted_fukui_analytical_and_sum_rules():
    _, res, _, _, u = _nonorthogonal_diatomic_result()
    delta = 2.0
    out = CalcOrbitalWeightedFukui(res, deltaEV=delta, scheme="lowdin")
    json.dumps(out, allow_nan=False)

    # Lowdin atom populations are exact column sums of U^2:
    # MO 0: (0.75, 0.25), MO 1: (0.90, 0.10), MO 2: (0.25, 0.75), MO 3: (0.10, 0.90)
    mu = 0.5 * (-8.0 + 1.0)
    assert out["mu_ev"] == pytest.approx(mu)
    w_occ = np.exp(-(((mu - np.array([-12.0, -8.0])) / delta) ** 2))
    w_occ /= w_occ.sum()
    w_vir = np.exp(-(((mu - np.array([1.0, 4.0])) / delta) ** 2))
    w_vir /= w_vir.sum()
    p_low = np.array([[0.75, 0.90, 0.25, 0.10], [0.25, 0.10, 0.75, 0.90]])
    expected_minus = p_low[:, :2] @ w_occ
    expected_plus = p_low[:, 2:] @ w_vir

    np.testing.assert_allclose([a["f_w_minus"] for a in out["atoms"]], expected_minus, atol=1e-12)
    np.testing.assert_allclose([a["f_w_plus"] for a in out["atoms"]], expected_plus, atol=1e-12)
    np.testing.assert_allclose(
        [a["f_w_radical"] for a in out["atoms"]], 0.5 * (expected_plus + expected_minus), atol=1e-12
    )
    np.testing.assert_allclose(
        [a["dual_w"] for a in out["atoms"]], expected_plus - expected_minus, atol=1e-12
    )
    assert sum(a["f_w_plus"] for a in out["atoms"]) == pytest.approx(1.0, abs=1e-12)
    assert sum(a["f_w_minus"] for a in out["atoms"]) == pytest.approx(1.0, abs=1e-12)
    assert sum(a["f_w_radical"] for a in out["atoms"]) == pytest.approx(1.0, abs=1e-12)
    assert sum(a["dual_w"] for a in out["atoms"]) == pytest.approx(0.0, abs=1e-12)
    assert out["missing"] == {}

    # Tiny deltaEV must remain numerically stable (no underflow NaN/ZeroDivision) and converge to HOMO/LUMO.
    tiny = CalcOrbitalWeightedFukui(res, deltaEV=1e-3, scheme="lowdin")
    np.testing.assert_allclose([a["f_w_minus"] for a in tiny["atoms"]], [0.90, 0.10], atol=1e-12)
    np.testing.assert_allclose([a["f_w_plus"] for a in tiny["atoms"]], [0.25, 0.75], atol=1e-12)

    # Mulliken scheme also satisfies exact sum rules.
    mul = CalcOrbitalWeightedFukui(res, deltaEV=delta, scheme="mulliken")
    assert sum(a["f_w_plus"] for a in mul["atoms"]) == pytest.approx(1.0, abs=1e-12)
    assert sum(a["f_w_minus"] for a in mul["atoms"]) == pytest.approx(1.0, abs=1e-12)
    assert sum(a["dual_w"] for a in mul["atoms"]) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize(
    "occupations, expected_reason",
    [
        ((0.0, 0.0, 0.0, 0.0), "homo_unavailable"),
        ((2.0, 2.0, 2.0, 2.0), "lumo_unavailable"),
    ],
)
def test_orbital_weighted_fukui_missing_homo_or_lumo(occupations, expected_reason):
    _, res, _, _, _ = _nonorthogonal_diatomic_result(occupations=occupations)
    out = CalcOrbitalWeightedFukui(res)
    json.dumps(out, allow_nan=False)
    assert out["mu_ev"] is None
    assert out["sum_errors"] is None
    assert out["missing"]["mu_ev"] == expected_reason
    assert out["missing"]["sum_errors"] == expected_reason
    if expected_reason == "homo_unavailable":
        assert out["homo_ev"] is None
        assert out["missing"]["homo_ev"] == "homo_unavailable"
    else:
        assert out["lumo_ev"] is None
        assert out["missing"]["lumo_ev"] == "lumo_unavailable"
    for row in out["atoms"]:
        for key in ("f_w_plus", "f_w_minus", "f_w_radical", "dual_w"):
            assert row[key] is None
            assert row["missing"][key] == expected_reason


def test_atomic_orbital_energies_analytical_and_edge_cases():
    _, res, _, _, _ = _nonorthogonal_diatomic_result()
    out = CalcAtomicOrbitalEnergies(res, scheme="lowdin", windowEV=2.0)
    json.dumps(out, allow_nan=False)

    # Occupied MOs: MO 0 (-12 eV, n=2, p=[0.75, 0.25]), MO 1 (-8 eV, n=2, p=[0.90, 0.10])
    # Virtual MOs:  MO 2 (+1 eV, p=[0.25, 0.75]), MO 3 (+4 eV, p=[0.10, 0.90])
    a0, a1 = out["atoms"]
    assert a0["occ_e_mean_ev"] == pytest.approx((0.75 * -12.0 + 0.90 * -8.0) / (0.75 + 0.90))
    assert a1["occ_e_mean_ev"] == pytest.approx((0.25 * -12.0 + 0.10 * -8.0) / (0.25 + 0.10))
    assert a0["vir_e_mean_ev"] == pytest.approx((0.25 * 1.0 + 0.10 * 4.0) / (0.25 + 0.10))
    assert a1["vir_e_mean_ev"] == pytest.approx((0.75 * 1.0 + 0.90 * 4.0) / (0.75 + 0.90))
    # Sign convention: bound valence orbitals give negative occ_e_mean_ev = -ALIE.
    assert a0["occ_e_mean_ev"] < 0.0
    assert "negative" in out["sign_convention"].lower()
    # PDOS windows: windowEV=2.0 includes only HOMO (-8 eV) in [-10, -8] and only LUMO (+1 eV) in [1, 3].
    assert out["homo_window_orb_idx"] == [1]
    assert out["lumo_window_orb_idx"] == [2]
    assert a0["pdos_homo_win_e"] == pytest.approx(2.0 * 0.90)
    assert a1["pdos_homo_win_e"] == pytest.approx(2.0 * 0.10)
    assert a0["pdos_lumo_win"] == pytest.approx(0.25)
    assert a1["pdos_lumo_win"] == pytest.approx(0.75)

    # Window large enough to capture all occupied/virtual orbitals sums to total counts.
    wide = CalcAtomicOrbitalEnergies(res, scheme="lowdin", windowEV=10.0)
    assert sum(a["pdos_homo_win_e"] for a in wide["atoms"]) == pytest.approx(4.0)
    assert sum(a["pdos_lumo_win"] for a in wide["atoms"]) == pytest.approx(2.0)

    # Zero-population atom on occupied subspace must return None + explicit missing reason.
    isolated = replace(
        res,
        orbitalCoefficients=tuple(map(tuple, np.eye(4))),
        overlapMatrix=tuple(map(tuple, np.eye(4))),
        orbitalOccupancies=(2.0, 2.0, 0.0, 0.0),
    )
    iso_out = CalcAtomicOrbitalEnergies(isolated, windowEV=0.0)
    assert iso_out["atoms"][1]["occ_e_mean_ev"] is None
    assert iso_out["atoms"][1]["missing"]["occ_e_mean_ev"] == "zero_occupied_population"
    assert iso_out["atoms"][0]["vir_e_mean_ev"] is None
    assert iso_out["atoms"][0]["missing"]["vir_e_mean_ev"] == "zero_virtual_population"

    # No virtual orbitals -> explicit missing.
    no_vir = CalcAtomicOrbitalEnergies(replace(res, orbitalOccupancies=(2.0, 2.0, 2.0, 2.0)))
    assert no_vir["missing"]["virtual"] == "no_virtual_orbitals"
    assert all(a["vir_e_mean_ev"] is None and a["pdos_lumo_win"] is None for a in no_vir["atoms"])


def test_bond_overlap_descriptors_analytical_trace_and_gfn0():
    mol, res, c, s, _ = _nonorthogonal_diatomic_result()
    out = CalcBondOverlapDescriptors(mol, res)
    json.dumps(out, allow_nan=False)

    assert len(out["bonds"]) == 1
    row = out["bonds"][0]
    assert row["bond_index"] == 0
    assert (row["begin_atom_index"], row["end_atom_index"]) == (0, 1)

    # Direct formula check
    n_occ = np.array([2.0, 2.0])
    e_occ = np.array([-12.0, -8.0])
    c_occ = c[:, :2]
    p_mat = (c_occ * n_occ[None, :]) @ c_occ.T
    w_mat = (c_occ * (n_occ * e_occ)[None, :]) @ c_occ.T
    ia, ib = [0, 1], [2, 3]
    expected_coop = 2.0 * np.sum((p_mat * s)[np.ix_(ia, ib)])
    expected_ecohp = 2.0 * np.sum((w_mat * s)[np.ix_(ia, ib)])
    assert row["coop_e"] == pytest.approx(expected_coop, abs=1e-12)
    assert row["ecohp_ev"] == pytest.approx(expected_ecohp, abs=1e-12)
    # Bonding sign convention: constructive overlap -> COOP > 0, ECOHP < 0.
    assert row["coop_e"] > 0.0
    assert row["ecohp_ev"] < 0.0

    # Exact trace conservation: on-atom + off-diagonal bond overlap = Tr(PS) and Tr(WS).
    on_atom_p = sum(np.sum((p_mat * s)[np.ix_(idx, idx)]) for idx in (ia, ib))
    on_atom_w = sum(np.sum((w_mat * s)[np.ix_(idx, idx)]) for idx in (ia, ib))
    assert on_atom_p + row["coop_e"] == pytest.approx(float(np.sum(n_occ)), abs=1e-12)
    assert on_atom_w + row["ecohp_ev"] == pytest.approx(float(np.sum(n_occ * e_occ)), abs=1e-12)

    # Polarity fields
    assert row["delta_q_e"] == pytest.approx(0.8)
    assert row["signed_delta_q_e"] == pytest.approx(-0.8)
    assert row["q_prod_e2"] == pytest.approx(-0.16)
    assert row["delta_cm5_q_e"] == pytest.approx(0.6)
    assert row["signed_delta_cm5_q_e"] == pytest.approx(-0.6)
    assert row["cm5_q_prod_e2"] == pytest.approx(-0.09)
    assert out["missing"] == {}

    # GFN0 has no CM5 charges -> explicit None + reason, never 0.
    gfn0 = replace(res, method="GFN0-xTB", chargeModel="gfn0_eeq", cm5Charges=())
    gfn0_out = CalcBondOverlapDescriptors(mol, gfn0)
    gfn0_row = gfn0_out["bonds"][0]
    assert gfn0_row["delta_cm5_q_e"] is None
    assert gfn0_row["cm5_q_prod_e2"] is None
    assert gfn0_row["missing"]["delta_cm5_q_e"] == "cm5_unavailable_for_gfn0"
    assert gfn0_out["missing"]["cm5_q_prod_e2"] == "cm5_unavailable_for_gfn0"

    # No occupied orbitals -> explicit None + reason for coop_e and ecohp_ev.
    empty_occ = replace(res, orbitalOccupancies=(0.0, 0.0, 0.0, 0.0), numElectrons=0)
    empty_out = CalcBondOverlapDescriptors(mol, empty_occ)
    assert empty_out["bonds"][0]["coop_e"] is None
    assert empty_out["bonds"][0]["ecohp_ev"] is None
    assert empty_out["bonds"][0]["missing"]["coop_e"] == "no_occupied_orbitals"


def test_ao_permutation_and_mo_phase_invariance():
    mol, res, c, s, _ = _nonorthogonal_diatomic_result()
    perm = np.array([2, 0, 3, 1])
    phases = np.array([1.0, -1.0, -1.0, 1.0])
    c_mod = c[perm, :] * phases[None, :]
    s_mod = s[np.ix_(perm, perm)]
    ao_mod = tuple(np.asarray(res.aoAtomMap)[perm].tolist())
    changed = replace(
        res,
        orbitalCoefficients=tuple(map(tuple, c_mod)),
        overlapMatrix=tuple(map(tuple, s_mod)),
        aoAtomMap=ao_mod,
    )

    for scheme in ("lowdin", "mulliken"):
        f1 = CalcOrbitalWeightedFukui(res, scheme=scheme)
        f2 = CalcOrbitalWeightedFukui(changed, scheme=scheme)
        for key in ("f_w_plus", "f_w_minus", "f_w_radical", "dual_w"):
            np.testing.assert_allclose(
                [a[key] for a in f1["atoms"]], [a[key] for a in f2["atoms"]], atol=1e-12
            )

        e1 = CalcAtomicOrbitalEnergies(res, scheme=scheme)
        e2 = CalcAtomicOrbitalEnergies(changed, scheme=scheme)
        for key in ("occ_e_mean_ev", "vir_e_mean_ev", "occ_e_frontier_ev", "vir_e_frontier_ev", "pdos_homo_win_e", "pdos_lumo_win"):
            np.testing.assert_allclose(
                [a[key] for a in e1["atoms"]], [a[key] for a in e2["atoms"]], atol=1e-12
            )

    b1 = CalcBondOverlapDescriptors(mol, res)
    b2 = CalcBondOverlapDescriptors(mol, changed)
    assert b1["bonds"][0]["coop_e"] == pytest.approx(b2["bonds"][0]["coop_e"], abs=1e-12)
    assert b1["bonds"][0]["ecohp_ev"] == pytest.approx(b2["bonds"][0]["ecohp_ev"], abs=1e-12)


def test_validation_and_unavailable_basis_errors():
    mol, res, _, _, _ = _nonorthogonal_diatomic_result()
    no_basis = replace(res, overlapMatrix=())
    with pytest.raises(DescriptorUnavailableError):
        CalcOrbitalWeightedFukui(no_basis)
    with pytest.raises(DescriptorUnavailableError):
        CalcAtomicOrbitalEnergies(no_basis)
    with pytest.raises(DescriptorUnavailableError):
        CalcBondOverlapDescriptors(mol, no_basis)

    with pytest.raises(ValueError, match="deltaEV"):
        CalcOrbitalWeightedFukui(res, deltaEV=0.0)
    with pytest.raises(ValueError, match="deltaEV"):
        CalcOrbitalWeightedFukui(res, deltaEV=-0.5)
    with pytest.raises(ValueError, match="deltaEV"):
        CalcOrbitalWeightedFukui(res, deltaEV=True)
    with pytest.raises(ValueError, match="deltaEV"):
        CalcOrbitalWeightedFukui(res, deltaEV=np.bool_(True))
    with pytest.raises(ValueError, match="deltaEV"):
        CalcOrbitalWeightedFukui(res, deltaEV="0.5")
    with pytest.raises(ValueError, match="scheme"):
        CalcOrbitalWeightedFukui(res, scheme="hirshfeld")
    with pytest.raises(ValueError, match="windowEV"):
        CalcAtomicOrbitalEnergies(res, windowEV=-1.0)
    with pytest.raises(ValueError, match="windowEV"):
        CalcAtomicOrbitalEnergies(res, windowEV=np.bool_(True))
    with pytest.raises(ValueError, match="windowEV"):
        CalcAtomicOrbitalEnergies(res, windowEV="2.0")

    # Geometry mismatch rejected on bond overlap descriptors
    moved = _rotate_translate_mol(mol, shift=(5.0, 0.0, 0.0))
    with pytest.raises(ValueError, match="coordinates"):
        CalcBondOverlapDescriptors(moved, res)


@pytest.mark.integration
def test_benzene_degenerate_rotation_invariance_real_xtb(library):
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    conf = Chem.Conformer(12)
    conf.Set3D(True)
    for i in range(6):
        for index, radius in ((i, 1.40), (i + 6, 2.48)):
            conf.SetAtomPosition(index, (radius * math.cos(i * math.pi / 3), radius * math.sin(i * math.pi / 3), 0.0))
    mol.AddConformer(conf)

    res = XTBCalculator(method="gfn2", xtbLibrary=library, keepOverlap=True).Calculate(mol)
    energies = np.asarray(res.orbitalEnergiesEV)
    occ = np.asarray(res.orbitalOccupancies)
    occ_idx = np.flatnonzero(occ > 1e-7)
    vir_idx = np.flatnonzero(occ <= 1e-7)
    homo_pair = [int(occ_idx[-2]), int(occ_idx[-1])]
    lumo_pair = [int(vir_idx[0]), int(vir_idx[1])]
    assert abs(energies[homo_pair[0]] - energies[homo_pair[1]]) < 1e-6
    assert abs(energies[lumo_pair[0]] - energies[lumo_pair[1]]) < 1e-6

    c_rot = np.asarray(res.orbitalCoefficients, dtype=float).copy()
    rot2 = np.array([[0.6, -0.8], [0.8, 0.6]])
    c_rot[:, homo_pair] = c_rot[:, homo_pair] @ rot2
    c_rot[:, lumo_pair] = c_rot[:, lumo_pair] @ rot2
    rotated_res = replace(res, orbitalCoefficients=tuple(map(tuple, c_rot)))

    for scheme in ("lowdin", "mulliken"):
        f_base = CalcOrbitalWeightedFukui(res, scheme=scheme)
        f_rot = CalcOrbitalWeightedFukui(rotated_res, scheme=scheme)
        for key in ("f_w_plus", "f_w_minus", "f_w_radical", "dual_w"):
            np.testing.assert_allclose(
                [a[key] for a in f_base["atoms"]], [a[key] for a in f_rot["atoms"]], atol=1e-11
            )
        # All 6 symmetric carbons must have equal orbital-weighted Fukui values.
        c_fplus = [f_base["atoms"][i]["f_w_plus"] for i in range(6)]
        np.testing.assert_allclose(c_fplus, np.mean(c_fplus), atol=1e-6)

        e_base = CalcAtomicOrbitalEnergies(res, scheme=scheme)
        e_rot = CalcAtomicOrbitalEnergies(rotated_res, scheme=scheme)
        for key in ("occ_e_mean_ev", "vir_e_mean_ev", "occ_e_frontier_ev", "vir_e_frontier_ev", "pdos_homo_win_e", "pdos_lumo_win"):
            np.testing.assert_allclose(
                [a[key] for a in e_base["atoms"]], [a[key] for a in e_rot["atoms"]], atol=1e-10
            )

    b_base = CalcBondOverlapDescriptors(mol, res)
    b_rot = CalcBondOverlapDescriptors(mol, rotated_res)
    np.testing.assert_allclose(
        [b["coop_e"] for b in b_base["bonds"]], [b["coop_e"] for b in b_rot["bonds"]], atol=1e-11
    )
    np.testing.assert_allclose(
        [b["ecohp_ev"] for b in b_base["bonds"]], [b["ecohp_ev"] for b in b_rot["bonds"]], atol=1e-10
    )


@pytest.mark.integration
def test_rigid_motion_and_atom_renumbering_invariance_real_xtb(library):
    mol = _embed_mol("CCO", seed=1)
    moved = _rotate_translate_mol(mol)
    calc = XTBCalculator(method="gfn2", xtbLibrary=library, keepOverlap=True)
    r1 = calc.Calculate(mol)
    r2 = calc.Calculate(moved)

    f1, f2 = CalcOrbitalWeightedFukui(r1), CalcOrbitalWeightedFukui(r2)
    for key in ("f_w_plus", "f_w_minus", "f_w_radical", "dual_w"):
        np.testing.assert_allclose([a[key] for a in f1["atoms"]], [a[key] for a in f2["atoms"]], atol=1e-6)

    e1, e2 = CalcAtomicOrbitalEnergies(r1), CalcAtomicOrbitalEnergies(r2)
    for key in ("occ_e_mean_ev", "vir_e_mean_ev", "occ_e_frontier_ev", "vir_e_frontier_ev", "pdos_homo_win_e", "pdos_lumo_win"):
        np.testing.assert_allclose([a[key] for a in e1["atoms"]], [a[key] for a in e2["atoms"]], atol=1e-5)

    b1, b2 = CalcBondOverlapDescriptors(mol, r1), CalcBondOverlapDescriptors(moved, r2)
    for key in ("coop_e", "ecohp_ev", "delta_q_e", "signed_delta_q_e", "q_prod_e2", "delta_cm5_q_e", "cm5_q_prod_e2"):
        np.testing.assert_allclose([b[key] for b in b1["bonds"]], [b[key] for b in b2["bonds"]], atol=1e-5)

    # Atom renumbering: permuted atom rows match reindexed original rows, and bond reversal flips signed_delta_q_e.
    order = list(reversed(range(mol.GetNumAtoms())))
    reindexed = Chem.RenumberAtoms(mol, order)
    r_perm = calc.Calculate(reindexed)
    f_perm = CalcOrbitalWeightedFukui(r_perm)
    np.testing.assert_allclose(
        [f1["atoms"][old]["f_w_minus"] for old in order],
        [a["f_w_minus"] for a in f_perm["atoms"]],
        atol=1e-6,
    )
    b_perm = CalcBondOverlapDescriptors(reindexed, r_perm)
    inv_map = {old: new for new, old in enumerate(order)}
    perm_by_pair = {
        (row["begin_atom_index"], row["end_atom_index"]): row for row in b_perm["bonds"]
    }
    for row in b1["bonds"]:
        na, nb = inv_map[row["begin_atom_index"]], inv_map[row["end_atom_index"]]
        if (na, nb) in perm_by_pair:
            other = perm_by_pair[(na, nb)]
            sign = 1.0
        else:
            other = perm_by_pair[(nb, na)]
            sign = -1.0
        assert other["coop_e"] == pytest.approx(row["coop_e"], abs=1e-5)
        assert other["ecohp_ev"] == pytest.approx(row["ecohp_ev"], abs=1e-5)
        assert other["delta_q_e"] == pytest.approx(row["delta_q_e"], abs=1e-5)
        assert other["signed_delta_q_e"] == pytest.approx(sign * row["signed_delta_q_e"], abs=1e-5)


@pytest.mark.integration
@pytest.mark.parametrize("method", ["GFN0-xTB", "GFN1-xTB", "GFN2-xTB"])
def test_native_and_compiled_parity_real_xtb(library, method):
    if method == "GFN0-xTB" and not _native_gfn0_build(library):
        pytest.skip("native GFN0 requires patched library")
    rdXTBTools = pytest.importorskip("eled.rdXTBTools")
    mol = _embed_mol("CCO", seed=1)
    native = XTBCalculator(method=method, xtbLibrary=library, keepOverlap=True).Calculate(mol)
    ok, compiled = rdXTBTools.RunMol(mol, method=method, xtbLibrary=library, keepOverlap=True)
    assert ok and compiled.overlapRetained

    f_nat = CalcOrbitalWeightedFukui(native)
    f_cmp = CalcOrbitalWeightedFukui(compiled)
    for key in ("f_w_plus", "f_w_minus", "f_w_radical", "dual_w"):
        np.testing.assert_allclose([a[key] for a in f_nat["atoms"]], [a[key] for a in f_cmp["atoms"]], atol=1e-6)

    e_nat = CalcAtomicOrbitalEnergies(native)
    e_cmp = CalcAtomicOrbitalEnergies(compiled)
    for key in ("occ_e_mean_ev", "vir_e_mean_ev", "occ_e_frontier_ev", "vir_e_frontier_ev", "pdos_homo_win_e", "pdos_lumo_win"):
        np.testing.assert_allclose([a[key] for a in e_nat["atoms"]], [a[key] for a in e_cmp["atoms"]], atol=1e-6)

    b_nat = CalcBondOverlapDescriptors(mol, native)
    b_cmp = CalcBondOverlapDescriptors(mol, compiled)
    for key in ("coop_e", "ecohp_ev", "delta_q_e", "signed_delta_q_e", "q_prod_e2"):
        np.testing.assert_allclose([b[key] for b in b_nat["bonds"]], [b[key] for b in b_cmp["bonds"]], atol=1e-5)
    if method == "GFN0-xTB":
        assert all(b["delta_cm5_q_e"] is None for b in b_nat["bonds"])
        assert all(b["delta_cm5_q_e"] is None for b in b_cmp["bonds"])
        assert b_nat["missing"]["delta_cm5_q_e"] == "cm5_unavailable_for_gfn0"
        assert b_cmp["missing"]["delta_cm5_q_e"] == "cm5_unavailable_for_gfn0"
    else:
        for key in ("delta_cm5_q_e", "signed_delta_cm5_q_e", "cm5_q_prod_e2"):
            np.testing.assert_allclose(
                [b[key] for b in b_nat["bonds"]], [b[key] for b in b_cmp["bonds"]], atol=1e-5
            )
