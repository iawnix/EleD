"""Analytical, invariance, withholding, and real-xTB tests for localized_orbitals."""
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
from eled.localized_orbitals import CalcLocalizedOrbitals
from eled.runtime import _find_xtb_library, _native_basis_build


def _diatomic_mol(smiles="[O][O]", distance=1.2):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (distance, 0.0, 0.0))
    mol.AddConformer(conf)
    return mol


def _synthetic_result(
    c: np.ndarray,
    s: np.ndarray,
    ao_atom: tuple[int, ...],
    energies: tuple[float, ...],
    occupations: tuple[float, ...],
    *,
    mol=None,
    unpaired: int | None = 0,
):
    if mol is None:
        mol = _diatomic_mol()
    natoms = mol.GetNumAtoms()
    meta = _input_metadata(mol, -1, 0, 0 if unpaired is None else unpaired, _xyz_from_mol(mol, -1))
    return mol, XTBResults(
        method="GFN2-xTB",
        version="test",
        totalEnergyHartree=-2.0,
        homoLumoGapEV=None,
        numOrbitals=len(energies),
        numElectrons=int(round(sum(occupations))),
        numUnpairedElectrons=unpaired,
        orbitalEnergiesEV=tuple(float(v) for v in energies),
        orbitalOccupancies=tuple(float(v) for v in occupations),
        atomicCharges=(0.0,) * natoms,
        dipoleAtomicUnits=(0.0, 0.0, 0.0),
        wibergBondOrders=tuple(tuple(0.0 for _ in range(natoms)) for _ in range(natoms)),
        orbitalCoefficients=tuple(map(tuple, c)),
        overlapMatrix=tuple(map(tuple, s)),
        aoAtomMap=ao_atom,
        aoShellMap=tuple(range(len(ao_atom))),
        shellAtomMap=ao_atom,
        shellAngularMomenta=(0,) * len(ao_atom),
        metadata=meta,
    )


def _embed_mol(smiles: str, seed: int = 1):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _rotate_translate_mol(mol, angle_z=0.65, angle_y=-0.52, shift=(2.1, -1.4, 0.8)):
    rotated = Chem.Mol(mol)
    conf = rotated.GetConformer()
    cz, sz = math.cos(angle_z), math.sin(angle_z)
    cy, sy = math.cos(angle_y), math.sin(angle_y)
    rz = np.array([[cz, -sz, 0.0], [sz, cz, 0.0], [0.0, 0.0, 1.0]])
    ry = np.array([[cy, 0.0, sy], [0.0, 1.0, 0.0], [-sy, 0.0, cy]])
    rot = rz @ ry
    shift_vec = np.asarray(shift, dtype=float)
    pos = np.asarray(conf.GetPositions(), dtype=float)
    transformed = pos @ rot.T + shift_vec[None, :]
    for i, xyz in enumerate(transformed):
        conf.SetAtomPosition(i, tuple(map(float, xyz)))
    return rotated, rot, shift_vec


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


@pytest.mark.parametrize("scheme", ["lowdin", "mulliken"])
def test_pipek_mezey_delocalized_pair_recovers_exact_localized_lone_pairs(scheme):
    # Two orthogonal AOs (one on atom 0, one on atom 1), mixed 45 deg in canonical MOs.
    # Initial P = 2 * (0.5^2 + 0.5^2) = 1.0; localized P = 2 * (1.0^2 + 0.0^2) = 2.0.
    inv_sqrt2 = 1.0 / math.sqrt(2.0)
    c = np.array([[inv_sqrt2, inv_sqrt2], [inv_sqrt2, -inv_sqrt2]])
    s = np.eye(2)
    mol, res = _synthetic_result(c, s, (0, 1), (-12.0, -8.0), (2.0, 2.0))
    out = CalcLocalizedOrbitals(mol, res, scheme=scheme)
    json.dumps(out, allow_nan=False)

    assert out["converged"] is True
    assert out["iterations"] >= 1
    assert out["initial_functional"] == pytest.approx(1.0, abs=1e-12)
    assert out["localization_functional"] == pytest.approx(2.0, abs=1e-12)
    assert out["n_localized_orbitals"] == 2
    assert sum(o["lmo_energy_ev"] for o in out["orbitals"]) == pytest.approx(-20.0, abs=1e-12)
    for o in out["orbitals"]:
        assert o["kind"] == "lone_pair"
        assert o["n_centers"] == pytest.approx(1.0, abs=1e-10)
        assert o["lmo_energy_ev"] == pytest.approx(-10.0, abs=1e-10)
    assert {o["dominant_atoms"][0] for o in out["orbitals"]} == {0, 1}
    assert out["orbitals"][0]["centroid_a"] == pytest.approx([0.0, 0.0, 0.0], abs=1e-10)
    assert out["orbitals"][1]["centroid_a"] == pytest.approx([1.2, 0.0, 0.0], abs=1e-10)
    for a in out["atoms"]:
        assert a["n_lone_pairs"] == 1
        assert a["n_localized_bonds"] == 0
        assert a["max_lone_pair_energy_ev"] == pytest.approx(-10.0, abs=1e-10)
        assert a["missing"] == {}


def test_single_occupied_bond_orbital_and_nonorthogonal_basis():
    # H2-like 1 occupied bonding MO and 1 virtual antibonding MO with non-zero overlap S_01 = 0.5.
    s = np.array([[1.0, 0.5], [0.5, 1.0]])
    c = np.array(
        [
            [1.0 / math.sqrt(3.0), 1.0],
            [1.0 / math.sqrt(3.0), -1.0],
        ]
    )
    mol = Chem.AddHs(Chem.MolFromSmiles("[H][H]"))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (0.74, 0.0, 0.0))
    mol.AddConformer(conf)
    _, res = _synthetic_result(c, s, (0, 1), (-11.0, 3.0), (2.0, 0.0), mol=mol)
    for scheme in ("lowdin", "mulliken"):
        out = CalcLocalizedOrbitals(mol, res, scheme=scheme)
        json.dumps(out, allow_nan=False)
        assert out["converged"] is True
        assert out["iterations"] == 0
        assert out["n_localized_orbitals"] == 1
        lmo = out["orbitals"][0]
        assert lmo["kind"] == "bond"
        assert lmo["dominant_atoms"] == [0, 1]
        assert lmo["n_centers"] == pytest.approx(2.0, abs=1e-12)
        assert lmo["lmo_energy_ev"] == pytest.approx(-11.0, abs=1e-12)
        assert lmo["centroid_a"] == pytest.approx([0.37, 0.0, 0.0], abs=1e-12)
        for a in out["atoms"]:
            assert a["n_lone_pairs"] == 0
            assert a["n_localized_bonds"] == 1
            assert a["max_lone_pair_energy_ev"] is None
            assert a["missing"]["max_lone_pair_energy_ev"] == "no_lone_pairs_on_atom"


@pytest.mark.parametrize(
    "unpaired, occupations, expected_reason",
    [
        (1, (2.0, 1.0), "open_shell_spin_orbitals_unavailable"),
        (None, (2.0, 0.0), "spin_state_unknown"),
        (0, (1.5, 0.5), "fractional_occupations_unsupported"),
        (0, (0.0, 0.0), "no_occupied_orbitals"),
    ],
)
def test_withholding_for_open_shell_unknown_spin_fractional_and_empty(
    unpaired, occupations, expected_reason
):
    mol, res = _synthetic_result(
        np.eye(2), np.eye(2), (0, 1), (-10.0, -5.0), occupations, unpaired=unpaired
    )
    out = CalcLocalizedOrbitals(mol, res)
    json.dumps(out, allow_nan=False)
    assert out["converged"] is None
    assert out["iterations"] is None
    assert out["localization_functional"] is None
    assert out["initial_functional"] is None
    assert out["n_localized_orbitals"] is None
    assert out["orbitals"] == []
    for key in (
        "localized_orbitals",
        "converged",
        "iterations",
        "localization_functional",
        "initial_functional",
        "n_localized_orbitals",
    ):
        assert out["missing"][key] == expected_reason
    for a in out["atoms"]:
        assert a["n_lone_pairs"] is None
        assert a["n_localized_bonds"] is None
        assert a["max_lone_pair_energy_ev"] is None
        assert a["missing"]["n_lone_pairs"] == expected_reason


def test_ao_permutation_mo_phase_invariance_and_validation():
    s = np.array(
        [
            [1.0, 0.0, 0.4, 0.0],
            [0.0, 1.0, 0.0, 0.2],
            [0.4, 0.0, 1.0, 0.0],
            [0.0, 0.2, 0.0, 1.0],
        ]
    )
    eig, vec = np.linalg.eigh(s)
    s_inv_half = (vec / np.sqrt(eig)) @ vec.T
    theta = 0.6
    u = np.array(
        [
            [math.cos(theta), -math.sin(theta), 0.0, 0.0],
            [math.sin(theta), math.cos(theta), 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ]
    )
    # Mix across atoms so localization is nontrivial
    mix = np.array(
        [
            [0.7, 0.0, 0.7, 0.0],
            [0.0, 0.8, 0.0, 0.6],
            [-0.7, 0.0, 0.7, 0.0],
            [0.0, -0.6, 0.0, 0.8],
        ]
    )
    mix /= np.linalg.norm(mix, axis=0, keepdims=True)
    c = s_inv_half @ (mix @ u)
    mol, res = _synthetic_result(c, s, (0, 0, 1, 1), (-15.0, -10.0, -6.0, 2.0), (2.0, 2.0, 2.0, 0.0))

    perm = np.array([2, 0, 3, 1])
    phases = np.array([1.0, -1.0, 1.0, -1.0])
    changed = replace(
        res,
        orbitalCoefficients=tuple(map(tuple, c[perm, :] * phases[None, :])),
        overlapMatrix=tuple(map(tuple, s[np.ix_(perm, perm)])),
        aoAtomMap=tuple(np.asarray(res.aoAtomMap)[perm].tolist()),
    )
    for scheme in ("lowdin", "mulliken"):
        a = CalcLocalizedOrbitals(mol, res, scheme=scheme)
        b = CalcLocalizedOrbitals(mol, changed, scheme=scheme)
        assert a["localization_functional"] == pytest.approx(b["localization_functional"], abs=1e-10)
        np.testing.assert_allclose(
            [o["lmo_energy_ev"] for o in a["orbitals"]],
            [o["lmo_energy_ev"] for o in b["orbitals"]],
            atol=1e-9,
        )

    with pytest.raises(DescriptorUnavailableError):
        CalcLocalizedOrbitals(mol, replace(res, overlapMatrix=()))
    with pytest.raises(ValueError, match="scheme"):
        CalcLocalizedOrbitals(mol, res, scheme="boys")
    with pytest.raises(ValueError, match="occupationThreshold"):
        CalcLocalizedOrbitals(mol, res, occupationThreshold=0.0)
    with pytest.raises(ValueError, match="occupationThreshold"):
        CalcLocalizedOrbitals(mol, res, occupationThreshold=np.bool_(True))
    with pytest.raises(ValueError, match="occupationThreshold"):
        CalcLocalizedOrbitals(mol, res, occupationThreshold="1.5")
    with pytest.raises(ValueError, match="maxIterations"):
        CalcLocalizedOrbitals(mol, res, maxIterations=0)
    with pytest.raises(ValueError, match="maxIterations"):
        CalcLocalizedOrbitals(mol, res, maxIterations=np.bool_(True))
    with pytest.raises(ValueError, match="tolerance"):
        CalcLocalizedOrbitals(mol, res, tolerance=-1e-8)
    with pytest.raises(ValueError, match="tolerance"):
        CalcLocalizedOrbitals(mol, res, tolerance="1e-8")


@pytest.mark.integration
@pytest.mark.parametrize("scheme", ["lowdin", "mulliken"])
def test_water_ammonia_ethanol_benzene_lmo_chemistry_real_xtb(library, scheme):
    calc = XTBCalculator(method="gfn2", xtbLibrary=library, keepOverlap=True)

    # Water: 2 lone pairs on O + 2 O-H bonds
    water = _embed_mol("O", seed=1)
    r_water = calc.Calculate(water)
    lmo_water = CalcLocalizedOrbitals(water, r_water, scheme=scheme)
    json.dumps(lmo_water, allow_nan=False)
    assert lmo_water["converged"] is True
    assert lmo_water["localization_functional"] >= lmo_water["initial_functional"]
    assert sum(o["kind"] == "lone_pair" for o in lmo_water["orbitals"]) == 2
    assert sum(o["kind"] == "bond" for o in lmo_water["orbitals"]) == 2
    o_row = lmo_water["atoms"][0]
    assert o_row["n_lone_pairs"] == 2
    assert o_row["n_localized_bonds"] == 2
    assert o_row["max_lone_pair_energy_ev"] is not None
    assert lmo_water["atoms"][1]["n_lone_pairs"] == 0
    assert lmo_water["atoms"][1]["n_localized_bonds"] == 1
    assert lmo_water["atoms"][1]["max_lone_pair_energy_ev"] is None
    # Exact energy trace conservation
    occ_e = np.asarray(r_water.orbitalEnergiesEV)[np.asarray(r_water.orbitalOccupancies) >= 1.5]
    assert sum(o["lmo_energy_ev"] for o in lmo_water["orbitals"]) == pytest.approx(float(occ_e.sum()), abs=1e-9)

    # Ammonia: 1 lone pair on N + 3 N-H bonds
    ammonia = _embed_mol("N", seed=1)
    r_nh3 = calc.Calculate(ammonia)
    lmo_nh3 = CalcLocalizedOrbitals(ammonia, r_nh3, scheme=scheme)
    assert lmo_nh3["converged"] is True
    assert sum(o["kind"] == "lone_pair" for o in lmo_nh3["orbitals"]) == 1
    assert sum(o["kind"] == "bond" for o in lmo_nh3["orbitals"]) == 3
    assert lmo_nh3["atoms"][0]["n_lone_pairs"] == 1
    assert lmo_nh3["atoms"][0]["n_localized_bonds"] == 3
    assert lmo_nh3["atoms"][0]["max_lone_pair_energy_ev"] > o_row["max_lone_pair_energy_ev"]

    # Ethanol: 2 lone pairs on O + 8 sigma bonds (C-C, C-O, 5 C-H, O-H)
    ethanol = _embed_mol("CCO", seed=1)
    r_etoh = calc.Calculate(ethanol)
    lmo_etoh = CalcLocalizedOrbitals(ethanol, r_etoh, scheme=scheme)
    assert lmo_etoh["converged"] is True
    assert sum(o["kind"] == "lone_pair" for o in lmo_etoh["orbitals"]) == 2
    assert sum(o["kind"] == "bond" for o in lmo_etoh["orbitals"]) == 8
    o_idx = [a.GetIdx() for a in ethanol.GetAtoms() if a.GetAtomicNum() == 8][0]
    assert lmo_etoh["atoms"][o_idx]["n_lone_pairs"] == 2
    assert lmo_etoh["atoms"][o_idx]["n_localized_bonds"] == 2
    for a in ethanol.GetAtoms():
        if a.GetAtomicNum() == 6:
            assert lmo_etoh["atoms"][a.GetIdx()]["n_lone_pairs"] == 0
            assert lmo_etoh["atoms"][a.GetIdx()]["n_localized_bonds"] == 4

    # Benzene: 12 sigma bonds (6 C-C + 6 C-H) + 3 delocalized aromatic pi orbitals, 0 lone pairs
    benzene = _embed_mol("c1ccccc1", seed=1)
    r_benz = calc.Calculate(benzene)
    lmo_benz = CalcLocalizedOrbitals(benzene, r_benz, scheme=scheme)
    assert lmo_benz["converged"] is True
    assert sum(o["kind"] == "lone_pair" for o in lmo_benz["orbitals"]) == 0
    assert sum(o["kind"] == "bond" for o in lmo_benz["orbitals"]) == 12
    assert sum(o["kind"] == "delocalized" for o in lmo_benz["orbitals"]) == 3
    for a in benzene.GetAtoms():
        row = lmo_benz["atoms"][a.GetIdx()]
        assert row["n_lone_pairs"] == 0
        assert row["n_localized_bonds"] == (3 if a.GetAtomicNum() == 6 else 1)


@pytest.mark.integration
def test_rigid_motion_degenerate_benzene_allyl_and_open_shell_real_xtb(library):
    calc = XTBCalculator(method="gfn2", xtbLibrary=library, keepOverlap=True)

    # Rigid motion invariance & centroid equivariance on ethanol
    ethanol = _embed_mol("CCO", seed=1)
    moved, rot, shift = _rotate_translate_mol(ethanol)
    r1 = calc.Calculate(ethanol)
    r2 = calc.Calculate(moved)
    l1 = CalcLocalizedOrbitals(ethanol, r1, scheme="lowdin")
    l2 = CalcLocalizedOrbitals(moved, r2, scheme="lowdin")
    assert l1["localization_functional"] == pytest.approx(l2["localization_functional"], abs=1e-5)
    np.testing.assert_allclose(
        [o["lmo_energy_ev"] for o in l1["orbitals"]],
        [o["lmo_energy_ev"] for o in l2["orbitals"]],
        atol=1e-4,
    )
    assert [a["n_lone_pairs"] for a in l1["atoms"]] == [a["n_lone_pairs"] for a in l2["atoms"]]
    assert [a["n_localized_bonds"] for a in l1["atoms"]] == [a["n_localized_bonds"] for a in l2["atoms"]]
    # Transform centroids from l1 and verify matching set in l2
    c1_transformed = np.asarray([o["centroid_a"] for o in l1["orbitals"]]) @ rot.T + shift[None, :]
    c2 = np.asarray([o["centroid_a"] for o in l2["orbitals"]])
    dists = np.linalg.norm(c1_transformed[:, None, :] - c2[None, :, :], axis=-1)
    np.testing.assert_allclose(dists.min(axis=1), 0.0, atol=1e-4)

    # Non-convergence with maxIterations=1 is explicitly flagged.
    unconv = CalcLocalizedOrbitals(ethanol, r1, maxIterations=1)
    assert unconv["converged"] is False
    assert unconv["iterations"] == 1
    assert unconv["missing"]["convergence"] == "max_iterations_reached"

    # Ideal D6h benzene: degenerate HOMO rotation preserves functional, energy trace, and counts.
    benz = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    conf = Chem.Conformer(12)
    conf.Set3D(True)
    for i in range(6):
        for index, radius in ((i, 1.40), (i + 6, 2.48)):
            conf.SetAtomPosition(index, (radius * math.cos(i * math.pi / 3), radius * math.sin(i * math.pi / 3), 0.0))
    benz.AddConformer(conf)
    r_benz = calc.Calculate(benz)
    occ_idx = np.flatnonzero(np.asarray(r_benz.orbitalOccupancies) > 1e-7)
    homo_pair = [int(occ_idx[-2]), int(occ_idx[-1])]
    c_rot = np.asarray(r_benz.orbitalCoefficients, dtype=float).copy()
    c_rot[:, homo_pair] = c_rot[:, homo_pair] @ np.array([[0.6, -0.8], [0.8, 0.6]])
    r_benz_rot = replace(r_benz, orbitalCoefficients=tuple(map(tuple, c_rot)))
    lb1 = CalcLocalizedOrbitals(benz, r_benz)
    lb2 = CalcLocalizedOrbitals(benz, r_benz_rot)
    assert lb1["converged"] and lb2["converged"]
    assert lb1["localization_functional"] == pytest.approx(lb2["localization_functional"], abs=1e-4)
    assert sum(o["kind"] == "bond" for o in lb1["orbitals"]) == 12
    assert sum(o["kind"] == "delocalized" for o in lb1["orbitals"]) == 3
    assert [a["n_localized_bonds"] for a in lb1["atoms"]] == [a["n_localized_bonds"] for a in lb2["atoms"]]

    # Allyl cation has a genuine 3-center 2-electron delocalized pi orbital.
    allyl = _embed_mol("C=C[CH2+]", seed=1)
    r_allyl = calc.Calculate(allyl)
    l_allyl = CalcLocalizedOrbitals(allyl, r_allyl)
    assert l_allyl["converged"] is True
    deloc = [o for o in l_allyl["orbitals"] if o["kind"] == "delocalized"]
    assert len(deloc) == 1
    assert set(deloc[0]["dominant_atoms"]) == {0, 1, 2}

    # Open-shell radical [CH3] is explicitly withheld.
    radical = _embed_mol("[CH3]", seed=1)
    r_rad = calc.Calculate(radical)
    l_rad = CalcLocalizedOrbitals(radical, r_rad)
    assert l_rad["converged"] is None
    assert l_rad["missing"]["localized_orbitals"] == "open_shell_spin_orbitals_unavailable"


@pytest.mark.integration
def test_native_and_compiled_parity_real_xtb(library):
    rdXTBTools = pytest.importorskip("eled.rdXTBTools")
    calc = XTBCalculator(method="gfn2", xtbLibrary=library, keepOverlap=True)

    # Native vs compiled parity on water, ammonia, ethanol, and benzene
    for smi in ("O", "N", "CCO", "c1ccccc1"):
        mol = _embed_mol(smi, seed=1)
        native = calc.Calculate(mol)
        ok, compiled = rdXTBTools.RunMol(mol, method="GFN2-xTB", xtbLibrary=library, keepOverlap=True)
        assert ok and compiled.overlapRetained
        for scheme in ("lowdin", "mulliken"):
            ln = CalcLocalizedOrbitals(mol, native, scheme=scheme)
            lc = CalcLocalizedOrbitals(mol, compiled, scheme=scheme)
            assert ln["localization_functional"] == pytest.approx(lc["localization_functional"], abs=1e-6)
            np.testing.assert_allclose(
                [o["lmo_energy_ev"] for o in ln["orbitals"]],
                [o["lmo_energy_ev"] for o in lc["orbitals"]],
                atol=1e-5,
            )
            assert [a["n_lone_pairs"] for a in ln["atoms"]] == [a["n_lone_pairs"] for a in lc["atoms"]]
            assert [a["n_localized_bonds"] for a in ln["atoms"]] == [a["n_localized_bonds"] for a in lc["atoms"]]
