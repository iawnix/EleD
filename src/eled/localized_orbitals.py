"""Pipek-Mezey localized molecular orbitals (LMOs) from retained xTB results.

Every function in this module analyzes an existing ``XTBResults`` (or compiled
``rdXTBTools.XTBResults``) without invoking xTB. Same-calculation overlap,
MO coefficients, and AO-to-atom mapping are required
(``XTBCalculator(...).Calculate(..., keepOverlap=True)``).
"""
from __future__ import annotations

import json
import math
from numbers import Integral, Real

import numpy as np

from .electronic_descriptors import (
    _as_result,
    _has_population_basis,
    _orbital_basis_arrays,
    _orbitals,
)
from .errors import DescriptorUnavailableError
from .mapping import _validate_mapping

_FRACTIONAL_TOL = 1e-2
_LONE_PAIR_POP_THRESHOLD = 0.85
_BOND_TOP2_POP_THRESHOLD = 0.80
_BOND_MIN_SECOND_POP = 0.10
_BOND_MAX_N_CENTERS = 2.5
_DELOCALIZED_ATOM_MIN_POP = 0.10


def _validate_lmo_args(
    scheme: str,
    occupationThreshold: float,
    maxIterations: int,
    tolerance: float,
) -> tuple[str, float, int, float]:
    if scheme not in {"lowdin", "mulliken"}:
        raise ValueError("scheme must be 'lowdin' or 'mulliken'")
    if (
        isinstance(occupationThreshold, (bool, np.bool_))
        or not isinstance(occupationThreshold, Real)
        or not np.isfinite(occupationThreshold)
        or occupationThreshold <= 0.0
        or occupationThreshold > 2.0
    ):
        raise ValueError("occupationThreshold must be finite and in (0, 2]")
    if (
        isinstance(maxIterations, (bool, np.bool_))
        or not isinstance(maxIterations, Integral)
        or maxIterations < 1
    ):
        raise ValueError("maxIterations must be a positive integer")
    if (
        isinstance(tolerance, (bool, np.bool_))
        or not isinstance(tolerance, Real)
        or not np.isfinite(tolerance)
        or tolerance <= 0.0
    ):
        raise ValueError("tolerance must be finite and positive")
    return scheme, float(occupationThreshold), int(maxIterations), float(tolerance)


def _charge_matrices(
    c_occ: np.ndarray,
    s: np.ndarray,
    ao_atom: np.ndarray,
    natoms: int,
    scheme: str,
) -> np.ndarray:
    """Build symmetric atom-charge matrices Q^A of shape (natoms, n_occ, n_occ)."""
    n_occ = c_occ.shape[1]
    q_mats = np.empty((natoms, n_occ, n_occ), dtype=float)
    if scheme == "lowdin":
        eig, vec = np.linalg.eigh(s)
        s_half = (vec * np.sqrt(eig)) @ vec.T
        c_orth = s_half @ c_occ
        for a in range(natoms):
            block = c_orth[ao_atom == a, :]
            q_mats[a] = block.T @ block
    else:
        sc_occ = s @ c_occ
        for a in range(natoms):
            mask = ao_atom == a
            raw = c_occ[mask, :].T @ sc_occ[mask, :]
            q_mats[a] = 0.5 * (raw + raw.T)
    return q_mats


def _pm_functional(q_mats: np.ndarray) -> float:
    diag = np.diagonal(q_mats, axis1=1, axis2=2)
    return float(np.sum(diag * diag))


def _jacobi_pipek_mezey(
    q_mats: np.ndarray,
    *,
    max_iterations: int,
    tolerance: float,
) -> tuple[np.ndarray, np.ndarray, bool, int, float, float]:
    """Maximize P = sum_i sum_A (Q^A_{ii})^2 by 2x2 Jacobi rotations."""
    q = q_mats.copy()
    n_occ = q.shape[1]
    u = np.eye(n_occ, dtype=float)
    initial_p = _pm_functional(q)
    if n_occ <= 1:
        return q, u, True, 0, initial_p, initial_p

    converged = False
    iterations = 0
    current_p = initial_p

    for sweep in range(1, max_iterations + 1):
        iterations = sweep
        max_gain = 0.0
        max_angle = 0.0
        for s_idx in range(n_occ - 1):
            for t_idx in range(s_idx + 1, n_occ):
                diff = q[:, s_idx, s_idx] - q[:, t_idx, t_idx]
                off = q[:, s_idx, t_idx]
                a_st = float(np.sum(off * off - 0.25 * diff * diff))
                b_st = float(np.sum(off * diff))
                r_st = math.hypot(a_st, b_st)
                if r_st <= 1e-15:
                    continue
                gain = (b_st * b_st) / (r_st - a_st) if a_st < 0.0 else (a_st + r_st)
                gamma = 0.25 * math.atan2(b_st, -a_st)
                if gain > max_gain:
                    max_gain = gain
                if abs(gamma) > max_angle:
                    max_angle = abs(gamma)
                if abs(gamma) <= 1e-15:
                    continue
                cos_g = math.cos(gamma)
                sin_g = math.sin(gamma)

                qs = q[:, :, s_idx].copy()
                qt = q[:, :, t_idx].copy()
                q[:, :, s_idx] = cos_g * qs + sin_g * qt
                q[:, :, t_idx] = -sin_g * qs + cos_g * qt

                rs = q[:, s_idx, :].copy()
                rt = q[:, t_idx, :].copy()
                q[:, s_idx, :] = cos_g * rs + sin_g * rt
                q[:, t_idx, :] = -sin_g * rs + cos_g * rt

                us = u[:, s_idx].copy()
                ut = u[:, t_idx].copy()
                u[:, s_idx] = cos_g * us + sin_g * ut
                u[:, t_idx] = -sin_g * us + cos_g * ut

        # Re-symmetrize to prevent floating-point drift across many sweeps.
        q = 0.5 * (q + np.swapaxes(q, 1, 2))
        new_p = _pm_functional(q)
        delta_p = abs(new_p - current_p)
        current_p = new_p
        if max_gain <= tolerance and delta_p <= tolerance and max_angle <= math.sqrt(tolerance):
            converged = True
            break

    return q, u, converged, iterations, initial_p, current_p


def _classify_lmo(atom_pops: np.ndarray) -> tuple[str, list[int], float]:
    """Classify an LMO as lone_pair, bond, or delocalized and compute participation ratio."""
    denom = float(np.sum(atom_pops * atom_pops))
    n_centers = float(1.0 / denom) if denom > 0.0 else float(len(atom_pops))
    # Deterministic order: descending population, ascending atom index on ties.
    order = sorted(range(len(atom_pops)), key=lambda idx: (-round(float(atom_pops[idx]), 12), idx))
    top1 = float(atom_pops[order[0]])
    if top1 >= _LONE_PAIR_POP_THRESHOLD:
        return "lone_pair", [int(order[0])], n_centers
    if (
        len(order) >= 2
        and n_centers <= _BOND_MAX_N_CENTERS
        and (top1 + float(atom_pops[order[1]])) >= _BOND_TOP2_POP_THRESHOLD
        and float(atom_pops[order[1]]) >= _BOND_MIN_SECOND_POP
    ):
        return "bond", [int(order[0]), int(order[1])], n_centers
    dominant = [int(idx) for idx in order if float(atom_pops[idx]) >= _DELOCALIZED_ATOM_MIN_POP]
    if not dominant:
        running = 0.0
        for idx in order:
            dominant.append(int(idx))
            running += float(atom_pops[idx])
            if running >= _BOND_TOP2_POP_THRESHOLD:
                break
    return "delocalized", dominant, n_centers


def _withheld_output(
    natoms: int,
    reason: str,
    *,
    scheme: str,
    occupation_threshold: float,
    max_iterations: int,
    tolerance: float,
    ortho_err: float,
) -> dict:
    atoms = [
        {
            "atom_index": int(a),
            "n_lone_pairs": None,
            "n_localized_bonds": None,
            "max_lone_pair_energy_ev": None,
            "missing": {
                "n_lone_pairs": reason,
                "n_localized_bonds": reason,
                "max_lone_pair_energy_ev": reason,
            },
        }
        for a in range(natoms)
    ]
    output = {
        "converged": None,
        "iterations": None,
        "localization_functional": None,
        "initial_functional": None,
        "n_localized_orbitals": None,
        "orbitals": [],
        "atoms": atoms,
        "scheme": scheme,
        "occupation_threshold": occupation_threshold,
        "max_iterations": max_iterations,
        "tolerance": tolerance,
        "orthonormality_error": ortho_err,
        "analysis_settings": _analysis_settings(scheme, occupation_threshold, max_iterations, tolerance),
        "units": _lmo_units(),
        "invariance": _lmo_invariance(),
        "missing": {
            "localized_orbitals": reason,
            "converged": reason,
            "iterations": reason,
            "localization_functional": reason,
            "initial_functional": reason,
            "n_localized_orbitals": reason,
        },
    }
    json.dumps(output, allow_nan=False)
    return output


def _analysis_settings(
    scheme: str,
    occupation_threshold: float,
    max_iterations: int,
    tolerance: float,
) -> dict:
    return {
        "method": "pipek_mezey_jacobi_2x2",
        "scheme": scheme,
        "occupation_threshold": occupation_threshold,
        "max_iterations": max_iterations,
        "tolerance": tolerance,
        "fractional_occupation_tolerance": _FRACTIONAL_TOL,
        "functional_definition": "P = sum_i sum_A (Q^A_{ii})^2",
        "n_centers_definition": "1 / sum_A p_A^2",
        "classification_rules": {
            "lone_pair": f"max_A p_A >= {_LONE_PAIR_POP_THRESHOLD}",
            "bond": (
                f"not lone_pair, n_centers <= {_BOND_MAX_N_CENTERS}, "
                f"top-2 sum >= {_BOND_TOP2_POP_THRESHOLD}, "
                f"and 2nd atom pop >= {_BOND_MIN_SECOND_POP}"
            ),
            "delocalized": "all remaining LMOs",
        },
        "centroid_definition": "sum_A p_A * R_A (population-weighted atomic position proxy)",
    }


def _lmo_units() -> dict[str, str]:
    return {
        "localization_functional": "dimensionless",
        "initial_functional": "dimensionless",
        "atom_populations": "dimensionless",
        "n_centers": "dimensionless",
        "lmo_energy_ev": "eV",
        "centroid_a": "Angstrom",
        "n_lone_pairs": "count",
        "n_localized_bonds": "count",
        "max_lone_pair_energy_ev": "eV",
    }


def _lmo_invariance() -> str:
    return (
        "localization_functional, orbital energy trace, and symmetry-unique lone-pair/bond "
        "counts are invariant to AO permutation, MO phase, and rigid translation/rotation "
        "(centroid_a transforms equivariantly with molecular coordinates); preserves original "
        "RDKit atom order; individual LMOs are defined up to sign and symmetry-equivalent "
        "degenerate localization solutions (e.g. Kekule resonance structures)"
    )


def CalcLocalizedOrbitals(
    mol,
    result,
    *,
    scheme: str = "lowdin",
    occupationThreshold: float = 1.5,
    maxIterations: int = 200,
    tolerance: float = 1e-8,
):
    """Pipek-Mezey localization of closed-shell occupied MOs via 2x2 Jacobi rotations.

    Maximizes ``P = sum_i sum_A (Q^A_{ii})^2`` on the doubly-occupied subspace,
    using either Lowdin orthogonalized orbitals (``C_orth = S^(1/2) C``) or
    symmetrized Mulliken atom-charge matrices
    ``Q^A = 0.5 * (C^T S_A C + C^T S_A^T C)``.

    Because xTB exposes spin-summed spatial MOs rather than separate alpha/beta
    spin orbitals, localization is explicitly withheld (with ``None`` values and
    explicit ``missing`` reasons) for open-shell calculations, unknown spin
    states, or strongly fractional orbital occupations where unitary occupied
    rotations do not preserve the density matrix.
    """
    scheme, occ_cut, max_iter, tol = _validate_lmo_args(
        scheme, occupationThreshold, maxIterations, tolerance
    )
    result = _as_result(result)
    _validate_mapping(mol, result)
    if not _has_population_basis(result):
        raise DescriptorUnavailableError(
            "localized orbitals require same-calculation overlap, coefficients, and AO-to-atom mapping; "
            "use XTBCalculator(...).Calculate(..., keepOverlap=True)"
        )

    energies, occupations, _, _ = _orbitals(result)
    if np.any(occupations > 2.0 + 1e-8):
        raise ValueError("spin-summed occupations must not exceed two")

    natoms = mol.GetNumAtoms()
    c, s, ao_atom = _orbital_basis_arrays(result, len(energies), natoms)
    eig, _ = np.linalg.eigh(s)
    if eig.min() <= max(1.0, eig.max()) * 1e-10:
        raise ValueError("overlap must be positive definite")
    metric = c.T @ s @ c
    ortho_err = float(np.max(np.abs(metric - np.eye(len(energies)))))
    if ortho_err > 1e-5:
        raise ValueError(f"coefficients are not S-orthonormal: error={ortho_err:g}")

    if result.numUnpairedElectrons is None:
        return _withheld_output(
            natoms,
            "spin_state_unknown",
            scheme=scheme,
            occupation_threshold=occ_cut,
            max_iterations=max_iter,
            tolerance=tol,
            ortho_err=ortho_err,
        )
    if result.numUnpairedElectrons:
        return _withheld_output(
            natoms,
            "open_shell_spin_orbitals_unavailable",
            scheme=scheme,
            occupation_threshold=occ_cut,
            max_iterations=max_iter,
            tolerance=tol,
            ortho_err=ortho_err,
        )

    fractional = np.any((occupations > _FRACTIONAL_TOL) & (occupations < 2.0 - _FRACTIONAL_TOL))
    if fractional:
        return _withheld_output(
            natoms,
            "fractional_occupations_unsupported",
            scheme=scheme,
            occupation_threshold=occ_cut,
            max_iterations=max_iter,
            tolerance=tol,
            ortho_err=ortho_err,
        )

    occ_idx = np.flatnonzero(occupations >= occ_cut)
    if not len(occ_idx):
        return _withheld_output(
            natoms,
            "no_occupied_orbitals",
            scheme=scheme,
            occupation_threshold=occ_cut,
            max_iterations=max_iter,
            tolerance=tol,
            ortho_err=ortho_err,
        )

    c_occ = c[:, occ_idx]
    e_occ = energies[occ_idx]
    q_init = _charge_matrices(c_occ, s, ao_atom, natoms, scheme)
    q_rot, u_mat, converged, iterations, initial_p, final_p = _jacobi_pipek_mezey(
        q_init, max_iterations=max_iter, tolerance=tol
    )

    # Diagonal of rotated Fock/eigenvalue matrix U^T diag(eps_occ) U.
    lmo_energies = np.sum((u_mat * u_mat) * e_occ[:, None], axis=0)
    atom_pops_all = np.diagonal(q_rot, axis1=1, axis2=2)  # shape (natoms, n_occ)
    positions = np.asarray(mol.GetConformer(result.metadata["conformer_id"]).GetPositions(), dtype=float)

    raw_rows = []
    for i in range(len(occ_idx)):
        pops = atom_pops_all[:, i]
        kind, dominant, n_centers = _classify_lmo(pops)
        centroid = pops @ positions
        raw_rows.append({
            "atom_populations": [float(v) for v in pops],
            "dominant_atoms": dominant,
            "n_centers": n_centers,
            "kind": kind,
            "lmo_energy_ev": float(lmo_energies[i]),
            "centroid_a": [float(v) for v in centroid],
        })

    # Sort deterministically by energy ascending, then dominant atoms.
    raw_rows.sort(key=lambda row: (round(row["lmo_energy_ev"], 10), row["dominant_atoms"]))
    orbitals = []
    for idx, row in enumerate(raw_rows):
        orbitals.append({"lmo_index": idx, **row})

    lp_energies: list[list[float]] = [[] for _ in range(natoms)]
    bond_counts = [0 for _ in range(natoms)]
    for row in orbitals:
        if row["kind"] == "lone_pair":
            owner = row["dominant_atoms"][0]
            lp_energies[owner].append(row["lmo_energy_ev"])
        elif row["kind"] == "bond":
            for atom_idx in row["dominant_atoms"]:
                bond_counts[atom_idx] += 1

    atoms = []
    missing_atoms: dict[str, dict[str, str]] = {}
    for a in range(natoms):
        lps = lp_energies[a]
        row_missing: dict[str, str] = {}
        if lps:
            max_lp_e = float(max(lps))
        else:
            max_lp_e = None
            row_missing["max_lone_pair_energy_ev"] = "no_lone_pairs_on_atom"
        if row_missing:
            missing_atoms[str(a)] = dict(row_missing)
        atoms.append({
            "atom_index": int(a),
            "n_lone_pairs": int(len(lps)),
            "n_localized_bonds": int(bond_counts[a]),
            "max_lone_pair_energy_ev": max_lp_e,
            "missing": row_missing,
        })

    missing: dict[str, object] = {}
    if not converged:
        missing["convergence"] = "max_iterations_reached"
    if missing_atoms:
        missing["atoms"] = missing_atoms

    output = {
        "converged": bool(converged),
        "iterations": int(iterations),
        "localization_functional": float(final_p),
        "initial_functional": float(initial_p),
        "n_localized_orbitals": int(len(orbitals)),
        "orbitals": orbitals,
        "atoms": atoms,
        "scheme": scheme,
        "occupation_threshold": occ_cut,
        "max_iterations": max_iter,
        "tolerance": tol,
        "orthonormality_error": ortho_err,
        "analysis_settings": _analysis_settings(scheme, occ_cut, max_iter, tol),
        "units": _lmo_units(),
        "invariance": _lmo_invariance(),
        "missing": missing,
    }
    json.dumps(output, allow_nan=False)
    return output
