"""Orbital-weighted reactivity, atomic orbital energy/PDOS, and bond overlap descriptors.

Every function in this module analyzes an existing retained ``XTBResults`` (or
compiled ``rdXTBTools.XTBResults``) without invoking xTB. Same-calculation
overlap, MO coefficients, and AO-to-atom mapping are required
(``XTBCalculator(...).Calculate(..., keepOverlap=True)``).
"""
from __future__ import annotations

import json
from numbers import Real

import numpy as np

from .electronic_descriptors import (
    _as_result,
    _charges,
    _cm5_charges,
    _has_population_basis,
    _nonnegative,
    _orbital_basis_arrays,
    _orbitals,
)
from .errors import DescriptorUnavailableError
from .mapping import _validate_mapping
from .populations import CalcOrbitalPopulations

_DENOM_TOL = 1e-12
_WINDOW_TOL_EV = 1e-12


def _validate_positive(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")
    return float(value)


def _validate_nonnegative(value, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be finite and nonnegative")
    _nonnegative(value, name)
    return float(value)


def _validate_scheme(scheme: str) -> str:
    if scheme not in {"lowdin", "mulliken"}:
        raise ValueError("scheme must be 'lowdin' or 'mulliken'")
    return scheme


def _validate_spin_summed_occupations(occupations: np.ndarray) -> None:
    if np.any(occupations > 2.0 + 1e-8):
        raise ValueError("spin-summed occupations must not exceed two")


def _normalized_gaussian_weights(
    energies: np.ndarray,
    indices: np.ndarray,
    center: float,
    width: float,
) -> np.ndarray:
    """Compute normalized Gaussian weights exp(-((center - e) / width)^2) stably."""
    scaled_sq = ((center - energies[indices]) / width) ** 2
    shift = float(np.min(scaled_sq))
    raw = np.exp(-(scaled_sq - shift))
    return raw / float(np.sum(raw))


def CalcOrbitalWeightedFukui(
    result,
    *,
    deltaEV: float = 0.5,
    scheme: str = "lowdin",
):
    """Pino-Rios / Martinez-Araya orbital-weighted condensed Fukui indices.

    Computes nucleophilic (``f_w_plus``), electrophilic (``f_w_minus``),
    radical (``f_w_radical``), and dual (``dual_w``) condensed descriptors from
    one single-point ``XTBResults`` using Gaussian orbital weights centered at
    the frontier midpoint ``mu = (eps_HOMO + eps_LUMO) / 2`` with width
    ``deltaEV > 0``:

    - ``w_i = exp(-((mu - eps_i) / deltaEV)^2)``
    - ``f_w_plus[A]  = sum_{i in vir} w_tilde_i^+ * p_{A,i}``
    - ``f_w_minus[A] = sum_{i in occ} w_tilde_i^- * p_{A,i}``
    - ``f_w_radical[A] = (f_w_plus[A] + f_w_minus[A]) / 2``
    - ``dual_w[A] = f_w_plus[A] - f_w_minus[A]``

    Sum rules hold by construction: ``sum_A f_w_plus = 1``, ``sum_A f_w_minus = 1``,
    ``sum_A f_w_radical = 1``, and ``sum_A dual_w = 0``. When either HOMO or
    LUMO is absent, all condensed values are ``None`` with explicit ``missing``
    reasons.
    """
    result = _as_result(result)
    delta = _validate_positive(deltaEV, "deltaEV")
    scheme = _validate_scheme(scheme)
    energies, occupations, occupied, virtual = _orbitals(result)
    _validate_spin_summed_occupations(occupations)
    pop = CalcOrbitalPopulations(result, scheme=scheme)
    atom_pop = np.asarray(pop["atom_pop"], dtype=float)
    natoms = atom_pop.shape[0]

    homo_ev = float(energies[occupied[-1]]) if len(occupied) else None
    lumo_ev = float(energies[virtual[0]]) if len(virtual) else None

    missing: dict[str, str] = {}
    if homo_ev is None or lumo_ev is None:
        reason = "homo_unavailable" if homo_ev is None else "lumo_unavailable"
        if homo_ev is None:
            missing["homo_ev"] = "homo_unavailable"
        if lumo_ev is None:
            missing["lumo_ev"] = "lumo_unavailable"
        for key in ("mu_orb_ev", "sum_errors", "f_w_plus", "f_w_minus", "f_w_radical", "dual_w"):
            missing[key] = reason
        atoms = [
            {
                "atom_index": int(a),
                "f_w_plus": None,
                "f_w_minus": None,
                "f_w_radical": None,
                "dual_w": None,
                "missing": {
                    "f_w_plus": reason,
                    "f_w_minus": reason,
                    "f_w_radical": reason,
                    "dual_w": reason,
                },
            }
            for a in range(natoms)
        ]
        mu_ev = None
        sum_errors = None
    else:
        mu_ev = 0.5 * (homo_ev + lumo_ev)
        w_minus = _normalized_gaussian_weights(energies, occupied, mu_ev, delta)
        w_plus = _normalized_gaussian_weights(energies, virtual, mu_ev, delta)
        f_minus = atom_pop[:, occupied] @ w_minus
        f_plus = atom_pop[:, virtual] @ w_plus
        f_radical = 0.5 * (f_plus + f_minus)
        dual = f_plus - f_minus
        sum_errors = {
            "f_w_plus": float(np.sum(f_plus) - 1.0),
            "f_w_minus": float(np.sum(f_minus) - 1.0),
            "f_w_radical": float(np.sum(f_radical) - 1.0),
            "dual_w": float(np.sum(dual)),
        }
        atoms = [
            {
                "atom_index": int(a),
                "f_w_plus": float(f_plus[a]),
                "f_w_minus": float(f_minus[a]),
                "f_w_radical": float(f_radical[a]),
                "dual_w": float(dual[a]),
                "missing": {},
            }
            for a in range(natoms)
        ]

    output = {
        "atoms": atoms,
        "mu_orb_ev": mu_ev,
        "homo_ev": homo_ev,
        "lumo_ev": lumo_ev,
        "delta_ev": delta,
        "scheme": scheme,
        "frontier_definition": "aufbau electron count",
        "n_occupied": int(len(occupied)),
        "n_virtual": int(len(virtual)),
        "sum_errors": sum_errors,
        "orthonormality_error": float(pop["orthonormality_error"]),
        "analysis_settings": {
            "method": "orbital_weighted_fukui",
            "reference": "Pino-Rios et al. (2017) / Martinez-Araya (2017)",
            "delta_ev": delta,
            "scheme": scheme,
            "frontier_definition": "aufbau electron count",
            "mu_definition": "(eps_HOMO + eps_LUMO) / 2",
            "weight_definition": "exp(-((mu - eps_i) / delta_ev)^2), normalized over occupied/virtual subspaces",
            "f_w_radical_definition": "(f_w_plus + f_w_minus) / 2",
            "dual_w_definition": "f_w_plus - f_w_minus",
        },
        "units": {
            "mu_orb_ev": "eV",
            "homo_ev": "eV",
            "lumo_ev": "eV",
            "delta_ev": "eV",
            "f_w_plus": "dimensionless",
            "f_w_minus": "dimensionless",
            "f_w_radical": "dimensionless",
            "dual_w": "dimensionless",
        },
        "invariance": (
            "invariant to AO permutation, MO phase, degenerate orbital subspace rotation, "
            "and rigid translation/rotation; preserves original calculation atom order"
        ),
        "missing": missing,
    }
    json.dumps(output, allow_nan=False)
    return output


def CalcAtomicOrbitalEnergies(
    result,
    *,
    scheme: str = "lowdin",
    windowEV: float = 2.0,
):
    """Per-atom condensed ALIE proxy, mean virtual energies, and PDOS window populations.

    Sign convention:
        ``occ_e_mean_ev`` is the occupation-weighted mean occupied MO energy on
        atom ``A`` in eV, ``sum_{i in occ}(n_i * p_{A,i} * eps_i) / sum_{i in occ}(n_i * p_{A,i})``.
        For bound valence orbitals (``eps_i < 0``), ``occ_e_mean_ev`` is
        negative and equals the **negative** of the condensed Average Local
        Ionization Energy (ALIE) proxy (``ALIE_A = -occ_e_mean_ev``). Less
        negative values correspond to more easily ionized sites.

    Frontier-weighted energies (``occ_e_frontier_ev``, ``vir_e_frontier_ev``)
    apply Gaussian frontier weights ``exp(-((eps_frontier - eps_i) / windowEV)^2)``
    anchored at HOMO / LUMO (reducing to the exact degenerate frontier subspace
    when ``windowEV == 0``).

    PDOS windows report occupied electron population ``pdos_homo_win_e`` (in e)
    over ``[eps_HOMO - windowEV, eps_HOMO]`` and virtual state weight
    ``pdos_lumo_win`` (dimensionless) over ``[eps_LUMO, eps_LUMO + windowEV]``.
    """
    result = _as_result(result)
    scheme = _validate_scheme(scheme)
    window = _validate_nonnegative(windowEV, "windowEV")
    energies, occupations, occupied, virtual = _orbitals(result)
    _validate_spin_summed_occupations(occupations)
    pop = CalcOrbitalPopulations(result, scheme=scheme)
    atom_pop = np.asarray(pop["atom_pop"], dtype=float)
    natoms = atom_pop.shape[0]

    homo_ev = float(energies[occupied[-1]]) if len(occupied) else None
    lumo_ev = float(energies[virtual[0]]) if len(virtual) else None

    if len(occupied):
        occ_e = energies[occupied]
        occ_n = occupations[occupied]
        occ_p = atom_pop[:, occupied]
        occ_weights = occ_p * occ_n[None, :]
        occ_denom = np.sum(occ_weights, axis=1)
        occ_num = occ_weights @ occ_e

        if window > 0.0:
            g_occ = np.exp(-(((homo_ev - occ_e) / window) ** 2))
        else:
            g_occ = (np.abs(occ_e - homo_ev) <= _WINDOW_TOL_EV).astype(float)
        occ_f_weights = occ_weights * g_occ[None, :]
        occ_f_denom = np.sum(occ_f_weights, axis=1)
        occ_f_num = occ_f_weights @ occ_e

        homo_win_mask = (occ_e >= homo_ev - window - _WINDOW_TOL_EV) & (occ_e <= homo_ev + _WINDOW_TOL_EV)
        homo_win_idx = [int(occupied[k]) for k in np.flatnonzero(homo_win_mask)]
        pdos_homo = np.sum(occ_weights[:, homo_win_mask], axis=1)
    else:
        occ_denom = occ_num = occ_f_denom = occ_f_num = pdos_homo = None
        homo_win_idx = []

    if len(virtual):
        vir_e = energies[virtual]
        vir_p = atom_pop[:, virtual]
        vir_denom = np.sum(vir_p, axis=1)
        vir_num = vir_p @ vir_e

        if window > 0.0:
            g_vir = np.exp(-(((vir_e - lumo_ev) / window) ** 2))
        else:
            g_vir = (np.abs(vir_e - lumo_ev) <= _WINDOW_TOL_EV).astype(float)
        vir_f_weights = vir_p * g_vir[None, :]
        vir_f_denom = np.sum(vir_f_weights, axis=1)
        vir_f_num = vir_f_weights @ vir_e

        lumo_win_mask = (vir_e >= lumo_ev - _WINDOW_TOL_EV) & (vir_e <= lumo_ev + window + _WINDOW_TOL_EV)
        lumo_win_idx = [int(virtual[k]) for k in np.flatnonzero(lumo_win_mask)]
        pdos_lumo = np.sum(vir_p[:, lumo_win_mask], axis=1)
    else:
        vir_denom = vir_num = vir_f_denom = vir_f_num = pdos_lumo = None
        lumo_win_idx = []

    def _safe_weighted_mean(num: float, denom: float, zero_reason: str):
        if denom <= _DENOM_TOL:
            reason = zero_reason if abs(denom) <= _DENOM_TOL else f"negative_{zero_reason.removeprefix('zero_')}"
            return None, reason
        return float(num / denom), None

    atoms = []
    missing_atoms: dict[str, dict[str, str]] = {}
    for a in range(natoms):
        row_missing: dict[str, str] = {}
        if len(occupied):
            occ_n_e = float(occ_denom[a])
            occ_e_mean, err = _safe_weighted_mean(
                float(occ_num[a]), float(occ_denom[a]), "zero_occupied_population"
            )
            if err is not None:
                row_missing["occ_e_mean_ev"] = err
            occ_e_front, err = _safe_weighted_mean(
                float(occ_f_num[a]), float(occ_f_denom[a]), "zero_frontier_occupied_population"
            )
            if err is not None:
                row_missing["occ_e_frontier_ev"] = err
            pdos_h = float(pdos_homo[a])
        else:
            occ_n_e = occ_e_mean = occ_e_front = pdos_h = None
            for key in ("occ_n_e", "occ_e_mean_ev", "occ_e_frontier_ev", "pdos_homo_win_e"):
                row_missing[key] = "no_occupied_orbitals"

        if len(virtual):
            vir_n = float(vir_denom[a])
            vir_e_mean, err = _safe_weighted_mean(
                float(vir_num[a]), float(vir_denom[a]), "zero_virtual_population"
            )
            if err is not None:
                row_missing["vir_e_mean_ev"] = err
            vir_e_front, err = _safe_weighted_mean(
                float(vir_f_num[a]), float(vir_f_denom[a]), "zero_frontier_virtual_population"
            )
            if err is not None:
                row_missing["vir_e_frontier_ev"] = err
            pdos_l = float(pdos_lumo[a])
        else:
            vir_n = vir_e_mean = vir_e_front = pdos_l = None
            for key in ("vir_n", "vir_e_mean_ev", "vir_e_frontier_ev", "pdos_lumo_win"):
                row_missing[key] = "no_virtual_orbitals"

        row = {
            "atom_index": int(a),
            "occ_n_e": occ_n_e,
            "vir_n": vir_n,
            "occ_e_mean_ev": occ_e_mean,
            "vir_e_mean_ev": vir_e_mean,
            "occ_e_frontier_ev": occ_e_front,
            "vir_e_frontier_ev": vir_e_front,
            "pdos_homo_win_e": pdos_h,
            "pdos_lumo_win": pdos_l,
            "missing": row_missing,
        }
        if row_missing:
            missing_atoms[str(a)] = dict(row_missing)
        atoms.append(row)

    missing: dict[str, object] = {}
    if not len(occupied):
        missing["occupied"] = "no_occupied_orbitals"
        missing["homo_ev"] = "no_occupied_orbitals"
    if not len(virtual):
        missing["virtual"] = "no_virtual_orbitals"
        missing["lumo_ev"] = "no_virtual_orbitals"
    if missing_atoms:
        missing["atoms"] = missing_atoms

    output = {
        "atoms": atoms,
        "homo_ev": homo_ev,
        "lumo_ev": lumo_ev,
        "window_ev": window,
        "homo_window_orb_idx": homo_win_idx,
        "lumo_window_orb_idx": lumo_win_idx,
        "scheme": scheme,
        "frontier_definition": "aufbau electron count",
        "orthonormality_error": float(pop["orthonormality_error"]),
        "sign_convention": (
            "occ_e_mean_ev is the occupation-weighted mean occupied MO energy in eV "
            "(typically negative for bound orbitals), equal to the negative of the "
            "condensed Average Local Ionization Energy (ALIE) proxy: ALIE_A = -occ_e_mean_ev."
        ),
        "analysis_settings": {
            "scheme": scheme,
            "frontier_definition": "aufbau electron count",
            "window_ev": window,
            "occ_e_mean_definition": "sum_{i in occ}(n_i * p_{A,i} * eps_i) / sum_{i in occ}(n_i * p_{A,i})",
            "vir_e_mean_definition": "sum_{i in vir}(p_{A,i} * eps_i) / sum_{i in vir}(p_{A,i})",
            "frontier_weight_definition": (
                "Gaussian decay exp(-((eps_frontier - eps_i) / window_ev)^2) anchored at HOMO/LUMO "
                "(exact degenerate frontier subspace when window_ev == 0)"
            ),
            "pdos_homo_win_definition": "sum_{i in occ, eps_HOMO - window_ev <= eps_i <= eps_HOMO}(n_i * p_{A,i})",
            "pdos_lumo_win_definition": "sum_{i in vir, eps_LUMO <= eps_i <= eps_LUMO + window_ev}(p_{A,i})",
        },
        "units": {
            "homo_ev": "eV",
            "lumo_ev": "eV",
            "window_ev": "eV",
            "occ_n_e": "e",
            "vir_n": "dimensionless",
            "occ_e_mean_ev": "eV",
            "vir_e_mean_ev": "eV",
            "occ_e_frontier_ev": "eV",
            "vir_e_frontier_ev": "eV",
            "pdos_homo_win_e": "e",
            "pdos_lumo_win": "dimensionless",
        },
        "invariance": (
            "invariant to AO permutation, MO phase, degenerate orbital subspace rotation, "
            "and rigid translation/rotation; preserves original calculation atom order"
        ),
        "missing": missing,
    }
    json.dumps(output, allow_nan=False)
    return output


def _overlap_basis(result, *, orthonormalityTolerance: float = 1e-5):
    if not _has_population_basis(result):
        raise DescriptorUnavailableError(
            "bond overlap descriptors require same-calculation overlap, coefficients, and AO-to-atom mapping; "
            "use XTBCalculator(...).Calculate(..., keepOverlap=True)"
        )
    energies, occupations, occupied, virtual = _orbitals(result)
    _validate_spin_summed_occupations(occupations)
    natoms = result.metadata.get("atom_count", len(result.atomicCharges))
    c, s, ao_atom = _orbital_basis_arrays(result, len(energies), natoms)
    eig, _ = np.linalg.eigh(s)
    if eig.min() <= max(1.0, eig.max()) * 1e-10:
        raise ValueError("overlap must be positive definite")
    metric = c.T @ s @ c
    err = float(np.max(np.abs(metric - np.eye(len(energies)))))
    if err > orthonormalityTolerance:
        raise ValueError(f"coefficients are not S-orthonormal: error={err:g}")
    return energies, occupations, occupied, virtual, c, s, ao_atom, err


def CalcBondOverlapDescriptors(
    mol,
    result,
):
    """Per-RDKit-bond COOP, eigenvalue-weighted COHP proxy, and charge polarity.

    For each RDKit bond ``(A, B)`` in original bond order:
    - ``coop_e``: total occupied Mulliken overlap population in e,
      ``2 * sum_i n_i * sum_{mu in A, nu in B} C_{mu,i} C_{nu,i} S_{mu,nu}``
      (positive for bonding, negative for antibonding).
    - ``ecohp_ev``: occupied eigenvalue-weighted Mulliken overlap population in eV,
      ``2 * sum_i n_i * eps_i * sum_{mu in A, nu in B} C_{mu,i} C_{nu,i} S_{mu,nu}``
      (negative for bonding stabilization, positive for antibonding). Note that
      ``Tr(W S) = sum_i n_i * eps_i`` using the energy-weighted
      density ``W = C diag(n * eps) C^T`` over all orbitals.
    - ``delta_q_e``: ``|q_A - q_B|`` in e.
    - ``signed_delta_q_e``: ``q_begin - q_end`` in e (RDKit bond orientation).
    - ``q_prod_e2``: ``q_A * q_B`` in e^2.
    - ``delta_cm5_q_e``, ``signed_delta_cm5_q_e``, ``cm5_q_prod_e2``: CM5
      polarity counterparts when available; ``None`` with explicit ``missing``
      reason (e.g. ``cm5_unavailable_for_gfn0``) otherwise.
    """
    result = _as_result(result)
    _validate_mapping(mol, result)
    energies, occupations, occupied, _, c, s, ao_atom, ortho_err = _overlap_basis(result)
    charges = _charges(result)
    if len(charges) != mol.GetNumAtoms():
        raise ValueError("bond overlap descriptors require atomic partial charges matching the molecule")
    cm5 = _cm5_charges(result)
    cm5_missing_reason = None
    if not len(cm5):
        cm5_missing_reason = (
            "cm5_unavailable_for_gfn0"
            if result.method == "GFN0-xTB"
            else "cm5_charges_unavailable"
        )

    occ_missing_reason = None if np.any(occupations > 0.0) else "no_occupied_orbitals"
    if occ_missing_reason is None:
        # All orbitals weighted by their (possibly fractional) occupations: the
        # exact density, independent of any occupied/virtual partition.
        p_mat = (c * occupations[None, :]) @ c.T
        w_mat = (c * (occupations * energies)[None, :]) @ c.T
        ps_elem = p_mat * s
        ws_elem = w_mat * s
    else:
        ps_elem = ws_elem = None

    atom_aos = [np.flatnonzero(ao_atom == a) for a in range(mol.GetNumAtoms())]
    rows = []
    for bond in mol.GetBonds():
        b_idx = int(bond.GetIdx())
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        row_missing: dict[str, str] = {}

        if occ_missing_reason is None:
            ia, ib = atom_aos[a], atom_aos[b]
            coop = float(2.0 * np.sum(ps_elem[np.ix_(ia, ib)]))
            ecohp = float(2.0 * np.sum(ws_elem[np.ix_(ia, ib)]))
        else:
            coop = ecohp = None
            row_missing["coop_e"] = occ_missing_reason
            row_missing["ecohp_ev"] = occ_missing_reason

        qa = float(charges[a])
        qb = float(charges[b])
        signed_dq = qa - qb
        dq = abs(signed_dq)
        q_prod = qa * qb

        if cm5_missing_reason is None:
            cm5_a = float(cm5[a])
            cm5_b = float(cm5[b])
            signed_dcm5 = cm5_a - cm5_b
            dcm5 = abs(signed_dcm5)
            cm5_prod = cm5_a * cm5_b
        else:
            signed_dcm5 = dcm5 = cm5_prod = None
            row_missing["delta_cm5_q_e"] = cm5_missing_reason
            row_missing["signed_delta_cm5_q_e"] = cm5_missing_reason
            row_missing["cm5_q_prod_e2"] = cm5_missing_reason

        rows.append({
            "bond_index": b_idx,
            "begin_atom_index": a,
            "end_atom_index": b,
            "coop_e": coop,
            "ecohp_ev": ecohp,
            "delta_q_e": float(dq),
            "signed_delta_q_e": float(signed_dq),
            "q_prod_e2": float(q_prod),
            "delta_cm5_q_e": float(dcm5) if dcm5 is not None else None,
            "signed_delta_cm5_q_e": float(signed_dcm5) if signed_dcm5 is not None else None,
            "cm5_q_prod_e2": float(cm5_prod) if cm5_prod is not None else None,
            "missing": row_missing,
        })

    missing: dict[str, str] = {}
    if occ_missing_reason is not None:
        missing["coop_e"] = occ_missing_reason
        missing["ecohp_ev"] = occ_missing_reason
    if cm5_missing_reason is not None:
        missing["delta_cm5_q_e"] = cm5_missing_reason
        missing["signed_delta_cm5_q_e"] = cm5_missing_reason
        missing["cm5_q_prod_e2"] = cm5_missing_reason

    output = {
        "bonds": rows,
        "charge_model": result.chargeModel,
        "frontier_definition": "aufbau electron count",
        "orthonormality_error": ortho_err,
        "sign_convention": {
            "coop_e": "positive for bonding overlap, negative for antibonding overlap",
            "ecohp_ev": "negative for bonding stabilization, positive for antibonding (eigenvalue-weighted COHP proxy)",
            "signed_delta_q_e": "q(begin_atom_index) - q(end_atom_index)",
        },
        "analysis_settings": {
            "frontier_definition": "aufbau electron count",
            "coop_definition": "2 * sum_i n_i * sum_{mu in A, nu in B} C_{mu,i} C_{nu,i} S_{mu,nu}",
            "ecohp_definition": "2 * sum_i n_i * eps_i * sum_{mu in A, nu in B} C_{mu,i} C_{nu,i} S_{mu,nu}",
            "charge_model": result.chargeModel,
        },
        "units": {
            "coop_e": "e",
            "ecohp_ev": "eV",
            "delta_q_e": "e",
            "signed_delta_q_e": "e",
            "q_prod_e2": "e^2",
            "delta_cm5_q_e": "e",
            "signed_delta_cm5_q_e": "e",
            "cm5_q_prod_e2": "e^2",
        },
        "invariance": (
            "coop_e, ecohp_ev, delta_q_e, q_prod_e2, delta_cm5_q_e, and cm5_q_prod_e2 "
            "are invariant to AO permutation, MO phase, degenerate occupied subspace rotation, "
            "bond direction reversal, and rigid translation/rotation; signed_delta_* antisymmetrizes "
            "under bond direction reversal; preserves original RDKit bond order"
        ),
        "missing": missing,
    }
    json.dumps(output, allow_nan=False)
    return output


_HARTREE_TO_EV = 27.211386245988


def CalcLinearResponseMatrix(
    result,
    *,
    scheme: str = "lowdin",
    gapToleranceEV: float = 1e-6,
    occupationTolerance: float = 1e-4,
):
    """Independent-particle atom-condensed linear response matrix chi_AB (e/eV, e/Ha).

    Evaluates the closed-shell Kohn-Sham / independent-particle linear response
    kernel condensed onto atoms (Geerlings et al., *Chem. Soc. Rev.* 2014, 43,
    4989):

        chi_AB^IP = -4 * sum_{i in occ} sum_{a in vir} (Q^A_{ia} * Q^B_{ia}) / (eps_a - eps_i)

    where ``Q^A_{ia}`` is the atom-resolved transition overlap matrix in either
    the symmetric orthogonalized Lowdin basis (default, ``Q^A_{ia} = sum_{mu in A}
    C_tilde_{mu,i} C_tilde_{mu,a}``) or symmetrized Mulliken partitioning
    (``Q^A_{ia} = 0.5 * sum_{mu in A} (C_{mu,i} (SC)_{mu,a} + (SC)_{mu,i} C_{mu,a})``).

    Properties:
    - Symmetric and negative semi-definite (``chi_AA <= 0``).
    - Exact electron-conservation sum rule: ``sum_B chi_AB = 0`` for every atom ``A``.
    - Explicitly withheld (``available=False``, ``None`` fields, explicit ``missing``
      reason) for open-shell, unknown-spin, fractionally occupied, or tiny/zero-gap
      calculations. Never invokes xTB.
    """
    result = _as_result(result)
    scheme = _validate_scheme(scheme)
    gap_tol = _validate_positive(gapToleranceEV, "gapToleranceEV")
    occ_tol = _validate_positive(occupationTolerance, "occupationTolerance")
    energies, occupations, occupied, virtual, c, s, ao_atom, ortho_err = _overlap_basis(result)
    natoms = int(result.metadata.get("atom_count", len(result.atomicCharges)))

    uhf = result.numUnpairedElectrons
    if uhf is None:
        uhf = result.metadata.get("uhf")

    withheld_reason: str | None = None
    if uhf is None:
        withheld_reason = "unknown_spin_state"
    elif int(uhf) != 0:
        withheld_reason = "open_shell_not_supported"
    elif np.any((occupations > occ_tol) & (np.abs(occupations - 2.0) > occ_tol)):
        withheld_reason = "fractional_occupation_not_supported"
    elif not len(occupied):
        withheld_reason = "no_occupied_orbitals"
    elif not len(virtual):
        withheld_reason = "no_virtual_orbitals"
    else:
        gap_ev = float(energies[virtual[0]] - energies[occupied[-1]])
        if gap_ev <= gap_tol:
            withheld_reason = "degenerate_homo_lumo_gap"

    units = {
        "chi_e_per_ev": "e/eV",
        "chi_e_per_ha": "e/Hartree",
        "sum_rule_max_error_e_per_ev": "e/eV",
        "homo_lumo_gap_ev": "eV",
        "chi_diag_e_per_ev": "e/eV",
        "chi_diag_e_per_ha": "e/Hartree",
        "neg_chi_diag_e_per_ev": "e/eV",
        "nonlocal_coupling_e_per_ev": "e/eV",
    }

    if withheld_reason is not None:
        missing = {
            "chi_e_per_ev": withheld_reason,
            "chi_e_per_ha": withheld_reason,
            "sum_rule_max_error_e_per_ev": withheld_reason,
            "atoms": withheld_reason,
        }
        output = {
            "available": False,
            "scheme": scheme,
            "chi_e_per_ev": None,
            "chi_e_per_ha": None,
            "sum_rule_max_error_e_per_ev": None,
            "homo_lumo_gap_ev": (
                float(energies[virtual[0]] - energies[occupied[-1]])
                if (len(occupied) and len(virtual))
                else None
            ),
            "n_occupied": int(len(occupied)),
            "n_virtual": int(len(virtual)),
            "orthonormality_error": ortho_err,
            "atoms": [],
            "analysis_settings": {
                "method": "independent_particle_condensed_linear_response",
                "reference": "Geerlings et al., Chem. Soc. Rev. 2014, 43, 4989",
                "scheme": scheme,
                "frontier_definition": "aufbau electron count",
                "gap_tolerance_ev": gap_tol,
                "occupation_tolerance": occ_tol,
                "definition": "chi_AB = -4 * sum_{i in occ, a in vir} Q^A_{ia} Q^B_{ia} / (eps_a - eps_i)",
                "coupled_response": False,
            },
            "units": units,
            "invariance": (
                "invariant to rigid translation/rotation, AO permutation, MO phase, "
                "and degenerate occupied/virtual subspace rotations; preserves original atom order"
            ),
            "missing": missing,
        }
        json.dumps(output, allow_nan=False)
        return output

    c_occ = c[:, occupied]
    c_vir = c[:, virtual]
    n_occ = len(occupied)
    n_vir = len(virtual)
    de_ev = energies[virtual][None, :] - energies[occupied][:, None]
    inv_de_flat = (1.0 / de_ev).ravel()

    q_flat = np.zeros((natoms, n_occ * n_vir), dtype=float)
    if scheme == "lowdin":
        eigvals, eigvecs = np.linalg.eigh(s)
        s_half = (eigvecs * np.sqrt(eigvals)) @ eigvecs.T
        ct_occ = s_half @ c_occ
        ct_vir = s_half @ c_vir
        for a in range(natoms):
            idx = np.flatnonzero(ao_atom == a)
            q_flat[a] = (ct_occ[idx].T @ ct_vir[idx]).ravel()
    else:
        sc_occ = s @ c_occ
        sc_vir = s @ c_vir
        for a in range(natoms):
            idx = np.flatnonzero(ao_atom == a)
            q_block = 0.5 * (c_occ[idx].T @ sc_vir[idx] + sc_occ[idx].T @ c_vir[idx])
            q_flat[a] = q_block.ravel()

    chi_ev = -4.0 * ((q_flat * inv_de_flat[None, :]) @ q_flat.T)
    chi_ev = 0.5 * (chi_ev + chi_ev.T)
    chi_ha = chi_ev * _HARTREE_TO_EV
    sum_rule_err = float(np.max(np.abs(np.sum(chi_ev, axis=1))))

    atoms = []
    for a in range(natoms):
        diag_ev = float(chi_ev[a, a])
        diag_ha = float(chi_ha[a, a])
        off_abs = float(np.sum(np.abs(chi_ev[a])) - abs(diag_ev))
        atoms.append({
            "atom_index": int(a),
            "chi_diag_e_per_ev": diag_ev,
            "chi_diag_e_per_ha": diag_ha,
            "neg_chi_diag_e_per_ev": float(-diag_ev),
            "nonlocal_coupling_e_per_ev": off_abs,
            "missing": {},
        })

    output = {
        "available": True,
        "scheme": scheme,
        "chi_e_per_ev": chi_ev.tolist(),
        "chi_e_per_ha": chi_ha.tolist(),
        "sum_rule_max_error_e_per_ev": sum_rule_err,
        "homo_lumo_gap_ev": float(energies[virtual[0]] - energies[occupied[-1]]),
        "n_occupied": int(n_occ),
        "n_virtual": int(n_vir),
        "orthonormality_error": ortho_err,
        "atoms": atoms,
        "analysis_settings": {
            "method": "independent_particle_condensed_linear_response",
            "reference": "Geerlings et al., Chem. Soc. Rev. 2014, 43, 4989",
            "scheme": scheme,
            "frontier_definition": "aufbau electron count",
            "gap_tolerance_ev": gap_tol,
            "occupation_tolerance": occ_tol,
            "definition": "chi_AB = -4 * sum_{i in occ, a in vir} Q^A_{ia} Q^B_{ia} / (eps_a - eps_i)",
            "coupled_response": False,
        },
        "units": units,
        "invariance": (
            "invariant to rigid translation/rotation, AO permutation, MO phase, "
            "and degenerate occupied/virtual subspace rotations; preserves original atom order"
        ),
        "missing": {},
    }
    json.dumps(output, allow_nan=False)
    return output
