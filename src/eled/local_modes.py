"""Cremer-Kraka local mode bond force constants and Grunenberg compliance constants.

Computes coordinate-invariant, normal-mode-decoupled local bond stretching force
constants ``k^a`` (in ``Hartree/Bohr^2`` and ``mDyn/Angstrom``) and relaxed
compliance constants ``C_ii = (B H^+ B^T)_ii`` (in ``Bohr^2/Hartree`` and
``Angstrom/mDyn``) following Konkoli & Cremer (1998), Zou et al. (*J. Chem.
Phys.* 2012, 137, 084114), and Grunenberg (*Chem. Soc. Rev.* 2017, 46, 5152).

- :func:`CalcLocalBondForceConstantsFromHessian`: pure analysis of an existing
  Cartesian Hessian (e.g. from ``CalcVibrationalThermochemistry(..., includeHessian=True)``)
  with zero extra xTB calls.
- :func:`CalcLocalBondForceConstants`: computes the analytical-gradient central
  finite-difference Cartesian Hessian via native ``libxtb`` and extracts local
  bond force constants for every RDKit bond plus any optional ``extraPairs``
  (e.g. hydrogen bonds, halogen bonds, or metal-ligand contacts).
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from typing import Any

import numpy as np

from .backends.native import ANGSTROM_TO_BOHR, _run_xtb_native
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .thermochemistry import (
    _AMU_TO_AU,
    _BOHR_TO_M,
    _HARTREE_TO_CM1,
    _HARTREE_TO_J,
    _atomic_masses_amu,
    _external_mode_basis,
    _positive_real,
)
from .xtb import _resolve_options

# 1 Hartree / Bohr^2 in mDyn / Angstrom (1 mDyn/A = 100 N/m).
HARTREE_BOHR2_TO_MDYN_A = (_HARTREE_TO_J / (_BOHR_TO_M ** 2)) / 100.0
LOCAL_MODES_PROTOCOL = 'eled-local-modes-v1'

__all__ = [
    'HARTREE_BOHR2_TO_MDYN_A',
    'CalcLocalBondForceConstantsFromHessian',
    'CalcLocalBondForceConstants',
]


def _unweighted_external_basis(coords_bohr: np.ndarray) -> np.ndarray:
    """Return orthonormal basis of rigid translations and rotations in unweighted Cartesian space."""
    natoms = len(coords_bohr)
    com = np.mean(coords_bohr, axis=0)
    rel = coords_bohr - com[None, :]
    cols = []
    for axis in range(3):
        v = np.zeros((natoms, 3), dtype=float)
        v[:, axis] = 1.0
        cols.append(v.reshape(-1))
    if natoms > 1:
        for axis in range(3):
            e = np.zeros(3, dtype=float)
            e[axis] = 1.0
            rot = np.cross(rel, e[None, :]).reshape(-1)
            norm = float(np.linalg.norm(rot))
            if norm > 1e-10:
                cols.append(rot / norm)
    ext = np.column_stack(cols)
    u, s, _ = np.linalg.svd(ext, full_matrices=False)
    rank = int(np.sum(s > 1e-10 * s[0]))
    return u[:, :rank]


def _validate_extra_pairs(extraPairs: Any, natoms: int) -> list[tuple[int, int]]:
    pairs: list[tuple[int, int]] = []
    for item in extraPairs:
        if len(item) != 2:
            raise ValueError('each entry in extraPairs must be a (atom_a, atom_b) pair')
        a, b = int(item[0]), int(item[1])
        if not (0 <= a < natoms and 0 <= b < natoms) or a == b:
            raise ValueError(f'extraPairs indices ({a}, {b}) out of range for {natoms} atoms')
        pairs.append((a, b))
    return pairs


def CalcLocalBondForceConstantsFromHessian(
    mol: Any,
    hessianHartreePerBohr2: Any,
    *,
    confId: int = -1,
    extraPairs: Any = (),
    eigenvalueTolerance: float = 1e-6,
    gradientHartreePerBohr: Any = None,
    stationaryGradientThresholdHaBohr: float = 1e-3,
    imaginaryThresholdCm1: float = 10.0,
) -> dict[str, Any]:
    """Compute Cremer-Kraka local force constants & compliance constants from a Cartesian Hessian.

    Never invokes xTB. ``hessianHartreePerBohr2`` must be a finite symmetric
    ``[3N, 3N]`` Cartesian Hessian in ``Hartree / Bohr^2`` matching ``mol`` atom
    ordering. Non-stationary geometries or imaginary modes are explicitly
    flagged via ``is_stationary_point`` and ``warnings``.
    """
    _xyz_from_mol(mol, confId)
    eig_tol = _positive_real(eigenvalueTolerance, 'eigenvalueTolerance')
    grad_thr = _positive_real(stationaryGradientThresholdHaBohr, 'stationaryGradientThresholdHaBohr')
    imag_thr = _positive_real(imaginaryThresholdCm1, 'imaginaryThresholdCm1', allow_zero=True)

    natoms = int(mol.GetNumAtoms())
    n3 = 3 * natoms
    hess = np.asarray(hessianHartreePerBohr2, dtype=float)
    if hess.shape != (n3, n3) or not np.isfinite(hess).all():
        raise ValueError(f'hessianHartreePerBohr2 must be a finite [{n3}, {n3}] array')
    if not np.allclose(hess, hess.T, atol=1e-6, rtol=1e-5):
        raise ValueError('hessianHartreePerBohr2 must be symmetric')
    hess = 0.5 * (hess + hess.T)

    conf = mol.GetConformer(confId)
    resolved_conf_id = int(conf.GetId())
    coords_ang = np.asarray(conf.GetPositions(), dtype=float)
    coords_bohr = coords_ang * ANGSTROM_TO_BOHR
    masses_amu = _atomic_masses_amu(mol)
    pairs = _validate_extra_pairs(extraPairs, natoms)

    # Check imaginary harmonic modes in mass-weighted coordinates (matching thermochemistry.py).
    q_vib, _, _, _ = _external_mode_basis(coords_bohr, masses_amu)
    n_imaginary = 0
    if q_vib.shape[1] > 0:
        inv_sqrt_m = np.repeat(1.0 / np.sqrt(masses_amu), 3)
        mw_hess = (hess * inv_sqrt_m[:, None]) * inv_sqrt_m[None, :]
        vib_eigs = np.linalg.eigvalsh(q_vib.T @ mw_hess @ q_vib)
        freqs_cm1 = [
            float(-math.sqrt(-ev / _AMU_TO_AU) * _HARTREE_TO_CM1) if ev < 0.0
            else float(math.sqrt(ev / _AMU_TO_AU) * _HARTREE_TO_CM1)
            for ev in vib_eigs
        ]
        n_imaginary = sum(1 for f in freqs_cm1 if f < -imag_thr)

    warnings: list[str] = []
    max_atom_grad: float | None = None
    if gradientHartreePerBohr is not None:
        g_arr = np.asarray(gradientHartreePerBohr, dtype=float)
        if g_arr.shape != (natoms, 3) or not np.isfinite(g_arr).all():
            raise ValueError(f'gradientHartreePerBohr must be a finite [{natoms}, 3] array')
        max_atom_grad = float(np.linalg.norm(g_arr, axis=1).max())
        if max_atom_grad > grad_thr:
            warnings.append(
                f'non_stationary_geometry: max atomic gradient {max_atom_grad:.4e} Ha/Bohr '
                f'exceeds threshold {grad_thr:.4e} Ha/Bohr'
            )
    if n_imaginary > 0:
        warnings.append(
            f'imaginary_modes_present: {n_imaginary} mode(s) below -{imag_thr:g} cm^-1; '
            'local force constants on non-minimum structures should be interpreted with caution'
        )
    is_stationary = (n_imaginary == 0) and (max_atom_grad is None or max_atom_grad <= grad_thr)

    # Project out unweighted rigid translations/rotations and compute exact internal pseudoinverse H_int^+.
    if natoms > 1:
        q_ext = _unweighted_external_basis(coords_bohr)
        proj = np.eye(n3) - q_ext @ q_ext.T
        h_int = proj @ hess @ proj
        h_int = 0.5 * (h_int + h_int.T)
        evals, evecs = np.linalg.eigh(h_int)
        inv_evals = np.where(np.abs(evals) > eig_tol, 1.0 / np.where(np.abs(evals) > eig_tol, evals, 1.0), 0.0)
        h_pinv = (evecs * inv_evals[None, :]) @ evecs.T
    else:
        h_pinv = np.zeros((n3, n3), dtype=float)

    def _eval_pair(a: int, b: int) -> dict[str, Any]:
        diff = coords_ang[b] - coords_ang[a]
        length = float(np.linalg.norm(diff))
        if not math.isfinite(length) or length <= 1e-6:
            raise ValueError(f'degenerate atom pair ({a}, {b})')
        unit = diff / length
        bvec = np.zeros(n3, dtype=float)
        bvec[3 * a : 3 * a + 3] = -unit
        bvec[3 * b : 3 * b + 3] = unit
        comp_au = float(bvec @ h_pinv @ bvec)
        mu_amu = float((masses_amu[a] * masses_amu[b]) / (masses_amu[a] + masses_amu[b]))
        row_missing: dict[str, str] = {}
        if abs(comp_au) <= 1e-12:
            k_au = k_mdyn = comp_mdyn = omega_cm1 = None
            for k in ('k_local_ha_bohr2', 'k_local_mdyn_a', 'compliance_a_mdyn', 'omega_local_cm1'):
                row_missing[k] = 'singular_zero_compliance'
        else:
            k_au = float(1.0 / comp_au)
            k_mdyn = float(k_au * HARTREE_BOHR2_TO_MDYN_A)
            comp_mdyn = float(comp_au / HARTREE_BOHR2_TO_MDYN_A)
            if k_au > 0.0:
                omega_cm1: float | None = float(math.sqrt(k_au / (mu_amu * _AMU_TO_AU)) * _HARTREE_TO_CM1)
            else:
                omega_cm1 = None
                row_missing['omega_local_cm1'] = 'nonpositive_local_force_constant'
        return {
            'bond_length_a': length,
            'reduced_mass_amu': mu_amu,
            'compliance_bohr2_ha': comp_au,
            'compliance_a_mdyn': comp_mdyn,
            'k_local_ha_bohr2': k_au,
            'k_local_mdyn_a': k_mdyn,
            'omega_local_cm1': omega_cm1,
            'missing': row_missing,
        }

    bond_rows = []
    for bond in mol.GetBonds():
        b_idx = int(bond.GetIdx())
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        bond_rows.append({
            'bond_index': b_idx,
            'begin_atom_index': a,
            'end_atom_index': b,
            **_eval_pair(a, b),
        })

    pair_rows = []
    for p_idx, (a, b) in enumerate(pairs):
        pair_rows.append({
            'pair_index': p_idx,
            'atom_a': a,
            'atom_b': b,
            **_eval_pair(a, b),
        })

    missing: dict[str, str] = {}
    if not bond_rows:
        missing['bonds'] = 'molecule_has_no_bonds'

    payload = {
        'n_bonds': len(bond_rows),
        'n_extra_pairs': len(pair_rows),
        'is_stationary_point': is_stationary,
        'n_imaginary_modes': int(n_imaginary),
        'max_atom_gradient_ha_bohr': max_atom_grad,
        'warnings': warnings,
        'bonds': bond_rows,
        'extra_pairs': pair_rows,
        'units': {
            'bond_length_a': 'angstrom',
            'reduced_mass_amu': 'amu',
            'compliance_bohr2_ha': 'Bohr^2/Hartree',
            'compliance_a_mdyn': 'angstrom/mDyn',
            'k_local_ha_bohr2': 'Hartree/Bohr^2',
            'k_local_mdyn_a': 'mDyn/angstrom',
            'omega_local_cm1': 'cm^-1',
            'max_atom_gradient_ha_bohr': 'Hartree/Bohr',
        },
        'analysis_settings': {
            'protocol': LOCAL_MODES_PROTOCOL,
            'reference': 'Konkoli & Cremer (1998); Zou et al. JCP 2012; Grunenberg Chem. Soc. Rev. 2017',
            'definition': 'C_ii = b_i^T H_int^+ b_i; k^a = 1 / C_ii; omega^a = (1/2pi*c)*sqrt(k^a / mu_AB)',
            'hartree_bohr2_to_mdyn_a': HARTREE_BOHR2_TO_MDYN_A,
            'eigenvalue_tolerance': eig_tol,
            'stationary_gradient_threshold_ha_bohr': grad_thr,
            'imaginary_threshold_cm1': imag_thr,
            'conformer_id': resolved_conf_id,
        },
        'invariance': (
            'compliance constants, local force constants and local frequencies are strictly '
            'coordinate-system, rigid-motion, normal-mode-coupling, and bond-reversal invariant; '
            'preserves RDKit bond order'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))


def CalcLocalBondForceConstants(
    mol: Any,
    *,
    confId: int = -1,
    extraPairs: Any = (),
    stepBohr: float = 0.005,
    eigenvalueTolerance: float = 1e-6,
    stationaryGradientThresholdHaBohr: float = 1e-3,
    imaginaryThresholdCm1: float = 10.0,
    options: XTBOptions | None = None,
    method: str = 'GFN2-xTB',
    charge: int | None = None,
    numUnpairedElectrons: int | None = None,
    xtbLibrary: Any = None,
    solvent: str | None = None,
    solvationModel: str = 'alpb',
    solvationGrid: int | None = None,
    maxIterations: int = 250,
    accuracy: float = 1.0,
    electronicTemperature: float = 300.0,
) -> dict[str, Any]:
    """Compute Cremer-Kraka local bond force constants & compliance constants via native xTB Hessian."""
    from rdkit import Chem

    step = _positive_real(stepBohr, 'stepBohr')
    if options is not None:
        if not isinstance(options, XTBOptions):
            raise TypeError('options must be an XTBOptions instance or None')
        opts = options
    else:
        opts = XTBOptions(
            confId=confId,
            method=method,
            charge=charge,
            numUnpairedElectrons=numUnpairedElectrons,
            xtbLibrary=xtbLibrary,
            maxIterations=maxIterations,
            accuracy=accuracy,
            electronicTemperature=electronicTemperature,
            solvent=solvent,
            solvationModel=solvationModel,
            solvationGrid=solvationGrid,
        )

    working = Chem.Mol(mol)
    _xyz_from_mol(working, opts.confId)
    eff_charge, eff_unpaired = _effective_charge_uhf(working, opts.charge, opts.numUnpairedElectrons)
    resolved = _resolve_options(replace(opts, charge=eff_charge, numUnpairedElectrons=eff_unpaired))
    natoms = working.GetNumAtoms()
    _validate_extra_pairs(extraPairs, natoms)

    run_kwargs = dict(
        confId=resolved.confId,
        method=resolved.method,
        charge=resolved.charge,
        numUnpairedElectrons=resolved.numUnpairedElectrons,
        xtbLibrary=resolved.xtbLibrary,
        maxIterations=resolved.maxIterations,
        accuracy=resolved.accuracy,
        electronicTemperature=resolved.electronicTemperature,
        solvent=resolved.solvent,
        solvationModel=resolved.solvationModel,
        solvationGrid=resolved.solvationGrid,
        _return_gradient=True,
    )

    eq_res = _run_xtb_native(working, **run_kwargs)
    grad0 = eq_res.GetGradient()
    coords_ang = np.asarray(working.GetConformer(resolved.confId).GetPositions(), dtype=float)
    n3 = 3 * natoms
    hess = np.zeros((n3, n3), dtype=float)
    step_ang = step / ANGSTROM_TO_BOHR
    if natoms > 1:
        for atom_idx in range(natoms):
            for axis in range(3):
                col = 3 * atom_idx + axis
                grads = []
                for sign in (+1.0, -1.0):
                    displaced = Chem.Mol(working)
                    conf = displaced.GetConformer(resolved.confId)
                    pos = coords_ang[atom_idx].copy()
                    pos[axis] += sign * step_ang
                    conf.SetAtomPosition(atom_idx, pos.tolist())
                    res_d = _run_xtb_native(displaced, **run_kwargs)
                    grads.append(res_d.GetGradient().reshape(-1))
                hess[:, col] = (grads[0] - grads[1]) / (2.0 * step)
        hess = 0.5 * (hess + hess.T)

    report = CalcLocalBondForceConstantsFromHessian(
        working,
        hess,
        confId=resolved.confId,
        extraPairs=extraPairs,
        eigenvalueTolerance=eigenvalueTolerance,
        gradientHartreePerBohr=grad0,
        stationaryGradientThresholdHaBohr=stationaryGradientThresholdHaBohr,
        imaginaryThresholdCm1=imaginaryThresholdCm1,
    )
    report['method'] = resolved.method
    report['backend'] = 'native'
    report['analysis_settings']['step_bohr'] = step
    report['analysis_settings']['solvent'] = resolved.solvent
    report['analysis_settings']['solvation_model'] = resolved.solvationModel if resolved.solvent else None
    return report
