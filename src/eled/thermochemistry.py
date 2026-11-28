"""Harmonic vibrational frequency analysis and RRHO/qRRHO thermochemistry.

Computes the molecular Cartesian Hessian by central finite differences of
analytical xTB nuclear gradients (``GetGradient()``), projects out external
translations and rotations in mass-weighted coordinates, and evaluates
harmonic frequencies, IR intensities, ZPE, thermal enthalpy, entropy, and
Gibbs free energy using ideal-gas RRHO (with optional Grimme quasi-RRHO
low-frequency entropy interpolation). Non-stationary geometries (large nuclear
gradients or imaginary vibrational modes) are explicitly flagged.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from numbers import Integral, Real
from typing import Any

import numpy as np

from .backends.native import ANGSTROM_TO_BOHR, _run_xtb_native
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .xtb import _resolve_options

# CODATA 2018 exact/standard physical constants.
_KB_SI = 1.380649e-23            # J / K
_H_SI = 6.62607015e-34           # J * s
_C_SI = 2.99792458e8             # m / s
_NA = 6.02214076e23              # 1 / mol
_R_SI = _NA * _KB_SI             # J / (mol * K)
_HARTREE_TO_J = 4.3597447222071e-18
_HARTREE_TO_J_MOL = _HARTREE_TO_J * _NA
_CAL_TO_J = 4.184
_HARTREE_TO_KCAL_MOL = _HARTREE_TO_J_MOL / (1000.0 * _CAL_TO_J)
_AMU_TO_KG = 1.66053906660e-27
_ME_KG = 9.1093837015e-31
_AMU_TO_AU = _AMU_TO_KG / _ME_KG
_BOHR_TO_M = 5.29177210903e-11
_ATM_TO_PA = 101325.0
_HARTREE_TO_CM1 = 219474.6313632
_E_CHARGE_SI = 1.602176634e-19
_EPS0_SI = 8.8541878128e-12
# Exact IUPAC Naperian IR intensity prefactor from (dmu/dQ) in e / sqrt(amu) to km/mol:
# I = (N_A * pi) / (3 * 4*pi*eps0 * c^2) * (e^2 / amu_kg) / 1000
_IR_INTENSITY_FACTOR_KM_MOL = (
    (_NA * math.pi)
    / (3.0 * 4.0 * math.pi * _EPS0_SI * (_C_SI ** 2))
    * ((_E_CHARGE_SI ** 2) / _AMU_TO_KG)
    / 1000.0
)


def _positive_real(value: Any, name: str, *, allow_zero: bool = False) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(float(value)):
        raise ValueError(f'{name} must be a finite real number')
    out = float(value)
    if (allow_zero and out < 0.0) or (not allow_zero and out <= 0.0):
        qual = 'nonnegative' if allow_zero else 'positive'
        raise ValueError(f'{name} must be {qual}')
    return out


def _positive_int(value: Any, name: str) -> int:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or int(value) < 1:
        raise ValueError(f'{name} must be a positive integer')
    return int(value)


def _atomic_masses_amu(mol: Any) -> np.ndarray:
    from rdkit import Chem
    pt = Chem.GetPeriodicTable()
    masses = []
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        mass = pt.GetMassForIsotope(atom.GetAtomicNum(), isotope) if isotope else pt.GetAtomicWeight(atom.GetAtomicNum())
        if not math.isfinite(mass) or mass <= 0.0:
            raise ValueError(f'invalid atomic mass for atom index {atom.GetIdx()}')
        masses.append(float(mass))
    return np.asarray(masses, dtype=float)


def _external_mode_basis(coords_bohr: np.ndarray, masses_amu: np.ndarray) -> tuple[np.ndarray, np.ndarray, bool, bool]:
    """Return orthonormal mass-weighted external basis, principal inertias (amu*Bohr^2), is_linear, is_atom."""
    natoms = len(masses_amu)
    sqrt_m = np.sqrt(masses_amu)
    total_mass = float(np.sum(masses_amu))
    com = np.sum(masses_amu[:, None] * coords_bohr, axis=0) / total_mass
    rel = coords_bohr - com[None, :]

    # Moment of inertia tensor in amu * Bohr^2.
    inertia_tensor = np.zeros((3, 3), dtype=float)
    for a in range(natoms):
        r = rel[a]
        r2 = float(np.dot(r, r))
        inertia_tensor += masses_amu[a] * (r2 * np.eye(3) - np.outer(r, r))
    principal_moments, principal_axes = np.linalg.eigh(inertia_tensor)
    principal_moments = np.maximum(principal_moments, 0.0)

    is_atom = bool(natoms == 1)
    max_i = float(principal_moments[-1]) if not is_atom else 0.0
    is_linear = bool((not is_atom) and (principal_moments[0] < 1.0e-4 * max(max_i, 1.0)))

    cols = []
    # 3 translational basis vectors.
    for axis in range(3):
        vec = np.zeros((natoms, 3), dtype=float)
        vec[:, axis] = sqrt_m
        cols.append(vec.reshape(-1))
    # Rotational basis vectors around principal axes with non-negligible inertia.
    if not is_atom:
        for j in range(3):
            if is_linear and j == 0:
                continue
            if principal_moments[j] <= 1.0e-10:
                continue
            axis_vec = principal_axes[:, j]
            rot = np.cross(rel, axis_vec[None, :]) * sqrt_m[:, None]
            norm = float(np.linalg.norm(rot))
            if norm > 1.0e-10:
                cols.append((rot / norm).reshape(-1))

    ext = np.column_stack(cols)
    # Orthonormalize via SVD so external basis is exact to machine precision.
    u, s, _ = np.linalg.svd(ext, full_matrices=True)
    rank = int(np.sum(s > 1.0e-10 * s[0]))
    q_ext = u[:, :rank]
    q_vib = u[:, rank:]
    return q_vib, principal_moments, is_linear, is_atom


def CalcVibrationalThermochemistry(
    mol: Any,
    *,
    confId: int = -1,
    options: XTBOptions | None = None,
    temperatureK: float = 298.15,
    pressureAtm: float = 1.0,
    method: str = 'GFN2-xTB',
    stepBohr: float = 0.005,
    imaginaryThresholdCm1: float = 10.0,
    stationaryGradientThresholdHaBohr: float = 1.0e-3,
    rotorCutoffCm1: float = 100.0,
    symmetryNumber: int = 1,
    charge: int | None = None,
    numUnpairedElectrons: int | None = None,
    xtbLibrary: Any = None,
    solvent: str | None = None,
    solvationModel: str = 'alpb',
    solvationGrid: int | None = None,
    maxIterations: int = 250,
    accuracy: float = 1.0,
    electronicTemperature: float = 300.0,
    includeHessian: bool = False,
) -> dict[str, Any]:
    """Compute harmonic vibrational frequencies and RRHO/qRRHO thermochemistry.

    Uses central finite differences of analytical xTB nuclear gradients
    (``GetGradient()``) at step ``stepBohr``. Non-stationary geometries (atom
    gradient norm exceeding ``stationaryGradientThresholdHaBohr`` or imaginary
    modes below ``-imaginaryThresholdCm1``) are explicitly flagged via
    ``is_stationary_point`` and ``warnings``.

    When ``options`` (an ``XTBOptions`` instance) is provided, it takes
    precedence over individual single-point calculation keyword arguments
    (``confId``, ``method``, ``charge``, ``numUnpairedElectrons``,
    ``xtbLibrary``, ``solvent``, ``solvationModel``, ``solvationGrid``,
    ``maxIterations``, ``accuracy``, ``electronicTemperature``).
    """
    from rdkit import Chem

    temp = _positive_real(temperatureK, 'temperatureK')
    press = _positive_real(pressureAtm, 'pressureAtm')
    step = _positive_real(stepBohr, 'stepBohr')
    if not isinstance(includeHessian, bool):
        raise ValueError('includeHessian must be bool')
    imag_thr = _positive_real(imaginaryThresholdCm1, 'imaginaryThresholdCm1', allow_zero=True)
    grad_thr = _positive_real(stationaryGradientThresholdHaBohr, 'stationaryGradientThresholdHaBohr')
    sthr = _positive_real(rotorCutoffCm1, 'rotorCutoffCm1', allow_zero=True)
    sigma = _positive_int(symmetryNumber, 'symmetryNumber')

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
    e_elec_ha = float(eq_res.totalEnergyHartree)
    grad0 = eq_res.GetGradient()
    natoms = working.GetNumAtoms()
    grad_norm = float(np.linalg.norm(grad0))
    max_atom_grad = float(np.linalg.norm(grad0, axis=1).max())

    coords_ang = np.asarray(working.GetConformer(resolved.confId).GetPositions(), dtype=float)
    coords_bohr = coords_ang * ANGSTROM_TO_BOHR
    masses_amu = _atomic_masses_amu(working)
    total_mass_amu = float(np.sum(masses_amu))

    # Central finite-difference Cartesian Hessian (Hartree / Bohr^2) and dipole derivatives (e).
    n3 = 3 * natoms
    hess = np.zeros((n3, n3), dtype=float)
    dip_grad = np.zeros((3, n3), dtype=float)
    step_ang = step / ANGSTROM_TO_BOHR

    if natoms > 1:
        for atom_idx in range(natoms):
            for axis in range(3):
                col = 3 * atom_idx + axis
                grads = []
                dips = []
                for sign in (+1.0, -1.0):
                    displaced = Chem.Mol(working)
                    conf = displaced.GetConformer(resolved.confId)
                    pos = coords_ang[atom_idx].copy()
                    pos[axis] += sign * step_ang
                    conf.SetAtomPosition(atom_idx, pos.tolist())
                    res_d = _run_xtb_native(displaced, **run_kwargs)
                    grads.append(res_d.GetGradient().reshape(-1))
                    dips.append(res_d.GetDipole())
                hess[:, col] = (grads[0] - grads[1]) / (2.0 * step)
                dip_grad[:, col] = (dips[0] - dips[1]) / (2.0 * step)
        hess = 0.5 * (hess + hess.T)

    q_vib, principal_moments_amu_bohr2, is_linear, is_atom = _external_mode_basis(coords_bohr, masses_amu)
    n_vib_modes = q_vib.shape[1]

    if n_vib_modes > 0:
        inv_sqrt_m = np.repeat(1.0 / np.sqrt(masses_amu), 3)
        mw_hess = (hess * inv_sqrt_m[:, None]) * inv_sqrt_m[None, :]
        mw_hess = 0.5 * (mw_hess + mw_hess.T)
        h_sub = q_vib.T @ mw_hess @ q_vib
        h_sub = 0.5 * (h_sub + h_sub.T)
        eigvals, eigvecs_sub = np.linalg.eigh(h_sub)
        freqs_cm1 = [
            float(math.copysign(math.sqrt(abs(ev) / _AMU_TO_AU) * _HARTREE_TO_CM1, ev))
            for ev in eigvals
        ]
        # Normal-mode Cartesian vectors in amu^(-1/2): L = M^(-1/2) Q_vib V.
        l_cart = inv_sqrt_m[:, None] * (q_vib @ eigvecs_sub)
        dmu_dq = dip_grad @ l_cart  # shape (3, n_vib_modes), units e / sqrt(amu)
        ir_intensities = [
            float(_IR_INTENSITY_FACTOR_KM_MOL * float(np.dot(dmu_dq[:, k], dmu_dq[:, k])))
            for k in range(n_vib_modes)
        ]
    else:
        freqs_cm1 = []
        ir_intensities = []

    imag_freqs = [nu for nu in freqs_cm1 if nu < -imag_thr]
    n_imaginary = len(imag_freqs)
    real_freqs = [nu for nu in freqs_cm1 if nu > 0.0]

    warnings: list[str] = []
    if max_atom_grad > grad_thr:
        warnings.append('non_stationary_geometry_large_gradient')
    if n_imaginary > 0:
        warnings.append('imaginary_vibrational_modes_present')
    if any(-imag_thr <= nu <= 0.0 for nu in freqs_cm1):
        warnings.append('near_zero_or_small_negative_modes_excluded_from_vibrational_thermo')
    is_stationary = (max_atom_grad <= grad_thr) and (n_imaginary == 0)

    # Translational thermochemistry (J/mol and J/(mol*K)).
    mass_kg = total_mass_amu * _AMU_TO_KG
    pressure_pa = press * _ATM_TO_PA
    q_trans = ((2.0 * math.pi * mass_kg * _KB_SI * temp) / (_H_SI ** 2)) ** 1.5 * (_KB_SI * temp / pressure_pa)
    u_trans_j = 1.5 * _R_SI * temp
    s_trans_j = _R_SI * (math.log(q_trans) + 2.5)

    # Rotational thermochemistry.
    inertias_si = principal_moments_amu_bohr2 * _AMU_TO_KG * (_BOHR_TO_M ** 2)
    rot_constants_cm1: list[float] = []
    if is_atom:
        u_rot_j = 0.0
        s_rot_j = 0.0
    elif is_linear:
        i_lin = float(inertias_si[-1])
        b_cm1 = _H_SI / (8.0 * (math.pi ** 2) * _C_SI * i_lin * 100.0)
        rot_constants_cm1 = [b_cm1]
        theta_rot = (_H_SI ** 2) / (8.0 * (math.pi ** 2) * i_lin * _KB_SI)
        q_rot = temp / (sigma * theta_rot)
        u_rot_j = _R_SI * temp
        s_rot_j = _R_SI * (math.log(q_rot) + 1.0)
    else:
        thetas = []
        for i_val in inertias_si:
            b_cm1 = _H_SI / (8.0 * (math.pi ** 2) * _C_SI * float(i_val) * 100.0)
            rot_constants_cm1.append(b_cm1)
            thetas.append((_H_SI ** 2) / (8.0 * (math.pi ** 2) * float(i_val) * _KB_SI))
        q_rot = (math.sqrt(math.pi) / sigma) * math.sqrt((temp ** 3) / (thetas[0] * thetas[1] * thetas[2]))
        u_rot_j = 1.5 * _R_SI * temp
        s_rot_j = _R_SI * (math.log(q_rot) + 1.5)

    # Electronic spin entropy.
    multiplicity = int(eq_res.numUnpairedElectrons or 0) + 1
    s_elec_j = _R_SI * math.log(multiplicity)

    # Vibrational thermochemistry (over positive real modes).
    b_av_si = float(np.mean(inertias_si)) if (not is_atom and float(np.mean(inertias_si)) > 0.0) else 1.0e-44
    zpe_j = 0.0
    u_vib_th_j = 0.0
    s_vib_j = 0.0
    for nu in real_freqs:
        theta_v = (_H_SI * _C_SI * 100.0 * nu) / _KB_SI
        x = theta_v / temp
        zpe_j += 0.5 * _R_SI * theta_v
        if x < 500.0:
            ex = math.exp(x)
            emx = math.exp(-x)
            u_vib_th_j += _R_SI * theta_v / (ex - 1.0)
            s_ho = _R_SI * (x / (ex - 1.0) - math.log(1.0 - emx))
        else:
            s_ho = 0.0
        if sthr > 0.0:
            w = 1.0 / (1.0 + (sthr / nu) ** 4)
            mu = _H_SI / (8.0 * (math.pi ** 2) * _C_SI * 100.0 * nu)
            mu_eff = (mu * b_av_si) / (mu + b_av_si)
            s_free = _R_SI * (0.5 + math.log(math.sqrt(8.0 * (math.pi ** 3) * mu_eff * _KB_SI * temp / (_H_SI ** 2))))
            s_vib_j += w * s_ho + (1.0 - w) * s_free
        else:
            s_vib_j += s_ho

    # Assemble totals in J/mol, Hartree, and kcal/mol.
    u_corr_j = zpe_j + u_vib_th_j + u_rot_j + u_trans_j
    h_corr_j = u_corr_j + _R_SI * temp
    s_tot_j = s_trans_j + s_rot_j + s_vib_j + s_elec_j
    ts_j = temp * s_tot_j
    g_corr_j = h_corr_j - ts_j

    def _j_to_ha(val: float) -> float:
        return val / _HARTREE_TO_J_MOL

    def _j_to_kcal(val: float) -> float:
        return val / (1000.0 * _CAL_TO_J)

    zpe_ha = _j_to_ha(zpe_j)
    u_corr_ha = _j_to_ha(u_corr_j)
    h_corr_ha = _j_to_ha(h_corr_j)
    ts_ha = _j_to_ha(ts_j)
    g_corr_ha = _j_to_ha(g_corr_j)

    total_u_ha = e_elec_ha + u_corr_ha
    total_h_ha = e_elec_ha + h_corr_ha
    total_g_ha = e_elec_ha + g_corr_ha

    payload = {
        'method': eq_res.method,
        'solvent': eq_res.solvent,
        'solvation_model': eq_res.solvationModel,
        'thermo_model': 'qrrho' if sthr > 0.0 else 'rrho',
        'temperature_k': temp,
        'pressure_atm': press,
        'symmetry_number': sigma,
        'is_atom': is_atom,
        'is_linear': is_linear,
        'is_stationary_point': is_stationary,
        'warnings': warnings,
        'gradient_norm_ha_bohr': grad_norm,
        'max_atom_gradient_ha_bohr': max_atom_grad,
        'n_atoms': natoms,
        'molecular_mass_amu': total_mass_amu,
        'n_vibrational_modes': n_vib_modes,
        'n_imaginary': n_imaginary,
        'frequencies_cm1': freqs_cm1,
        'imaginary_frequencies_cm1': imag_freqs,
        'ir_intensities_km_mol': ir_intensities,
        'rotational_constants_cm1': rot_constants_cm1,
        'electronic_energy_ha': e_elec_ha,
        'zpe_ha': zpe_ha,
        'zpe_kcalmol': _j_to_kcal(zpe_j),
        'thermal_energy_correction_ha': u_corr_ha,
        'thermal_energy_correction_kcalmol': _j_to_kcal(u_corr_j),
        'thermal_enthalpy_correction_ha': h_corr_ha,
        'thermal_enthalpy_correction_kcalmol': _j_to_kcal(h_corr_j),
        'entropy_j_mol_k': s_tot_j,
        'entropy_cal_mol_k': s_tot_j / _CAL_TO_J,
        'entropy_components_j_mol_k': {
            'translational': s_trans_j,
            'rotational': s_rot_j,
            'vibrational': s_vib_j,
            'electronic': s_elec_j,
        },
        'ts_ha': ts_ha,
        'ts_kcalmol': _j_to_kcal(ts_j),
        'gibbs_correction_ha': g_corr_ha,
        'gibbs_correction_kcalmol': _j_to_kcal(g_corr_j),
        'total_energy_u_ha': total_u_ha,
        'total_energy_u_kcalmol': total_u_ha * _HARTREE_TO_KCAL_MOL,
        'total_enthalpy_ha': total_h_ha,
        'total_enthalpy_kcalmol': total_h_ha * _HARTREE_TO_KCAL_MOL,
        'total_gibbs_free_energy_ha': total_g_ha,
        'total_gibbs_free_energy_kcalmol': total_g_ha * _HARTREE_TO_KCAL_MOL,
        'settings': {
            'step_bohr': step,
            'imaginary_threshold_cm1': imag_thr,
            'stationary_gradient_threshold_ha_bohr': grad_thr,
            'rotor_cutoff_cm1': sthr,
            'symmetry_number': sigma,
            'temperature_k': temp,
            'pressure_atm': press,
            'hessian_method': 'central_finite_difference_analytical_gradient',
        },
        'provenance': {
            'conformer_id': eq_res.conformerId,
            'charge': eq_res.charge,
            'num_unpaired_electrons': eq_res.numUnpairedElectrons,
            'geometry_sha256': eq_res.metadata.get('geometry_sha256'),
            'ordered_structure_sha256': eq_res.metadata.get('ordered_structure_sha256'),
            'calculation_settings': eq_res.metadata.get('calculation_settings'),
        },
        'units': {
            'temperature_k': 'K',
            'pressure_atm': 'atm',
            'gradient_norm_ha_bohr': 'Hartree/Bohr',
            'max_atom_gradient_ha_bohr': 'Hartree/Bohr',
            'molecular_mass_amu': 'amu',
            'frequencies_cm1': 'cm^-1',
            'imaginary_frequencies_cm1': 'cm^-1',
            'ir_intensities_km_mol': 'km/mol',
            'rotational_constants_cm1': 'cm^-1',
            'electronic_energy_ha': 'Hartree',
            'zpe_ha': 'Hartree',
            'zpe_kcalmol': 'kcal/mol',
            'thermal_energy_correction_ha': 'Hartree',
            'thermal_energy_correction_kcalmol': 'kcal/mol',
            'thermal_enthalpy_correction_ha': 'Hartree',
            'thermal_enthalpy_correction_kcalmol': 'kcal/mol',
            'entropy_j_mol_k': 'J/(mol*K)',
            'entropy_cal_mol_k': 'cal/(mol*K)',
            'ts_ha': 'Hartree',
            'ts_kcalmol': 'kcal/mol',
            'gibbs_correction_ha': 'Hartree',
            'gibbs_correction_kcalmol': 'kcal/mol',
            'total_energy_u_ha': 'Hartree',
            'total_energy_u_kcalmol': 'kcal/mol',
            'total_enthalpy_ha': 'Hartree',
            'total_enthalpy_kcalmol': 'kcal/mol',
            'total_gibbs_free_energy_ha': 'Hartree',
            'total_gibbs_free_energy_kcalmol': 'kcal/mol',
        },
        'invariance': (
            'harmonic frequencies, IR intensities, rotational constants and thermochemical '
            'quantities are rigid-translation, rotation and atom-permutation invariant'
        ),
        'missing': {},
    }
    if includeHessian:
        payload['hessian_ha_bohr2'] = hess.tolist()
        payload['units']['hessian_ha_bohr2'] = 'Hartree/Bohr^2'
    return json.loads(json.dumps(payload, allow_nan=False))
