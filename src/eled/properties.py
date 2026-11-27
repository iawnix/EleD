"""Extended molecular and atomic electronic properties from the libxtb properties ABI.

Every descriptor function in this module is a pure analysis of an existing
single-point result (either native ``XTBResults`` or compiled
``rdXTBTools.XTBResults``). Unavailable or method-inapplicable quantities are
returned as ``None`` with explicit ``missing`` reasons and are never encoded as
zero.
"""
from __future__ import annotations

import json
import math
from typing import Any

import numpy as np

from .electronic_descriptors import AU_DIPOLE_TO_DEBYE, _as_result
from .errors import DescriptorUnavailableError
from .mapping import _validate_mapping
from .result import XTBResults

BOHR_TO_ANGSTROM = 0.529177210903
BOHR3_TO_ANGSTROM3 = BOHR_TO_ANGSTROM ** 3
HARTREE_TO_KCAL_MOL = 627.5094740631
HARTREE_TO_EV = 27.211386245988

_QUADRUPOLE_COMPONENTS = ('xx', 'xy', 'yy', 'xz', 'yz', 'zz')


def _is_compiled_result(result: Any) -> bool:
    return (
        type(result).__module__ == 'eled.rdXTBTools'
        and type(result).__name__ == 'XTBResults'
    )


def _as_property_result(result: Any) -> XTBResults:
    """Normalize native or compiled XTBResults while preserving extended properties."""
    if isinstance(result, XTBResults):
        if not result.success:
            detail = result.metadata.get('error_message', 'calculation failed')
            raise ValueError(f'calculation failed: {detail}')
        return result
    base = _as_result(result)
    if not _is_compiled_result(result):
        return base

    # Copy extended property fields directly from compiled rdXTBTools.XTBResults.
    extra: dict[str, Any] = {}
    for attr in (
        'solvent',
        'solvationModel',
        'solvationGrid',
        'electronicEnergyHartree',
        'repulsionEnergyHartree',
        'electrostaticEnergyHartree',
        'anisotropicElectrostaticEnergyHartree',
        'anisotropicXCEnergyHartree',
        'dispersionEnergyHartree',
        'halogenBondEnergyHartree',
        'solvationEnergyHartree',
        'solvationBornEnergyHartree',
        'solvationSASAEnergyHartree',
        'solvationHBEnergyHartree',
        'solvationShiftEnergyHartree',
        'fermiLevelAlphaHartree',
        'fermiLevelBetaHartree',
        'fermiLevelAlphaEV',
        'fermiLevelBetaEV',
        'molecularPolarizabilityAU',
        'molecularC6AU',
        'molecularC8AU',
    ):
        if hasattr(result, attr):
            extra[attr] = getattr(result, attr)

    if getattr(result, 'energyDecompositionAvailable', False):
        flags = 1
        if extra.get('electronicEnergyHartree') is not None:
            flags |= 2
        if extra.get('repulsionEnergyHartree') is not None:
            flags |= 4
        if extra.get('electrostaticEnergyHartree') is not None:
            flags |= 8
        if extra.get('anisotropicElectrostaticEnergyHartree') is not None:
            flags |= 16
        if extra.get('anisotropicXCEnergyHartree') is not None:
            flags |= 32
        if extra.get('dispersionEnergyHartree') is not None:
            flags |= 64
        if extra.get('halogenBondEnergyHartree') is not None:
            flags |= 128
        if extra.get('solvationEnergyHartree') is not None:
            flags |= 256
        if extra.get('fermiLevelAlphaHartree') is not None:
            flags |= 512
        if extra.get('solvationBornEnergyHartree') is not None:
            flags |= 1024
        extra['energyDecompositionFlags'] = flags

    if getattr(result, 'polarizabilityAvailable', False):
        extra['atomicPolarizabilitiesAU'] = tuple(map(float, result.GetAtomicPolarizabilities()))
        extra['atomicC6AU'] = tuple(map(float, result.GetAtomicC6()))
    if getattr(result, 'multipolesAvailable', False):
        extra['atomicDipolesAU'] = tuple(map(tuple, result.GetAtomicDipoles()))
        extra['atomicQuadrupolesAU'] = tuple(map(tuple, result.GetAtomicQuadrupoles()))
        extra['molecularQuadrupoleAU'] = tuple(map(float, result.GetMolecularQuadrupole()))
    if getattr(result, 'nuclearESPAvailable', False):
        extra['nuclearESPAU'] = tuple(map(float, result.GetNuclearESP()))

    from dataclasses import replace
    metadata = dict(base.metadata)
    if extra.get('solvent') is not None:
        metadata['solvent'] = {
            'solvent': extra['solvent'],
            'model': extra.get('solvationModel'),
            'grid': extra.get('solvationGrid'),
        }
    if getattr(result, 'propertiesExportVerified', False):
        metadata['properties_api'] = 1
    return replace(base, metadata=metadata, **extra)


def _unpack_symmetric_3x3(packed: np.ndarray) -> np.ndarray:
    """Unpack [xx, xy, yy, xz, yz, zz] into a symmetric 3x3 matrix."""
    xx, xy, yy, xz, yz, zz = (float(v) for v in packed)
    return np.array(
        [[xx, xy, xz], [xy, yy, yz], [xz, yz, zz]],
        dtype=float,
    )


def _quadrupole_invariants(packed: np.ndarray) -> dict[str, Any]:
    mat = _unpack_symmetric_3x3(packed)
    eigvals = np.linalg.eigvalsh(mat)
    frob = float(np.linalg.norm(mat, ord='fro'))
    anisotropy = float(math.sqrt(1.5) * frob)
    trace = float(np.trace(mat))
    return {
        'eigenvalues_ebohr2': [float(v) for v in eigvals],
        'frobenius_ebohr2': frob,
        'anisotropy_ebohr2': anisotropy,
        'trace_ebohr2': trace,
    }


def _json_safe(payload: dict[str, Any]) -> dict[str, Any]:
    return json.loads(json.dumps(payload, allow_nan=False))


def CalcEnergyDecomposition(result: Any) -> dict[str, Any]:
    """Return xTB energy decomposition components (Hartree) and Fermi levels.

    Method-inapplicable terms (e.g. ``e_aes_ha`` / ``e_axc_ha`` for GFN0/1,
    ``e_xb_ha`` for GFN2, ``g_solv_ha`` in gas phase, or Fermi levels when
    ``electronicTemperature <= 0.1 K``) are returned as ``None`` with explicit
    entries in ``missing``.
    """
    res = _as_property_result(result)
    method = res.method
    missing: dict[str, str] = {}
    available = res.energyDecompositionAvailable

    def _missing_reason(field: str) -> str:
        if not available:
            return 'properties_extension_unavailable'
        if field in ('e_aes_ha', 'e_axc_ha'):
            return 'not_applicable_to_gfn0' if method == 'GFN0-xTB' else 'not_applicable_to_gfn1'
        if field == 'e_xb_ha':
            return 'not_applicable_to_gfn2'
        if field == 'g_solv_ha':
            return 'gas_phase_no_solvation'
        if field in ('g_born_ha', 'g_sasa_ha', 'g_hb_ha', 'g_shift_ha'):
            if res.solvationEnergyHartree is None:
                return 'gas_phase_no_solvation'
            return 'not_separately_returned_by_gfn0'
        if field.startswith('fermi_'):
            return 'fermi_smearing_disabled_low_etemp'
        return 'not_exposed_by_selected_libxtb'

    raw_fields = {
        'e_total_ha': float(res.totalEnergyHartree) if (available and res.totalEnergyHartree is not None) else None,
        'e_elec_ha': res.electronicEnergyHartree if available else None,
        'e_rep_ha': res.repulsionEnergyHartree if available else None,
        'e_es_ha': res.electrostaticEnergyHartree if available else None,
        'e_aes_ha': res.anisotropicElectrostaticEnergyHartree if available else None,
        'e_axc_ha': res.anisotropicXCEnergyHartree if available else None,
        'e_disp_ha': res.dispersionEnergyHartree if available else None,
        'e_xb_ha': res.halogenBondEnergyHartree if available else None,
        'g_solv_ha': res.solvationEnergyHartree if available else None,
        'g_born_ha': res.solvationBornEnergyHartree if available else None,
        'g_sasa_ha': res.solvationSASAEnergyHartree if available else None,
        'g_hb_ha': res.solvationHBEnergyHartree if available else None,
        'g_shift_ha': res.solvationShiftEnergyHartree if available else None,
        'fermi_alpha_ha': res.fermiLevelAlphaHartree if available else None,
        'fermi_beta_ha': res.fermiLevelBetaHartree if available else None,
        'fermi_alpha_ev': res.fermiLevelAlphaEV if available else None,
        'fermi_beta_ev': res.fermiLevelBetaEV if available else None,
    }
    for key, value in raw_fields.items():
        if value is None:
            missing[key] = _missing_reason(key)
        elif not math.isfinite(value):
            raise ValueError(f'nonfinite energy decomposition component: {key}')

    xb_def = (
        'halogen_bond_correction'
        if method == 'GFN1-xTB'
        else 'short_range_bond_correction'
        if method == 'GFN0-xTB'
        else None
    )
    payload = {
        'available': available,
        'method': method,
        'solvent': res.solvent,
        'solvation_model': res.solvationModel,
        **raw_fields,
        'xb_term_definition': xb_def,
        'units': {
            'e_total_ha': 'Hartree',
            'e_elec_ha': 'Hartree',
            'e_rep_ha': 'Hartree',
            'e_es_ha': 'Hartree',
            'e_aes_ha': 'Hartree',
            'e_axc_ha': 'Hartree',
            'e_disp_ha': 'Hartree',
            'e_xb_ha': 'Hartree',
            'g_solv_ha': 'Hartree',
            'g_born_ha': 'Hartree',
            'g_sasa_ha': 'Hartree',
            'g_hb_ha': 'Hartree',
            'g_shift_ha': 'Hartree',
            'fermi_alpha_ha': 'Hartree',
            'fermi_beta_ha': 'Hartree',
            'fermi_alpha_ev': 'eV',
            'fermi_beta_ev': 'eV',
        },
        'invariance': 'all energy components and Fermi levels are rigid-motion and atom-permutation invariant',
        'missing': missing,
    }
    return _json_safe(payload)


def CalcMultipoleDescriptors(result: Any, *, mol: Any = None) -> dict[str, Any]:
    """Return atomic and molecular multipole descriptors with rotational invariants.

    GFN2 exposes self-consistent CAMM dipoles/quadrupoles; GFN1 exposes post-SCF
    Mulliken CAMM dipoles/quadrupoles. GFN0 does not compute atomic multipoles
    and returns ``available=False`` with ``None`` fields and explicit ``missing``
    reasons.
    """
    res = _as_property_result(result)
    if mol is not None:
        _validate_mapping(mol, res)

    multipole_units = {
        'dip_ebohr': 'e*Bohr',
        'dip_norm_ebohr': 'e*Bohr',
        'dip_norm_debye': 'Debye',
        'quad_ebohr2': 'e*Bohr^2',
        'quad_eigenvalues_ebohr2': 'e*Bohr^2',
        'quad_frobenius_ebohr2': 'e*Bohr^2',
        'quad_anisotropy_ebohr2': 'e*Bohr^2',
        'quad_trace_ebohr2': 'e*Bohr^2',
        'mol_quadrupole_ebohr2': 'e*Bohr^2',
        'mol_quadrupole_eigenvalues_ebohr2': 'e*Bohr^2',
        'mol_quadrupole_frobenius_ebohr2': 'e*Bohr^2',
        'mol_quadrupole_anisotropy_ebohr2': 'e*Bohr^2',
        'mol_quadrupole_trace_ebohr2': 'e*Bohr^2',
        'atom_dip_norm_max_ebohr': 'e*Bohr',
        'atom_dip_norm_mean_ebohr': 'e*Bohr',
        'atom_quad_frobenius_max_ebohr2': 'e*Bohr^2',
        'atom_quad_frobenius_mean_ebohr2': 'e*Bohr^2',
    }
    if not res.multipolesAvailable:
        reason = (
            'not_computed_by_gfn0'
            if res.method == 'GFN0-xTB'
            else 'properties_extension_unavailable'
        )
        missing = {
            'multipole_model': reason,
            'mol_quadrupole_ebohr2': reason,
            'mol_quadrupole_eigenvalues_ebohr2': reason,
            'mol_quadrupole_frobenius_ebohr2': reason,
            'mol_quadrupole_anisotropy_ebohr2': reason,
            'mol_quadrupole_trace_ebohr2': reason,
            'atom_dip_norm_max_ebohr': reason,
            'atom_dip_norm_mean_ebohr': reason,
            'atom_quad_frobenius_max_ebohr2': reason,
            'atom_quad_frobenius_mean_ebohr2': reason,
            'atoms': reason,
        }
        return _json_safe({
            'available': False,
            'method': res.method,
            'multipole_model': None,
            'component_order': list(_QUADRUPOLE_COMPONENTS),
            'mol_quadrupole_ebohr2': None,
            'mol_quadrupole_eigenvalues_ebohr2': None,
            'mol_quadrupole_frobenius_ebohr2': None,
            'mol_quadrupole_anisotropy_ebohr2': None,
            'mol_quadrupole_trace_ebohr2': None,
            'atom_dip_norm_max_ebohr': None,
            'atom_dip_norm_mean_ebohr': None,
            'atom_quad_frobenius_max_ebohr2': None,
            'atom_quad_frobenius_mean_ebohr2': None,
            'atoms': [],
            'units': multipole_units,
            'invariance': (
                'atomic dipole norms, quadrupole eigenvalues, Frobenius norms and anisotropies '
                'are rigid-translation and rotation invariant; molecular quadrupole depends on '
                'origin for systems with nonzero charge or dipole'
            ),
            'missing': missing,
        })

    dip_arr = res.GetAtomicDipoles()
    qp_arr = res.GetAtomicQuadrupoles()
    mol_qp = res.GetMolecularQuadrupole()
    if dip_arr.ndim != 2 or dip_arr.shape[1] != 3 or not np.isfinite(dip_arr).all():
        raise ValueError('atomic dipoles must be finite [natoms, 3]')
    if qp_arr.ndim != 2 or qp_arr.shape != (dip_arr.shape[0], 6) or not np.isfinite(qp_arr).all():
        raise ValueError('atomic quadrupoles must be finite [natoms, 6]')
    if mol_qp.shape != (6,) or not np.isfinite(mol_qp).all():
        raise ValueError('molecular quadrupole must be finite [6]')

    mol_inv = _quadrupole_invariants(mol_qp)
    atom_rows = []
    dip_norms = []
    quad_frobs = []
    for idx in range(dip_arr.shape[0]):
        dip_vec = dip_arr[idx]
        dnorm = float(np.linalg.norm(dip_vec))
        qinv = _quadrupole_invariants(qp_arr[idx])
        dip_norms.append(dnorm)
        quad_frobs.append(qinv['frobenius_ebohr2'])
        atom_rows.append({
            'atom_index': idx,
            'dip_ebohr': [float(v) for v in dip_vec],
            'dip_norm_ebohr': dnorm,
            'dip_norm_debye': dnorm * AU_DIPOLE_TO_DEBYE,
            'quad_ebohr2': [float(v) for v in qp_arr[idx]],
            'quad_eigenvalues_ebohr2': qinv['eigenvalues_ebohr2'],
            'quad_frobenius_ebohr2': qinv['frobenius_ebohr2'],
            'quad_anisotropy_ebohr2': qinv['anisotropy_ebohr2'],
            'quad_trace_ebohr2': qinv['trace_ebohr2'],
        })

    model = 'gfn2_scc_camm' if res.method == 'GFN2-xTB' else 'gfn1_post_scf_camm'
    return _json_safe({
        'available': True,
        'method': res.method,
        'multipole_model': model,
        'component_order': list(_QUADRUPOLE_COMPONENTS),
        'mol_quadrupole_ebohr2': [float(v) for v in mol_qp],
        'mol_quadrupole_eigenvalues_ebohr2': mol_inv['eigenvalues_ebohr2'],
        'mol_quadrupole_frobenius_ebohr2': mol_inv['frobenius_ebohr2'],
        'mol_quadrupole_anisotropy_ebohr2': mol_inv['anisotropy_ebohr2'],
        'mol_quadrupole_trace_ebohr2': mol_inv['trace_ebohr2'],
        'atom_dip_norm_max_ebohr': float(max(dip_norms)),
        'atom_dip_norm_mean_ebohr': float(np.mean(dip_norms)),
        'atom_quad_frobenius_max_ebohr2': float(max(quad_frobs)),
        'atom_quad_frobenius_mean_ebohr2': float(np.mean(quad_frobs)),
        'atoms': atom_rows,
        'units': multipole_units,
        'invariance': (
            'atomic dipole norms, quadrupole eigenvalues, Frobenius norms and anisotropies '
            'are rigid-translation and rotation invariant; Cartesian vector/tensor components '
            'rotate with the coordinate frame; molecular quadrupole depends on origin for '
            'systems with nonzero charge or dipole'
        ),
        'missing': {},
    })


def CalcNuclearESP(result: Any, *, mol: Any = None) -> dict[str, Any]:
    """Return exact electrostatic potential at nuclei (EPN, V_N) in Hartree/e and kcal/mol.

    Evaluates the valence electronic plus other-nuclei electrostatic potential at
    each nuclear position in original RDKit atom order, omitting the singular
    self-nuclear term ``Z_A / |r - R_A|``.
    """
    res = _as_property_result(result)
    if mol is not None:
        _validate_mapping(mol, res)

    epn_units = {
        'epn_au': 'Hartree/e',
        'epn_kcalmol': 'kcal/(mol*e)',
        'epn_min_au': 'Hartree/e',
        'epn_max_au': 'Hartree/e',
        'epn_mean_au': 'Hartree/e',
        'epn_sd_au': 'Hartree/e',
        'epn_min_kcalmol': 'kcal/(mol*e)',
        'epn_max_kcalmol': 'kcal/(mol*e)',
        'epn_mean_kcalmol': 'kcal/(mol*e)',
        'epn_sd_kcalmol': 'kcal/(mol*e)',
    }
    if not res.nuclearESPAvailable:
        reason = 'properties_extension_unavailable'
        return _json_safe({
            'available': False,
            'method': res.method,
            'epn_min_au': None,
            'epn_max_au': None,
            'epn_mean_au': None,
            'epn_sd_au': None,
            'epn_min_kcalmol': None,
            'epn_max_kcalmol': None,
            'epn_mean_kcalmol': None,
            'epn_sd_kcalmol': None,
            'atoms': [],
            'definition': (
                'Exact valence electronic plus other-nuclei electrostatic potential at each '
                'nuclear position omitting the singular self-nuclear term Z_A / r_iA'
            ),
            'units': epn_units,
            'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices',
            'missing': {
                'epn_min_au': reason,
                'epn_max_au': reason,
                'epn_mean_au': reason,
                'epn_sd_au': reason,
                'epn_min_kcalmol': reason,
                'epn_max_kcalmol': reason,
                'epn_mean_kcalmol': reason,
                'epn_sd_kcalmol': reason,
                'atoms': reason,
            },
        })

    values = res.GetNuclearESP()
    if values.ndim != 1 or not len(values) or not np.isfinite(values).all():
        raise ValueError('nuclear ESP values must be a nonempty finite 1D array')

    kcal = values * HARTREE_TO_KCAL_MOL
    rows = [
        {
            'atom_index': idx,
            'epn_au': float(values[idx]),
            'epn_kcalmol': float(kcal[idx]),
        }
        for idx in range(len(values))
    ]
    return _json_safe({
        'available': True,
        'method': res.method,
        'epn_min_au': float(np.min(values)),
        'epn_max_au': float(np.max(values)),
        'epn_mean_au': float(np.mean(values)),
        'epn_sd_au': float(np.std(values, ddof=0)),
        'epn_min_kcalmol': float(np.min(kcal)),
        'epn_max_kcalmol': float(np.max(kcal)),
        'epn_mean_kcalmol': float(np.mean(kcal)),
        'epn_sd_kcalmol': float(np.std(kcal, ddof=0)),
        'atoms': rows,
        'definition': (
            'Exact valence electronic plus other-nuclei electrostatic potential at each '
            'nuclear position omitting the singular self-nuclear term Z_A / r_iA'
        ),
        'units': epn_units,
        'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices',
        'missing': {},
    })


def CalcPolarizabilityDescriptors(result: Any, *, mol: Any = None) -> dict[str, Any]:
    """Return D4 static dipole polarizabilities and dispersion C6/C8 descriptors.

    Available for GFN2-xTB (self-consistent D4) and GFN0-xTB (D4 with EN
    charges). GFN1-xTB uses DFT-D3 dispersion without D4 polarizabilities and
    returns ``available=False`` with ``None`` fields and explicit ``missing``
    reasons.
    """
    res = _as_property_result(result)
    if mol is not None:
        _validate_mapping(mol, res)

    ff_note = (
        'Finite-field electronic polarizability is not supported by libxtb 6.7.1 GFN0/1/2 C API '
        '(external electric fields are restricted to tblite/GFN-FF/PTB, and xtb_setExternalCharges '
        'couples only to monopole shell charges rather than the electronic dipole operator). '
        'Use analytical D4 mol_pol_bohr3 and atomic alpha_bohr3 instead.'
    )

    if not res.polarizabilityAvailable:
        reason = (
            'not_computed_by_gfn1_d3_dispersion'
            if res.method == 'GFN1-xTB'
            else 'properties_extension_unavailable'
        )
        return _json_safe({
            'available': False,
            'method': res.method,
            'polarizability_model': None,
            'mol_pol_bohr3': None,
            'mol_pol_angstrom3': None,
            'mol_c6_ha_bohr6': None,
            'mol_c8_ha_bohr8': None,
            'atoms': [],
            'finite_field_supported': False,
            'finite_field_note': ff_note,
            'units': {
                'mol_pol_bohr3': 'Bohr^3',
                'mol_pol_angstrom3': 'Angstrom^3',
                'mol_c6_ha_bohr6': 'Hartree*Bohr^6',
                'mol_c8_ha_bohr8': 'Hartree*Bohr^8',
                'alpha_bohr3': 'Bohr^3',
                'alpha_angstrom3': 'Angstrom^3',
                'c6_aa_ha_bohr6': 'Hartree*Bohr^6',
            },
            'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices',
            'missing': {
                'polarizability_model': reason,
                'mol_pol_bohr3': reason,
                'mol_pol_angstrom3': reason,
                'mol_c6_ha_bohr6': reason,
                'mol_c8_ha_bohr8': reason,
                'atoms': reason,
            },
        })

    alpha = res.GetAtomicPolarizabilities()
    c6 = res.GetAtomicC6()
    if alpha.ndim != 1 or c6.shape != alpha.shape or not (np.isfinite(alpha).all() and np.isfinite(c6).all()):
        raise ValueError('atomic polarizabilities and C6 must be finite matching 1D arrays')
    molpol = float(res.molecularPolarizabilityAU)
    molc6 = float(res.molecularC6AU) if res.molecularC6AU is not None else None
    molc8 = float(res.molecularC8AU) if res.molecularC8AU is not None else None
    missing: dict[str, str] = {}
    if molc6 is None:
        missing['mol_c6_ha_bohr6'] = 'not_exposed_by_selected_libxtb'
    if molc8 is None:
        missing['mol_c8_ha_bohr8'] = 'not_exposed_by_selected_libxtb'

    rows = [
        {
            'atom_index': idx,
            'alpha_bohr3': float(alpha[idx]),
            'alpha_angstrom3': float(alpha[idx] * BOHR3_TO_ANGSTROM3),
            'c6_aa_ha_bohr6': float(c6[idx]),
        }
        for idx in range(len(alpha))
    ]
    model = 'gfn2_d4_self_consistent' if res.method == 'GFN2-xTB' else 'gfn0_d4_en_charges'
    return _json_safe({
        'available': True,
        'method': res.method,
        'polarizability_model': model,
        'mol_pol_bohr3': molpol,
        'mol_pol_angstrom3': molpol * BOHR3_TO_ANGSTROM3,
        'mol_c6_ha_bohr6': molc6,
        'mol_c8_ha_bohr8': molc8,
        'atoms': rows,
        'finite_field_supported': False,
        'finite_field_note': ff_note,
        'units': {
            'mol_pol_bohr3': 'Bohr^3',
            'mol_pol_angstrom3': 'Angstrom^3',
            'mol_c6_ha_bohr6': 'Hartree*Bohr^6',
            'mol_c8_ha_bohr8': 'Hartree*Bohr^8',
            'alpha_bohr3': 'Bohr^3',
            'alpha_angstrom3': 'Angstrom^3',
            'c6_aa_ha_bohr6': 'Hartree*Bohr^6',
        },
        'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices',
        'missing': missing,
    })


def CalcFiniteFieldPolarizability(mol: Any, *, options: Any = None, fieldAU: float = 0.001) -> dict[str, Any]:
    """Raise DescriptorUnavailableError explaining why libxtb 6.7.1 lacks GFN0/1/2 finite-field response."""
    raise DescriptorUnavailableError(
        'Finite-field polarizability is unavailable in libxtb 6.7.1 for GFN0/1/2-xTB '
        '(external electric fields are not supported by the C API without tblite/PTB, and '
        'point-charge embedding couples only to monopole shell charges, missing intra-atomic '
        'and out-of-plane electronic polarization). Use CalcPolarizabilityDescriptors(result) '
        'for analytical D4 molecular and atomic static polarizabilities.'
    )
