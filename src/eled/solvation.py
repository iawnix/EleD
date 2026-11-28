"""Implicit solvation (ALPB/GBSA) descriptors at fixed molecular geometry.

Computes matched gas-phase and solvated single points with identical method,
geometry, charge, spin, and numerical settings to quantify vertical solvation
free energy, dipole enhancement, frontier-orbital shifts, and atomic charge
redistribution.
"""
from __future__ import annotations

from dataclasses import replace
import json
from typing import Any

import numpy as np

from .calculator import XTBCalculator
from .electronic_descriptors import AU_DIPOLE_TO_DEBYE
from .options import XTBOptions, _validate_solvation_settings

HARTREE_TO_KCAL_MOL = 627.5094740631


def CalcSolvationDescriptors(
    mol: Any,
    solvent: str,
    *,
    confId: int = -1,
    method: str = 'GFN2-xTB',
    solvationModel: str = 'alpb',
    solvationGrid: int | None = None,
    charge: int | None = None,
    numUnpairedElectrons: int | None = None,
    options: XTBOptions | None = None,
    xtbLibrary: Any = None,
    maxIterations: int = 250,
    accuracy: float = 1.0,
    electronicTemperature: float = 300.0,
) -> dict[str, Any]:
    """Compute matched gas-phase vs implicit-solvation descriptors at fixed geometry.

    Parameters
    ----------
    mol : RDKit Mol
        Explicit-hydrogen 3D RDKit molecule (never modified).
    solvent : str
        xTB solvent name (e.g. ``'water'``, ``'dmso'``, ``'acetonitrile'``, ``'chloroform'``).
    solvationModel : {'alpb', 'gbsa'}
        Implicit solvation model (default ``'alpb'``).
    solvationGrid : int | None
        Optional Lebedev angular grid size (default ``None`` -> xTB default ``230``).
    options : XTBOptions | None
        Optional base calculation options; ``solvent``, ``solvationModel``, and
        ``solvationGrid`` override any solvation fields on ``options``.
    """
    if solvent is None:
        raise ValueError('solvent must be a nonempty solvent name string')
    solv_name, solv_model, solv_grid = _validate_solvation_settings(
        solvent, solvationModel, solvationGrid
    )
    if options is not None:
        if not isinstance(options, XTBOptions):
            raise TypeError('options must be an XTBOptions instance or None')
        gas_opts = replace(options, solvent=None, solvationModel='alpb', solvationGrid=None)
        solv_opts = replace(
            options, solvent=solv_name, solvationModel=solv_model, solvationGrid=solv_grid
        )
    else:
        gas_opts = XTBOptions(
            confId=confId,
            method=method,
            charge=charge,
            numUnpairedElectrons=numUnpairedElectrons,
            xtbLibrary=xtbLibrary,
            maxIterations=maxIterations,
            accuracy=accuracy,
            electronicTemperature=electronicTemperature,
            solvent=None,
        )
        solv_opts = replace(
            gas_opts, solvent=solv_name, solvationModel=solv_model, solvationGrid=solv_grid
        )

    gas_res = XTBCalculator(**gas_opts.ToDict()).Calculate(mol)
    solv_res = XTBCalculator(**solv_opts.ToDict()).Calculate(mol)

    e_gas = float(gas_res.totalEnergyHartree)
    e_solv = float(solv_res.totalEnergyHartree)
    delta_g_ha = e_solv - e_gas
    delta_g_kcal = delta_g_ha * HARTREE_TO_KCAL_MOL

    missing: dict[str, str] = {}
    g_direct = solv_res.solvationEnergyHartree
    if g_direct is None:
        missing['g_solv_direct_ha'] = 'properties_extension_unavailable'
        missing['g_solv_direct_kcalmol'] = 'properties_extension_unavailable'
        missing['solute_polarization_energy_ha'] = 'properties_extension_unavailable'
        missing['solute_polarization_energy_kcalmol'] = 'properties_extension_unavailable'
        g_direct_kcal = None
        pol_energy_ha = None
        pol_energy_kcal = None
    else:
        g_direct_kcal = float(g_direct * HARTREE_TO_KCAL_MOL)
        pol_energy_ha = float(delta_g_ha - g_direct)
        pol_energy_kcal = float(pol_energy_ha * HARTREE_TO_KCAL_MOL)

    sub_terms = {
        'g_born_ha': solv_res.solvationBornEnergyHartree,
        'g_sasa_ha': solv_res.solvationSASAEnergyHartree,
        'g_hb_ha': solv_res.solvationHBEnergyHartree,
        'g_shift_ha': solv_res.solvationShiftEnergyHartree,
    }
    for key, value in sub_terms.items():
        if value is None:
            missing[key] = (
                'not_separately_returned_by_gfn0'
                if solv_res.method == 'GFN0-xTB'
                else 'properties_extension_unavailable'
            )

    dip_gas = gas_res.GetDipole()
    dip_solv = solv_res.GetDipole()
    dip_gas_debye = float(np.linalg.norm(dip_gas) * AU_DIPOLE_TO_DEBYE)
    dip_solv_debye = float(np.linalg.norm(dip_solv) * AU_DIPOLE_TO_DEBYE)
    delta_dip_debye = float(dip_solv_debye - dip_gas_debye)
    delta_dip_vec_debye = float(np.linalg.norm(dip_solv - dip_gas) * AU_DIPOLE_TO_DEBYE)

    gas_homo, solv_homo = gas_res.homoEV, solv_res.homoEV
    gas_lumo, solv_lumo = gas_res.lumoEV, solv_res.lumoEV
    gas_gap, solv_gap = gas_res.gapEV, solv_res.gapEV

    def _diff(a: float | None, b: float | None, name: str, reason: str) -> float | None:
        if a is None or b is None:
            missing[name] = reason
            return None
        return float(b - a)

    if gas_homo is None:
        missing['gas_homo_ev'] = 'no_occupied_orbital_under_threshold'
    if solv_homo is None:
        missing['solv_homo_ev'] = 'no_occupied_orbital_under_threshold'
    if gas_lumo is None:
        missing['gas_lumo_ev'] = 'no_virtual_orbital_under_threshold'
    if solv_lumo is None:
        missing['solv_lumo_ev'] = 'no_virtual_orbital_under_threshold'
    if gas_gap is None:
        missing['gas_gap_ev'] = 'frontier_orbital_unavailable'
    if solv_gap is None:
        missing['solv_gap_ev'] = 'frontier_orbital_unavailable'

    delta_homo = _diff(gas_homo, solv_homo, 'delta_homo_ev', 'frontier_orbital_unavailable')
    delta_lumo = _diff(gas_lumo, solv_lumo, 'delta_lumo_ev', 'frontier_orbital_unavailable')
    delta_gap = _diff(gas_gap, solv_gap, 'delta_gap_ev', 'frontier_orbital_unavailable')

    q_gas = gas_res.GetAtomicCharges()
    q_solv = solv_res.GetAtomicCharges()
    dq = q_solv - q_gas
    has_cm5 = bool(len(gas_res.cm5Charges) and len(solv_res.cm5Charges))
    cm5_gas = gas_res.GetCM5Charges() if has_cm5 else None
    cm5_solv = solv_res.GetCM5Charges() if has_cm5 else None

    cm5_reason = (
        'cm5_undefined_for_gfn0_eeq'
        if solv_res.method == 'GFN0-xTB'
        else 'cm5_charges_unavailable'
    )
    if not has_cm5:
        missing['cm5_charges'] = cm5_reason

    atom_rows = []
    for idx in range(len(q_gas)):
        row: dict[str, Any] = {
            'atom_index': idx,
            'q_gas_e': float(q_gas[idx]),
            'q_solv_e': float(q_solv[idx]),
            'solv_delta_q_e': float(dq[idx]),
        }
        if has_cm5:
            row['cm5_q_gas_e'] = float(cm5_gas[idx])
            row['cm5_q_solv_e'] = float(cm5_solv[idx])
            row['solv_delta_cm5_q_e'] = float(cm5_solv[idx] - cm5_gas[idx])
            row['missing'] = {}
        else:
            row['cm5_q_gas_e'] = None
            row['cm5_q_solv_e'] = None
            row['solv_delta_cm5_q_e'] = None
            row['missing'] = {
                'cm5_q_gas_e': cm5_reason,
                'cm5_q_solv_e': cm5_reason,
                'solv_delta_cm5_q_e': cm5_reason,
            }
        atom_rows.append(row)

    payload = {
        'method': solv_res.method,
        'solvent': solv_res.solvent,
        'solvation_model': solv_res.solvationModel,
        'solvation_grid': solv_res.solvationGrid,
        'charge_model': solv_res.chargeModel,
        'e_gas_ha': e_gas,
        'e_solv_ha': e_solv,
        'delta_g_solv_ha': delta_g_ha,
        'delta_g_solv_kcalmol': delta_g_kcal,
        'g_solv_direct_ha': g_direct,
        'g_solv_direct_kcalmol': g_direct_kcal,
        'solute_polarization_energy_ha': pol_energy_ha,
        'solute_polarization_energy_kcalmol': pol_energy_kcal,
        **sub_terms,
        'gas_dipole_debye': dip_gas_debye,
        'solv_dipole_debye': dip_solv_debye,
        'delta_dipole_debye': delta_dip_debye,
        'delta_dipole_vec_norm_debye': delta_dip_vec_debye,
        'gas_homo_ev': gas_homo,
        'solv_homo_ev': solv_homo,
        'delta_homo_ev': delta_homo,
        'gas_lumo_ev': gas_lumo,
        'solv_lumo_ev': solv_lumo,
        'delta_lumo_ev': delta_lumo,
        'gas_gap_ev': gas_gap,
        'solv_gap_ev': solv_gap,
        'delta_gap_ev': delta_gap,
        'charge_shift_rms_e': float(np.sqrt(np.mean(dq ** 2))),
        'charge_shift_max_abs_e': float(np.max(np.abs(dq))),
        'charge_shift_l1_e': float(np.sum(np.abs(dq))),
        'atoms': atom_rows,
        'provenance': {
            'conformer_id': solv_res.conformerId,
            'charge': solv_res.charge,
            'num_unpaired_electrons': solv_res.numUnpairedElectrons,
            'geometry_sha256': solv_res.metadata.get('geometry_sha256'),
            'ordered_structure_sha256': solv_res.metadata.get('ordered_structure_sha256'),
            'calculation_settings': solv_res.metadata.get('calculation_settings'),
        },
        'units': {
            'e_gas_ha': 'Hartree',
            'e_solv_ha': 'Hartree',
            'delta_g_solv_ha': 'Hartree',
            'delta_g_solv_kcalmol': 'kcal/mol',
            'g_solv_direct_ha': 'Hartree',
            'g_solv_direct_kcalmol': 'kcal/mol',
            'solute_polarization_energy_ha': 'Hartree',
            'solute_polarization_energy_kcalmol': 'kcal/mol',
            'g_born_ha': 'Hartree',
            'g_sasa_ha': 'Hartree',
            'g_hb_ha': 'Hartree',
            'g_shift_ha': 'Hartree',
            'gas_dipole_debye': 'Debye',
            'solv_dipole_debye': 'Debye',
            'delta_dipole_debye': 'Debye',
            'delta_dipole_vec_norm_debye': 'Debye',
            'gas_homo_ev': 'eV',
            'solv_homo_ev': 'eV',
            'delta_homo_ev': 'eV',
            'gas_lumo_ev': 'eV',
            'solv_lumo_ev': 'eV',
            'delta_lumo_ev': 'eV',
            'gas_gap_ev': 'eV',
            'solv_gap_ev': 'eV',
            'delta_gap_ev': 'eV',
            'charge_shift_rms_e': 'e',
            'charge_shift_max_abs_e': 'e',
            'charge_shift_l1_e': 'e',
            'q_gas_e': 'e',
            'q_solv_e': 'e',
            'solv_delta_q_e': 'e',
            'cm5_q_gas_e': 'e',
            'cm5_q_solv_e': 'e',
            'solv_delta_cm5_q_e': 'e',
        },
        'invariance': (
            'rigid-translation invariant; rotationally invariant up to Lebedev surface-grid '
            'quadrature discretization; atom-permutation invariant with remapped indices; '
            'dipole magnitude is origin dependent for charged systems'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))
