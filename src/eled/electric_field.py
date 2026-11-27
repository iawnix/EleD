"""Bond-aligned local electric field descriptors from in-memory libxtb ESP.

Evaluates the molecular electric field ``E(r_mid) = -grad V(r_mid)`` at every
RDKit bond midpoint ``r_mid = (R_begin + R_end) / 2`` by 6-point central finite
differences of the exact in-memory ``libxtb`` electrostatic potential
``V(r)`` (Hartree/e) in **one** batched single-point calculation (``6 * n_bonds``
stencil points).

Sign and orientation convention:
- ``E = -grad V`` points from higher electrostatic potential (positive charge)
  toward lower electrostatic potential (negative charge), i.e. the force per
  unit positive test charge.
- ``u = (R_end - R_begin) / ||R_end - R_begin||`` follows RDKit's bond
  orientation (``begin_atom_index -> end_atom_index``).
- ``e_parallel_au = E . u`` is signed and antisymmetrizes upon bond reversal.
- ``e_perp_au = ||E - (E . u) u||`` and ``e_norm_au = ||E||`` are nonnegative
  and invariant to rigid 3D motion and bond direction reversal.
- Conversion: ``1 a.u. = E_h / (e * a_0) = 5142.2067476326 MV/cm``.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from numbers import Real
from typing import Any

import numpy as np

from .calculator import XTBCalculator
from .electronic_descriptors import _as_result
from .errors import XTBEnvironmentError, XTBRunError
from .esp import ANGSTROM_TO_BOHR, CalcESP
from .input import _effective_charge_uhf, _xyz_from_mol
from .mapping import _validate_mapping
from .options import XTBOptions

# CODATA 2018 exact atomic unit of electric field in MV/cm (1 a.u. = 5.1422067476326e11 V/m).
AU_EFIELD_TO_MV_CM = 5142.2067476326
BOND_EFIELD_PROTOCOL = 'eled-bond-electric-field-v1'

__all__ = [
    'AU_EFIELD_TO_MV_CM',
    'CalcBondElectricFieldFromESP',
    'CalcBondElectricField',
]


def _validate_step(stepAngstrom: Any) -> float:
    if (
        isinstance(stepAngstrom, (bool, np.bool_))
        or not isinstance(stepAngstrom, Real)
        or not math.isfinite(stepAngstrom)
        or stepAngstrom <= 0.0
    ):
        raise ValueError('stepAngstrom must be finite and positive')
    step = float(stepAngstrom)
    if step > 0.1:
        raise ValueError('stepAngstrom must be <= 0.1 Angstrom to stay well away from nuclear singularities')
    return step


def _bond_stencil(mol: Any, conf_id: int, step_angstrom: float) -> tuple[list[dict[str, Any]], np.ndarray]:
    """Build bond metadata and [6 * n_bonds, 3] central-difference stencil points (Angstrom)."""
    coords = np.asarray(mol.GetConformer(conf_id).GetPositions(), dtype=float)
    bonds_meta: list[dict[str, Any]] = []
    stencil_rows: list[np.ndarray] = []
    axes = np.eye(3, dtype=float) * step_angstrom

    for bond in mol.GetBonds():
        b_idx = int(bond.GetIdx())
        begin = int(bond.GetBeginAtomIdx())
        end = int(bond.GetEndAtomIdx())
        r_begin = coords[begin]
        r_end = coords[end]
        vec = r_end - r_begin
        length = float(np.linalg.norm(vec))
        if not math.isfinite(length) or length <= 1e-6:
            raise ValueError(f'bond {b_idx} ({begin}-{end}) has degenerate length {length}')
        unit = vec / length
        midpoint = 0.5 * (r_begin + r_end)
        if step_angstrom >= 0.25 * length:
            raise ValueError(
                f'stepAngstrom ({step_angstrom}) is too large relative to bond {b_idx} length ({length:.4f} A)'
            )
        bonds_meta.append({
            'bond_index': b_idx,
            'begin_atom_index': begin,
            'end_atom_index': end,
            'bond_length_a': length,
            'midpoint_angstrom': midpoint,
            'unit_vector': unit,
        })
        for axis in range(3):
            stencil_rows.append(midpoint + axes[axis])
            stencil_rows.append(midpoint - axes[axis])

    if not stencil_rows:
        return bonds_meta, np.empty((0, 3), dtype=float)
    return bonds_meta, np.asarray(stencil_rows, dtype=float)


def _field_units() -> dict[str, str]:
    return {
        'bond_length_a': 'angstrom',
        'midpoint_angstrom': 'angstrom',
        'e_vec_au': 'Hartree/(e*Bohr)',
        'e_parallel_au': 'Hartree/(e*Bohr)',
        'e_parallel_mv_cm': 'MV/cm',
        'e_perp_au': 'Hartree/(e*Bohr)',
        'e_perp_mv_cm': 'MV/cm',
        'e_norm_au': 'Hartree/(e*Bohr)',
        'e_norm_mv_cm': 'MV/cm',
        'step_angstrom': 'angstrom',
    }


def CalcBondElectricFieldFromESP(
    mol: Any,
    result: Any,
    *,
    stepAngstrom: float = 1e-3,
) -> dict[str, Any]:
    """Compute bond-midpoint electric fields from an existing ESP-sampled result.

    Requires ``result`` to contain ESP samples at the exact ``6 * n_bonds``
    central-difference stencil coordinates generated for ``(mol, confId,
    stepAngstrom)``. Never invokes xTB.
    """
    step = _validate_step(stepAngstrom)
    res = _as_result(result)
    _validate_mapping(mol, res)
    conf_id = int(res.metadata['conformer_id'])
    bonds_meta, expected_points = _bond_stencil(mol, conf_id, step)

    if not bonds_meta:
        payload = {
            'n_bonds': 0,
            'bonds': [],
            'method': res.method,
            'backend': res.backend,
            'step_angstrom': step,
            'units': _field_units(),
            'analysis_settings': {
                'protocol': BOND_EFIELD_PROTOCOL,
                'stencil': '6_point_central_finite_difference',
                'step_angstrom': step,
                'au_to_mv_cm': AU_EFIELD_TO_MV_CM,
                'sign_convention': 'E = -grad(V); e_parallel = E . (R_end - R_begin)/|R_end - R_begin|',
                'potential_source': 'libxtb_total_esp_valence_plus_nuclei',
            },
            'invariance': (
                'e_norm and e_perp are rigid-motion and bond-reversal invariant; '
                'e_parallel is rigid-motion invariant and antisymmetrizes under bond reversal; '
                'e_vec rotates with the Cartesian frame; preserves RDKit bond order'
            ),
            'missing': {'bonds': 'molecule_has_no_bonds'},
        }
        return json.loads(json.dumps(payload, allow_nan=False))

    esp_vals = np.asarray(CalcESP(res), dtype=float)
    esp_pts = np.asarray(res.GetESPPoints(), dtype=float)
    if esp_pts.shape != expected_points.shape or not np.allclose(esp_pts, expected_points, atol=1e-10, rtol=0.0):
        raise ValueError(
            'result ESP sampling points do not match the expected 6*n_bonds central-difference '
            'stencil for (mol, confId, stepAngstrom); use CalcBondElectricField(mol, ...)'
        )

    two_h_bohr = 2.0 * step * ANGSTROM_TO_BOHR
    rows: list[dict[str, Any]] = []
    for idx, meta in enumerate(bonds_meta):
        base = 6 * idx
        ex = -(esp_vals[base + 0] - esp_vals[base + 1]) / two_h_bohr
        ey = -(esp_vals[base + 2] - esp_vals[base + 3]) / two_h_bohr
        ez = -(esp_vals[base + 4] - esp_vals[base + 5]) / two_h_bohr
        e_vec = np.array([ex, ey, ez], dtype=float)
        if not np.isfinite(e_vec).all():
            raise ValueError(f'nonfinite electric field on bond {meta["bond_index"]}')
        unit = meta['unit_vector']
        e_par = float(np.dot(e_vec, unit))
        e_perp_vec = e_vec - e_par * unit
        e_perp = float(np.linalg.norm(e_perp_vec))
        e_norm = float(np.linalg.norm(e_vec))
        rows.append({
            'bond_index': meta['bond_index'],
            'begin_atom_index': meta['begin_atom_index'],
            'end_atom_index': meta['end_atom_index'],
            'bond_length_a': meta['bond_length_a'],
            'midpoint_angstrom': [float(v) for v in meta['midpoint_angstrom']],
            'e_vec_au': [float(v) for v in e_vec],
            'e_parallel_au': e_par,
            'e_parallel_mv_cm': float(e_par * AU_EFIELD_TO_MV_CM),
            'e_perp_au': e_perp,
            'e_perp_mv_cm': float(e_perp * AU_EFIELD_TO_MV_CM),
            'e_norm_au': e_norm,
            'e_norm_mv_cm': float(e_norm * AU_EFIELD_TO_MV_CM),
            'missing': {},
        })

    payload = {
        'n_bonds': len(rows),
        'bonds': rows,
        'method': res.method,
        'backend': res.backend,
        'step_angstrom': step,
        'units': _field_units(),
        'analysis_settings': {
            'protocol': BOND_EFIELD_PROTOCOL,
            'stencil': '6_point_central_finite_difference',
            'step_angstrom': step,
            'au_to_mv_cm': AU_EFIELD_TO_MV_CM,
            'sign_convention': 'E = -grad(V); e_parallel = E . (R_end - R_begin)/|R_end - R_begin|',
            'potential_source': 'libxtb_total_esp_valence_plus_nuclei',
            'conformer_id': conf_id,
            'stencil_point_count': int(len(expected_points)),
        },
        'invariance': (
            'e_norm and e_perp are rigid-motion and bond-reversal invariant; '
            'e_parallel is rigid-motion invariant and antisymmetrizes under bond reversal; '
            'e_vec rotates with the Cartesian frame; preserves RDKit bond order'
        ),
        'missing': {},
    }
    return json.loads(json.dumps(payload, allow_nan=False))


def CalcBondElectricField(
    mol: Any,
    *,
    stepAngstrom: float = 1e-3,
    options: XTBOptions | None = None,
    interface: str = 'native',
) -> dict[str, Any]:
    """Evaluate bond-aligned local electric fields in one batched xTB single point.

    Runs a single xTB calculation with ``6 * n_bonds`` ESP stencil points (or
    no ESP points when the molecule has 0 bonds) and returns the central-difference
    electric field ``E = -grad V`` projected along every RDKit bond. Supports both
    ``interface='native'`` and ``interface='compiled'``.
    """
    from rdkit import Chem

    step = _validate_step(stepAngstrom)
    if interface not in ('native', 'compiled'):
        raise ValueError("interface must be 'native' or 'compiled'")
    snapshot = Chem.Mol(mol)
    base = options or XTBOptions()
    _xyz_from_mol(snapshot, base.confId)
    charge, unpaired = _effective_charge_uhf(snapshot, base.charge, base.numUnpairedElectrons)
    conf_id = int(snapshot.GetConformer(base.confId).GetId())
    _, stencil = _bond_stencil(snapshot, conf_id, step)

    selected = replace(
        base,
        confId=conf_id,
        charge=charge,
        numUnpairedElectrons=unpaired,
        espPointsAngstrom=tuple(tuple(map(float, row)) for row in stencil),
    )
    if interface == 'native':
        result = XTBCalculator(**selected.ToDict()).Calculate(snapshot)
    else:
        from . import rdXTBTools

        ok, raw = rdXTBTools.RunMol(
            snapshot,
            confId=selected.confId,
            method=selected.method,
            charge=selected.charge,
            numUnpairedElectrons=selected.numUnpairedElectrons,
            keepOrbitalCoefficients=selected.keepOverlap,
            xtbLibrary=selected.xtbLibrary or '',
            maxIterations=selected.maxIterations,
            accuracy=selected.accuracy,
            electronicTemperature=selected.electronicTemperature,
            keepOverlap=selected.keepOverlap,
            espPointsAngstrom=stencil if len(stencil) else None,
            solvent=selected.solvent,
            solvationModel=selected.solvationModel,
            solvationGrid=selected.solvationGrid,
        )
        if not ok:
            error_type = XTBRunError if raw.failureKind == 'calculation' else XTBEnvironmentError
            raise error_type(raw.errorMessage)
        result = raw

    report = CalcBondElectricFieldFromESP(snapshot, result, stepAngstrom=step)
    report['analysis_settings']['interface'] = interface
    return report
