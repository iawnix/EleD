"""Sampled molecular ESP and constrained ESP / single-stage RESP charge fits.

Coordinates are Angstrom, potential is Hartree/e and fitted charges are in e.
No electronic-structure calculation is performed by these analysis functions.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from hashlib import sha256
import json
from numbers import Integral, Real
from typing import Any, Mapping

import numpy as np

from .errors import DescriptorUnavailableError
from .esp_input import _normalize_esp_points

ANGSTROM_TO_BOHR = 1.8897261254578281


def _real(value, name, positive=False):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not np.isfinite(value):
        raise ValueError(f'{name} must be finite and real')
    if positive and value <= 0:
        raise ValueError(f'{name} must be positive')
    return float(value)


def _vector(value, name, length):
    raw = np.asarray(value)
    if raw.dtype.kind not in 'iuf' or raw.shape != (length,):
        raise ValueError(f'{name} must be a real vector of length {length}')
    out = np.asarray(raw, dtype=float)
    if not np.isfinite(out).all():
        raise ValueError(f'{name} must be finite')
    return out


def CalcESP(result):
    """Return sampled potential in Hartree/e, in ``result.espPointsAngstrom`` order."""
    from .electronic_descriptors import _as_result
    result = _as_result(result)
    _, values = _esp_samples(result)
    return tuple(map(float, values))


def _esp_samples(result):
    points = _normalize_esp_points(result.espPointsAngstrom)
    stored = result.espValuesAU
    if stored is None:
        stored = ()
    values = _vector(stored, 'ESP values', len(points))
    if not points:
        raise DescriptorUnavailableError('ESP was not sampled; pass espPointsAngstrom when calculating')
    return np.asarray(points), values


def GenerateESPGrid(mol, *, confId=-1, scales=(1.4, 1.6, 1.8, 2.0), pointsPerAtom=100):
    """Generate deterministic exposed Fibonacci vdW shell points in Angstrom."""
    from rdkit import Chem
    from .input import _xyz_from_mol
    _xyz_from_mol(mol, confId)
    if (isinstance(pointsPerAtom, bool) or not isinstance(pointsPerAtom, Integral)
            or not 6 <= pointsPerAtom <= 100000):
        raise ValueError('pointsPerAtom must be an integer between 6 and 100000')
    scales = tuple(_real(v, 'scale', True) for v in scales)
    if not scales or any(v <= 1 for v in scales) or len(set(scales)) != len(scales):
        raise ValueError('scales must be distinct values greater than 1')
    coords = np.asarray(mol.GetConformer(confId).GetPositions())
    radii = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    if not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError('positive vdW radii are required for every atom')
    k = np.arange(pointsPerAtom, dtype=float)
    z = 1 - 2 * (k + 0.5) / pointsPerAtom
    angle = k * (np.pi * (3 - np.sqrt(5)))
    r = np.sqrt(1 - z * z)
    directions = np.column_stack((r * np.cos(angle), r * np.sin(angle), z))
    blocks = []
    for scale in scales:
        for atom, center in enumerate(coords):
            points = center + scale * radii[atom] * directions
            distances = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)
            exposed = np.all(distances >= scale * radii[None, :] - 1e-10, axis=1)
            blocks.extend(map(tuple, points[exposed]))
    if not blocks:
        raise ValueError('ESP grid has no exposed points')
    return _normalize_esp_points(blocks)


@dataclass(frozen=True)
class ESPChargeFit:
    """Owned fit result; atom order matches the source calculation."""

    chargeModel: str
    charges: tuple[float, ...]
    fittedPotentialAU: tuple[float, ...]
    residualsAU: tuple[float, ...]
    weightedRMSAU: float
    maxAbsResidualAU: float
    constraintError: float
    dataRank: int
    freeParameters: int
    conditionNumber: float | None
    iterations: int
    objective: float
    settings: Mapping[str, Any]
    provenance: Mapping[str, Any]

    def GetCharges(self) -> np.ndarray:
        """Return an owned NumPy copy of fitted atomic charges in ``e``."""
        return np.array(self.charges, dtype=float, copy=True)

    def ToDict(self) -> dict[str, Any]:
        """Return a JSON-compatible record with explicit unit-suffixed keys."""
        payload = {
            'charge_model': self.chargeModel,
            'charges_e': list(self.charges),
            'fitted_potential_au': list(self.fittedPotentialAU),
            'residuals_au': list(self.residualsAU),
            'weighted_rms_au': self.weightedRMSAU,
            'max_abs_residual_au': self.maxAbsResidualAU,
            'constraint_error_e': self.constraintError,
            'data_rank': self.dataRank,
            'free_parameters': self.freeParameters,
            'condition_number': self.conditionNumber,
            'iterations': self.iterations,
            'objective': self.objective,
            'settings': dict(self.settings),
            'provenance': dict(self.provenance),
        }
        return json.loads(json.dumps(payload, allow_nan=False))


def _constraints(n, charge, groups, fixed):
    rows, rhs = [np.ones(n)], [charge]

    def index(value):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, Integral) or not 0 <= value < n:
            raise ValueError('constraint atom indices must be integers in the original atom order')
        return int(value)

    normalized = []
    for group in groups:
        group = tuple(index(i) for i in group)
        if len(group) < 2 or len(set(group)) != len(group):
            raise ValueError('equivalenceGroups need at least two distinct atom indices')
        normalized.append(group)
        for i in group[1:]:
            row = np.zeros(n)
            row[group[0]] = 1
            row[i] = -1
            rows.append(row)
            rhs.append(0.0)
    if fixed is None:
        fixed = {}
    if not isinstance(fixed, Mapping):
        raise ValueError('fixedCharges must map atom indices to charges in e')
    normalized_fixed = {}
    for key, value in fixed.items():
        i, q = index(key), _real(value, 'fixed charge')
        row = np.zeros(n)
        row[i] = 1
        rows.append(row)
        rhs.append(q)
        normalized_fixed[i] = q
    return np.array(rows), np.array(rhs), normalized, normalized_fixed


def _fit(
    result,
    *,
    restraintStrength,
    restraintWidth,
    restrainHydrogens,
    weights,
    equivalenceGroups,
    fixedCharges,
    maxIterations,
    tolerance,
    rcond,
):
    from .electronic_descriptors import _as_result
    result = _as_result(result)
    points, values = _esp_samples(result)
    coords = np.asarray(_normalize_esp_points(result.metadata.get('coordinates_angstrom')))
    numbers = result.metadata.get('atomic_numbers', ())
    n = len(coords)
    if not n or numbers is None or np.ndim(numbers) != 1 or len(numbers) != n:
        raise ValueError('charge fitting requires source geometry and atomic numbers')
    if result.metadata.get('atom_count', n) != n:
        raise ValueError('source atom count does not match fitting coordinates')
    if any(isinstance(z, (bool, np.bool_)) or not isinstance(z, Integral) or not 1 <= z <= 118 for z in numbers):
        raise ValueError('atomic numbers must be integers in [1, 118]')
    charge = _real(result.metadata.get('charge'), 'source total charge')
    a = _real(restraintStrength, 'restraintStrength')
    b = _real(restraintWidth, 'restraintWidth', True)
    tol = _real(tolerance, 'tolerance', True)
    cutoff = _real(rcond, 'rcond', True)
    if a < 0 or cutoff >= 1:
        raise ValueError('restraintStrength must be nonnegative and rcond must be below 1')
    if not isinstance(restrainHydrogens, bool):
        raise ValueError('restrainHydrogens must be bool')
    if isinstance(maxIterations, bool) or not isinstance(maxIterations, Integral) or maxIterations < 1:
        raise ValueError('maxIterations must be a positive integer')
    w = np.ones(len(points)) if weights is None else _vector(weights, 'weights', len(points))
    if np.any(w < 0) or not np.any(w > 0):
        raise ValueError('weights must be nonnegative with a positive total')
    distances = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2) * ANGSTROM_TO_BOHR
    if not np.isfinite(distances).all() or np.any(distances < 1e-10):
        raise ValueError('sampling points must have finite nonzero distances from all nuclei')
    design = 1 / distances
    c, d, groups, fixed = _constraints(n, charge, equivalenceGroups, fixedCharges)
    u, s, vh = np.linalg.svd(c, full_matrices=True)
    rank = int(np.count_nonzero(s > cutoff * s[0]))
    q0 = vh[:rank].T @ ((u[:, :rank].T @ d) / s[:rank])
    if np.max(np.abs(c @ q0 - d)) > 1e-9:
        raise ValueError('inconsistent total-charge, equivalence or fixed-charge constraints')
    null = vh[rank:].T
    free = null.shape[1]
    sqrt_w = np.sqrt(w)
    matrix = sqrt_w[:, None] * (design @ null)
    target = sqrt_w * (values - design @ q0)
    solution, _, data_rank, singular = np.linalg.lstsq(matrix, target, rcond=cutoff)
    data_rank = int(data_rank)
    if a == 0 and data_rank < free:
        raise ValueError('ESP fit is rank deficient; add independent sampling points or constraints')
    condition = (
        float(singular[0] / singular[-1])
        if free and data_rank == free
        else 1.0
        if not free
        else None
    )
    q = q0 + null @ solution
    penalty = np.full(n, a)
    if not restrainHydrogens:
        penalty[np.asarray(numbers) == 1] = 0
    iterations = 1
    if a > 0 and free:
        for iterations in range(1, maxIterations + 1):
            diagonal = np.sqrt(penalty / np.sqrt(q * q + b * b))
            augmented = np.vstack((matrix, diagonal[:, None] * null))
            rhs = np.concatenate((target, -diagonal * q0))
            solution, _, fit_rank, _ = np.linalg.lstsq(augmented, rhs, rcond=cutoff)
            if fit_rank < free:
                raise ValueError('RESP restraints and constraints do not determine all free charges')
            updated = q0 + null @ solution
            change = float(np.max(np.abs(updated - q)))
            q = updated
            if change <= tol:
                break
        else:
            raise RuntimeError(f'RESP fit did not converge in {maxIterations} iterations')
    predicted = design @ q
    residual = predicted - values
    objective = 0.5 * float(np.dot(w, residual ** 2)) + float(np.dot(penalty, np.sqrt(q * q + b * b) - b))
    constraint_error = float(np.max(np.abs(c @ q - d)))
    if not np.isfinite(q).all() or not np.isfinite(objective) or constraint_error > 1e-8:
        raise ValueError('charge fit produced invalid values or violated constraints')
    settings = dict(
        protocol='eled-resp-single-stage-v1' if a else 'eled-esp-v1',
        restraint_strength=a,
        restraint_width_e=b,
        restrain_hydrogens=restrainHydrogens,
        equivalence_groups=groups,
        fixed_charges_e=fixed,
        total_charge_e=charge,
        weights=w.tolist(),
        max_iterations=int(maxIterations),
        tolerance_e=tol,
        rcond=cutoff,
        objective='0.5*sum(w*(Aq-V)^2)+sum(a_i*(sqrt(q_i^2+b^2)-b))',
        point_count=len(points),
        coordinate_unit='angstrom',
        potential_unit='hartree/e',
        charge_unit='e',
    )
    inputs = dict(
        points_angstrom=points.tolist(),
        potential_au=values.tolist(),
        coordinates_angstrom=coords.tolist(),
        atomic_numbers=list(map(int, numbers)),
        settings=settings,
    )
    digest = sha256(json.dumps(inputs, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    provenance = dict(
        input_sha256=digest,
        method=result.method,
        backend=result.backend,
        engine=result.metadata.get('engine'),
        parameter_file=result.metadata.get('parameter_file'),
        calculation_settings=result.metadata.get('calculation_settings'),
    )
    return ESPChargeFit(
        'resp' if a else 'esp',
        tuple(map(float, q)),
        tuple(map(float, predicted)),
        tuple(map(float, residual)),
        float(np.sqrt(np.dot(w, residual ** 2) / w.sum())),
        float(np.max(np.abs(residual))),
        constraint_error,
        data_rank,
        free,
        condition,
        iterations,
        objective,
        settings,
        deepcopy(provenance),
    )


def FitESPCharges(result, *, weights=None, equivalenceGroups=(), fixedCharges=None, rcond=1e-12):
    """Fit point charges to sampled ESP with exact total/equality/fixed constraints."""
    return _fit(
        result,
        restraintStrength=0.0,
        restraintWidth=0.1,
        restrainHydrogens=False,
        weights=weights,
        equivalenceGroups=equivalenceGroups,
        fixedCharges=fixedCharges,
        maxIterations=1,
        tolerance=1e-8,
        rcond=rcond,
    )


def FitRESPCharges(
    result,
    *,
    restraintStrength=0.0005,
    restraintWidth=0.1,
    restrainHydrogens=False,
    weights=None,
    equivalenceGroups=(),
    fixedCharges=None,
    maxIterations=1000,
    tolerance=1e-8,
    rcond=1e-12,
):
    """Single-stage hyperbolic RESP fit, not the full two-stage Amber protocol."""
    _real(restraintStrength, 'restraintStrength', True)
    return _fit(
        result,
        restraintStrength=restraintStrength,
        restraintWidth=restraintWidth,
        restrainHydrogens=restrainHydrogens,
        weights=weights,
        equivalenceGroups=equivalenceGroups,
        fixedCharges=fixedCharges,
        maxIterations=maxIterations,
        tolerance=tolerance,
        rcond=rcond,
    )
