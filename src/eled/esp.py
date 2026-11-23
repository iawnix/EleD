"""Sampled molecular ESP and constrained ESP / single-stage RESP charge fits.

Coordinates are Angstrom, potential is Hartree/e and fitted charges are in e.
No electronic-structure calculation is performed by these analysis functions.
"""
from copy import deepcopy
from dataclasses import asdict, dataclass
from hashlib import sha256
import json
from numbers import Integral, Real
from typing import Mapping

import numpy as np

from .errors import DescriptorUnavailableError
from .esp_input import ESPPoints

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
    """Return sampled potential in Hartree/e, in result.esp_points_angstrom order.

    Request esp_points_angstrom when calculating; an unsampled result raises
    DescriptorUnavailableError. Compiled XTBResults are accepted as well.
    """
    from .electronic_descriptors import _as_result
    result = _as_result(result)
    _, values = _esp_samples(result)
    return tuple(map(float, values))


def _esp_samples(result):
    points = ESPPoints(result.esp_points_angstrom)
    stored = result.esp_values_au
    if stored is None:
        stored = ()
    values = _vector(stored, 'ESP values', len(points))
    if not points:
        raise DescriptorUnavailableError('ESP was not sampled; pass esp_points_angstrom when calculating')
    return np.asarray(points), values


def GenerateESPGrid(mol, *, conf_id=-1, scales=(1.4, 1.6, 1.8, 2.0), points_per_atom=100):
    """Generate deterministic exposed Fibonacci vdW shell points in Angstrom.

    RDKit vdW radii, explicit hydrogen atoms, fixed lab-frame directions. Points
    inside another atom's scaled vdW sphere on the same shell are removed.
    This sampling rule is neither CHELPG nor Merz-Kollman; record its settings.
    """
    from rdkit import Chem
    from .input import _xyz_from_mol
    _xyz_from_mol(mol, conf_id)
    if (isinstance(points_per_atom, bool) or not isinstance(points_per_atom, Integral)
            or not 6 <= points_per_atom <= 100000):
        raise ValueError('points_per_atom must be an integer between 6 and 100000')
    scales = tuple(_real(v, 'scale', True) for v in scales)
    if not scales or any(v <= 1 for v in scales) or len(set(scales)) != len(scales):
        raise ValueError('scales must be distinct values greater than 1')
    coords = np.asarray(mol.GetConformer(conf_id).GetPositions())
    radii = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    if not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError('positive vdW radii are required for every atom')
    k = np.arange(points_per_atom, dtype=float)
    z = 1 - 2*(k+.5)/points_per_atom
    angle = k*(np.pi*(3-np.sqrt(5)))
    r = np.sqrt(1-z*z)
    directions = np.column_stack((r*np.cos(angle), r*np.sin(angle), z))
    blocks = []
    for scale in scales:
        for atom, center in enumerate(coords):
            points = center + scale*radii[atom]*directions
            # Per-atom blocks bound temporary memory to points_per_atom * natoms.
            distances = np.linalg.norm(points[:, None, :]-coords[None, :, :], axis=2)
            exposed = np.all(distances >= scale*radii[None, :]-1e-10, axis=1)
            blocks.extend(map(tuple, points[exposed]))
    if not blocks:
        raise ValueError('ESP grid has no exposed points')
    return ESPPoints(blocks)


@dataclass(frozen=True)
class ESPChargeFit:
    """Owned fit result; atom order matches the source calculation.

    condition_number is None for rank-deficient data (RESP may still be unique
    because of restraints). No unconverged fit is returned as successful.
    """
    charge_model: str
    charges_e: tuple[float, ...]
    fitted_potential_au: tuple[float, ...]
    residuals_au: tuple[float, ...]  # fitted minus reference
    weighted_rms_au: float
    max_abs_residual_au: float
    constraint_error_e: float
    data_rank: int
    free_parameters: int
    condition_number: float | None
    iterations: int
    objective: float
    settings: Mapping
    provenance: Mapping

    def to_dict(self):
        """Return a JSON-compatible record including all fit settings and input hash."""
        return json.loads(json.dumps(asdict(self), allow_nan=False))


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
            raise ValueError('equivalence groups need at least two distinct atom indices')
        normalized.append(group)
        for i in group[1:]:
            row = np.zeros(n); row[group[0]] = 1; row[i] = -1
            rows.append(row); rhs.append(0.)
    if fixed is None:
        fixed = {}
    if not isinstance(fixed, Mapping):
        raise ValueError('fixed_charges must map atom indices to charges in e')
    normalized_fixed = {}
    for key, value in fixed.items():
        i, q = index(key), _real(value, 'fixed charge')
        row = np.zeros(n); row[i] = 1
        rows.append(row); rhs.append(q); normalized_fixed[i] = q
    return np.array(rows), np.array(rhs), normalized, normalized_fixed


def _fit(result, *, restraint_strength, restraint_width, restrain_hydrogens,
         weights, equivalence_groups, fixed_charges, max_iterations, tolerance, rcond):
    from .electronic_descriptors import _as_result
    result = _as_result(result)
    points, values = _esp_samples(result)
    coords = np.asarray(ESPPoints(result.metadata.get('coordinates_angstrom')))
    numbers = result.metadata.get('atomic_numbers', ())
    n = len(coords)
    if not n or numbers is None or np.ndim(numbers) != 1 or len(numbers) != n:
        raise ValueError('charge fitting requires source geometry and atomic numbers')
    if result.metadata.get('atom_count', n) != n:
        raise ValueError('source atom count does not match fitting coordinates')
    if any(isinstance(z, (bool, np.bool_)) or not isinstance(z, Integral) or not 1 <= z <= 118 for z in numbers):
        raise ValueError('atomic numbers must be integers in [1, 118]')
    charge = _real(result.metadata.get('charge'), 'source total charge')
    a = _real(restraint_strength, 'restraint_strength')
    b = _real(restraint_width, 'restraint_width', True)
    tol = _real(tolerance, 'tolerance', True)
    cutoff = _real(rcond, 'rcond', True)
    if a < 0 or cutoff >= 1:
        raise ValueError('restraint_strength must be nonnegative and rcond must be below 1')
    if not isinstance(restrain_hydrogens, bool):
        raise ValueError('restrain_hydrogens must be bool')
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral) or max_iterations < 1:
        raise ValueError('max_iterations must be a positive integer')
    w = np.ones(len(points)) if weights is None else _vector(weights, 'weights', len(points))
    if np.any(w < 0) or not np.any(w > 0):
        raise ValueError('weights must be nonnegative with a positive total')
    distances = np.linalg.norm(points[:, None, :]-coords[None, :, :], axis=2)*ANGSTROM_TO_BOHR
    if not np.isfinite(distances).all() or np.any(distances < 1e-10):
        raise ValueError('sampling points must have finite nonzero distances from all nuclei')
    design = 1/distances
    c, d, groups, fixed = _constraints(n, charge, equivalence_groups, fixed_charges)
    # Solve constraints first, then fit only in their null space. This avoids
    # squaring the condition number of the electrostatic design matrix.
    u, s, vh = np.linalg.svd(c, full_matrices=True)
    rank = int(np.count_nonzero(s > cutoff*s[0]))
    q0 = vh[:rank].T @ ((u[:, :rank].T @ d)/s[:rank])
    if np.max(np.abs(c@q0-d)) > 1e-9:
        raise ValueError('inconsistent total-charge, equivalence or fixed-charge constraints')
    null = vh[rank:].T
    free = null.shape[1]
    sqrt_w = np.sqrt(w)
    matrix = sqrt_w[:, None]*(design@null)
    target = sqrt_w*(values-design@q0)
    solution, _, data_rank, singular = np.linalg.lstsq(matrix, target, rcond=cutoff)
    data_rank = int(data_rank)
    if a == 0 and data_rank < free:
        raise ValueError('ESP fit is rank deficient; add independent sampling points or constraints')
    condition = (float(singular[0]/singular[-1]) if free and data_rank == free else
                 1.0 if not free else None)
    q = q0 + null@solution
    penalty = np.full(n, a)
    if not restrain_hydrogens:
        penalty[np.asarray(numbers) == 1] = 0
    iterations = 1
    if a > 0 and free:
        for iterations in range(1, max_iterations+1):
            diagonal = np.sqrt(penalty/np.sqrt(q*q+b*b))
            augmented = np.vstack((matrix, diagonal[:, None]*null))
            rhs = np.concatenate((target, -diagonal*q0))
            solution, _, fit_rank, _ = np.linalg.lstsq(augmented, rhs, rcond=cutoff)
            if fit_rank < free:
                raise ValueError('RESP restraints and constraints do not determine all free charges')
            updated = q0 + null@solution
            change = float(np.max(np.abs(updated-q)))
            q = updated
            if change <= tol:
                break
        else:
            raise RuntimeError(f'RESP fit did not converge in {max_iterations} iterations')
    predicted = design@q
    residual = predicted-values
    objective = .5*float(np.dot(w, residual**2)) + float(np.dot(penalty, np.sqrt(q*q+b*b)-b))
    constraint_error = float(np.max(np.abs(c@q-d)))
    if not np.isfinite(q).all() or not np.isfinite(objective) or constraint_error > 1e-8:
        raise ValueError('charge fit produced invalid values or violated constraints')
    settings = dict(protocol='eled-resp-single-stage-v1' if a else 'eled-esp-v1',
                    restraint_strength=a, restraint_width_e=b, restrain_hydrogens=restrain_hydrogens,
                    equivalence_groups=groups, fixed_charges_e=fixed, total_charge_e=charge,
                    weights=w.tolist(), max_iterations=int(max_iterations), tolerance_e=tol, rcond=cutoff,
                    objective='0.5*sum(w*(Aq-V)^2)+sum(a_i*(sqrt(q_i^2+b^2)-b))',
                    point_count=len(points), coordinate_unit='angstrom', potential_unit='hartree/e', charge_unit='e')
    inputs = dict(points_angstrom=points.tolist(), potential_au=values.tolist(),
                  coordinates_angstrom=coords.tolist(), atomic_numbers=list(map(int, numbers)), settings=settings)
    digest = sha256(json.dumps(inputs, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
    provenance = dict(input_sha256=digest, method=result.method, backend=result.backend,
                      engine=result.metadata.get('engine'), parameter_file=result.metadata.get('parameter_file'),
                      calculation_settings=result.metadata.get('CalculationSettings'))
    return ESPChargeFit('resp' if a else 'esp', tuple(map(float, q)), tuple(map(float, predicted)),
                        tuple(map(float, residual)), float(np.sqrt(np.dot(w, residual**2)/w.sum())),
                        float(np.max(np.abs(residual))), constraint_error, data_rank, free, condition,
                        iterations, objective, settings, deepcopy(provenance))


def FitESPCharges(result, *, weights=None, equivalence_groups=(), fixed_charges=None, rcond=1e-12):
    """Fit point charges to sampled ESP with exact total/equality/fixed constraints.

    weights multiply squared residuals (not residuals); rank-deficient fits raise.
    The total charge and indexed atom coordinates come from result provenance.
    """
    return _fit(result, restraint_strength=0., restraint_width=.1, restrain_hydrogens=False,
                weights=weights, equivalence_groups=equivalence_groups, fixed_charges=fixed_charges,
                max_iterations=1, tolerance=1e-8, rcond=rcond)


def FitRESPCharges(result, *, restraint_strength=0.0005, restraint_width=0.1,
                   restrain_hydrogens=False, weights=None, equivalence_groups=(), fixed_charges=None,
                   max_iterations=1000, tolerance=1e-8, rcond=1e-12):
    """Single-stage hyperbolic RESP fit, not the full two-stage Amber protocol.

    Minimize 0.5*sum(w*(Aq-V)^2) + sum(a_i*(sqrt(q_i^2+b^2)-b)).
    a_i=restraint_strength except unrestrained H; b=restraint_width in e.
    Weights are not normalized: rescaling them changes relative restraint strength.
    User supplies chemical equivalences; no automatic atom regrouping is applied.
    Nonconvergence raises RuntimeError; insufficient independent information raises ValueError.
    """
    _real(restraint_strength, 'restraint_strength', True)
    return _fit(result, restraint_strength=restraint_strength, restraint_width=restraint_width,
                restrain_hydrogens=restrain_hydrogens, weights=weights,
                equivalence_groups=equivalence_groups, fixed_charges=fixed_charges,
                max_iterations=max_iterations, tolerance=tolerance, rcond=rcond)
