"""Real-space wavefunction descriptors: density, ALIE, LEA, and 0.001 a.u. isosurface.

Uses EleD's reviewed ``realspace_extension`` ABI v1 in ``libxtb`` to evaluate
exact Cartesian GTO valence electron density ``rho(r)`` (``e/Bohr^3``),
Average Local Ionization Energy ``ALIE(r)`` (``eV``), virtual orbital density
(``1/Bohr^3``), and Local Electron Affinity ``LEA(r)`` (``eV``) at arbitrary
3D points without additional SCF calculations.

- :func:`CalcRealSpaceProperties`: pure analysis of an existing result with
  sampled real-space points (never calls xTB). Points below ``densityCutoffAU``
  report ``None`` rather than zero.
- :func:`CalcRealSpaceSurfaceDescriptors`: area-weighted molecular and per-atom
  ALIE/LEA surface statistics on any :class:`ESPSurfaceGrid` (vdW Fibonacci
  shell or 0.001 a.u. electron-density isosurface) or explicit point set.
- :func:`CalcDensityIsosurface`: finds the exact ``isovalueAU`` (default
  ``0.001 e/Bohr^3``) electron-density isosurface in **one** native ``libxtb``
  single-point calculation via Voronoi-partitioned batched radial bisection,
  returning molecular volume, surface area, sphericity, surface grid, and
  optional isosurface Politzer-Murray ESP and ALIE/LEA surface descriptors.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from numbers import Integral, Real
from typing import Any

import numpy as np

from .backends.native import ANGSTROM_TO_BOHR, _run_xtb_native
from .electronic_descriptors import _as_result
from .errors import DescriptorUnavailableError
from .esp import _vector
from .esp_descriptors import (
    CalcESPSurfaceDescriptors,
    ESPSurfaceGrid,
    _assign_nearest_vdw_atoms,
    _require_mol,
    _validate_atom_indices,
    _validate_weights,
    _vdw_radii,
)
from .esp_input import _normalize_esp_points
from .input import _effective_charge_uhf, _xyz_from_mol
from .mapping import _validate_mapping
from .options import XTBOptions
from .xtb import _resolve_options

SURFACE_MODEL_DENSITY_ISOSURFACE = 'electron_density_isosurface'
REALSPACE_SURFACE_PROTOCOL = 'eled-realspace-surface-v1'
ISOSURFACE_PROTOCOL = 'eled-density-isosurface-v1'

__all__ = [
    'SURFACE_MODEL_DENSITY_ISOSURFACE',
    'CalcRealSpaceProperties',
    'CalcRealSpaceSurfaceDescriptors',
    'CalcDensityIsosurface',
]


def _validate_nonnegative_real(value: Any, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value) or value < 0.0:
        raise ValueError(f'{name} must be a finite nonnegative real number')
    return float(value)


def _validate_positive_real(value: Any, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value) or value <= 0.0:
        raise ValueError(f'{name} must be a finite positive real number')
    return float(value)


def _realspace_samples(result):
    points = _normalize_esp_points(result.realSpacePointsAngstrom)
    density = result.realSpaceDensityAU
    alie = result.realSpaceALIEEV
    vir_density = result.realSpaceVirtualDensityAU
    lea = result.realSpaceLEAEV
    if not points or not density:
        raise DescriptorUnavailableError(
            'real-space properties were not sampled; pass realSpacePointsAngstrom when calculating'
        )
    n = len(points)
    return (
        np.asarray(points, dtype=float),
        _vector(density, 'real-space density', n),
        _vector(alie, 'real-space ALIE', n),
        _vector(vir_density, 'real-space virtual density', n),
        _vector(lea, 'real-space LEA', n),
    )


def CalcRealSpaceProperties(
    result: Any,
    *,
    densityCutoffAU: float = 1e-18,
) -> dict[str, Any]:
    """Return sampled real-space density, ALIE, virtual density, and LEA without engine calls.

    Points where ``density <= densityCutoffAU`` have ``alie_ev = None``; points
    where ``vir_density <= densityCutoffAU`` have ``lea_ev = None`` (undefined
    ratios are never encoded as zero).
    """
    cutoff = _validate_nonnegative_real(densityCutoffAU, 'densityCutoffAU')
    res = _as_result(result)
    points, density, alie, vir_density, lea = _realspace_samples(res)

    alie_valid = density > cutoff
    lea_valid = vir_density > cutoff
    alie_out = [float(alie[i]) if alie_valid[i] else None for i in range(len(density))]
    lea_out = [float(lea[i]) if lea_valid[i] else None for i in range(len(vir_density))]

    n_alie = int(np.count_nonzero(alie_valid))
    n_lea = int(np.count_nonzero(lea_valid))
    missing: dict[str, str] = {}
    if n_alie < len(density):
        missing['alie_ev'] = 'density_below_cutoff_at_some_points'
    if n_lea < len(vir_density):
        missing['lea_ev'] = 'virtual_density_below_cutoff_at_some_points'

    payload = {
        'n_points': int(len(density)),
        'n_alie_defined': n_alie,
        'n_lea_defined': n_lea,
        'points_angstrom': points.tolist(),
        'density_au': [float(v) for v in density],
        'alie_ev': alie_out,
        'vir_density_au': [float(v) for v in vir_density],
        'lea_ev': lea_out,
        'density_cutoff_au': cutoff,
        'method': res.method,
        'backend': res.backend,
        'units': {
            'points_angstrom': 'angstrom',
            'density_au': 'e/Bohr^3',
            'alie_ev': 'eV',
            'vir_density_au': '1/Bohr^3',
            'lea_ev': 'eV',
            'density_cutoff_au': 'e/Bohr^3',
        },
        'analysis_settings': {
            'source': 'libxtb_realspace_cartesian_gto_v1',
            'density_cutoff_au': cutoff,
            'alie_definition': 'sum_{i in occ} f_i * (-eps_i) * |psi_i(r)|^2 / rho(r)',
            'lea_definition': 'sum_{a in vir} (2 - f_a) * (-eps_a) * |psi_a(r)|^2 / rho_vir(r)',
        },
        'invariance': (
            'invariant to rigid translation/rotation when sampling points co-transform with the molecule, '
            'MO phase, and degenerate orbital subspace rotations'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))


def _weighted_scalar_stats(values: np.ndarray, weights: np.ndarray) -> dict[str, float]:
    w_sum = float(np.sum(weights))
    mean = float(np.sum(weights * values) / w_sum)
    var = float(np.sum(weights * (values - mean) ** 2) / w_sum)
    mad = float(np.sum(weights * np.abs(values - mean)) / w_sum)
    return {
        'min': float(np.min(values)),
        'max': float(np.max(values)),
        'mean': mean,
        'sigma2': var,
        'pi': mad,
    }


def CalcRealSpaceSurfaceDescriptors(
    mol: Any,
    result: Any,
    *,
    atomIndices: Any = None,
    weights: Any = None,
    surfaceGrid: ESPSurfaceGrid | None = None,
    densityCutoffAU: float = 1e-18,
) -> dict[str, Any]:
    """Compute molecular and per-atom surface ALIE and LEA statistics without engine calls.

    Follows the same grid/ownership/area-weight contract as
    :func:`eled.CalcESPSurfaceDescriptors`. Undefined points (density or virtual
    density ``<= densityCutoffAU``) and buried atoms are reported as ``None``
    with explicit ``missing`` reasons.
    """
    _require_mol(mol)
    cutoff = _validate_nonnegative_real(densityCutoffAU, 'densityCutoffAU')
    res = _as_result(result)
    _validate_mapping(mol, res)
    points, density, alie, vir_density, lea = _realspace_samples(res)

    natoms = int(mol.GetNumAtoms())
    npoints = int(len(density))

    if surfaceGrid is not None:
        if not isinstance(surfaceGrid, ESPSurfaceGrid):
            raise TypeError('surfaceGrid must be an ESPSurfaceGrid instance')
        if atomIndices is not None or weights is not None:
            raise ValueError('pass either surfaceGrid or explicit atomIndices/weights, not both')
        grid_points = surfaceGrid.GetPoints()
        if grid_points.shape != points.shape or not np.allclose(grid_points, points, atol=1e-10, rtol=0.0):
            raise ValueError('surfaceGrid coordinates do not match result real-space sampling points')
        assigned = _validate_atom_indices(surfaceGrid.atomIndices, npoints, natoms)
        w = _vector(surfaceGrid.areaWeightsAngstrom2, 'weights', npoints)
        atom_assignment = 'surface_grid'
        weight_mode = 'surface_grid_area_a2'
        surface_model = surfaceGrid.surfaceModel
    else:
        if atomIndices is not None:
            assigned = _validate_atom_indices(atomIndices, npoints, natoms)
            atom_assignment = 'explicit'
        else:
            conf_id = int(res.metadata['conformer_id'])
            coords = np.asarray(mol.GetConformer(conf_id).GetPositions(), dtype=float)
            radii = _vdw_radii(mol)
            assigned = _assign_nearest_vdw_atoms(points, coords, radii)
            atom_assignment = 'nearest_vdw_normalized'
        if weights is None:
            w = np.ones(npoints, dtype=float)
            weight_mode = 'uniform'
        else:
            w = _validate_weights(weights, npoints)
            weight_mode = 'explicit'
        surface_model = 'sampled_realspace_points'

    if np.any(w < 0) or not np.any(w > 0):
        raise ValueError('weights must be nonnegative with a positive total')

    active = w > 0.0
    w_active = w[active]
    d_active = density[active]
    w_total = float(np.sum(w_active))
    d_mean = float(np.sum(w_active * d_active) / w_total)

    alie_mask = active & (density > cutoff)
    lea_mask = active & (vir_density > cutoff)
    n_alie = int(np.count_nonzero(alie_mask))
    n_lea = int(np.count_nonzero(lea_mask))

    mol_missing: dict[str, str] = {}
    if n_alie:
        st = _weighted_scalar_stats(alie[alie_mask], w[alie_mask])
        alie_min, alie_max, alie_mean, alie_s2, alie_pi = (
            st['min'], st['max'], st['mean'], st['sigma2'], st['pi'],
        )
    else:
        alie_min = alie_max = alie_mean = alie_s2 = alie_pi = None
        for key in ('alie_min_ev', 'alie_max_ev', 'alie_mean_ev', 'alie_sigma2_ev2', 'alie_pi_ev'):
            mol_missing[key] = 'density_below_cutoff'

    if n_lea:
        st = _weighted_scalar_stats(lea[lea_mask], w[lea_mask])
        lea_min, lea_max, lea_mean, lea_s2, lea_pi = (
            st['min'], st['max'], st['mean'], st['sigma2'], st['pi'],
        )
    else:
        lea_min = lea_max = lea_mean = lea_s2 = lea_pi = None
        for key in ('lea_min_ev', 'lea_max_ev', 'lea_mean_ev', 'lea_sigma2_ev2', 'lea_pi_ev'):
            mol_missing[key] = 'virtual_density_below_cutoff'

    atoms_rows: list[dict[str, Any]] = []
    atom_missing: dict[str, dict[str, str]] = {}
    for atom_idx in range(natoms):
        assigned_mask = assigned == atom_idx
        active_atom = assigned_mask & active
        n_atom_active = int(np.count_nonzero(active_atom))
        row_missing: dict[str, str] = {}
        if n_atom_active == 0:
            reason = (
                'zero_surface_weight_assigned'
                if np.any(assigned_mask)
                else 'no_surface_points_assigned'
            )
            row = {
                'atom_index': atom_idx,
                'n_points': 0,
                'n_alie_points': 0,
                'n_lea_points': 0,
                'alie_min_ev': None,
                'alie_max_ev': None,
                'alie_mean_ev': None,
                'alie_sigma2_ev2': None,
                'lea_min_ev': None,
                'lea_max_ev': None,
                'lea_mean_ev': None,
                'lea_sigma2_ev2': None,
            }
            for key in (
                'alie_min_ev', 'alie_max_ev', 'alie_mean_ev', 'alie_sigma2_ev2',
                'lea_min_ev', 'lea_max_ev', 'lea_mean_ev', 'lea_sigma2_ev2',
            ):
                row_missing[key] = reason
        else:
            a_alie_mask = active_atom & (density > cutoff)
            a_lea_mask = active_atom & (vir_density > cutoff)
            n_a_alie = int(np.count_nonzero(a_alie_mask))
            n_a_lea = int(np.count_nonzero(a_lea_mask))
            if n_a_alie:
                sa = _weighted_scalar_stats(alie[a_alie_mask], w[a_alie_mask])
                aa_min, aa_max, aa_mean, aa_s2 = sa['min'], sa['max'], sa['mean'], sa['sigma2']
            else:
                aa_min = aa_max = aa_mean = aa_s2 = None
                for key in ('alie_min_ev', 'alie_max_ev', 'alie_mean_ev', 'alie_sigma2_ev2'):
                    row_missing[key] = 'density_below_cutoff'
            if n_a_lea:
                sl = _weighted_scalar_stats(lea[a_lea_mask], w[a_lea_mask])
                al_min, al_max, al_mean, al_s2 = sl['min'], sl['max'], sl['mean'], sl['sigma2']
            else:
                al_min = al_max = al_mean = al_s2 = None
                for key in ('lea_min_ev', 'lea_max_ev', 'lea_mean_ev', 'lea_sigma2_ev2'):
                    row_missing[key] = 'virtual_density_below_cutoff'
            row = {
                'atom_index': atom_idx,
                'n_points': n_atom_active,
                'n_alie_points': n_a_alie,
                'n_lea_points': n_a_lea,
                'alie_min_ev': aa_min,
                'alie_max_ev': aa_max,
                'alie_mean_ev': aa_mean,
                'alie_sigma2_ev2': aa_s2,
                'lea_min_ev': al_min,
                'lea_max_ev': al_max,
                'lea_mean_ev': al_mean,
                'lea_sigma2_ev2': al_s2,
            }
        row['missing'] = row_missing
        if row_missing:
            atom_missing[str(atom_idx)] = dict(row_missing)
        atoms_rows.append(row)

    missing: dict[str, Any] = dict(mol_missing)
    if atom_missing:
        missing['atoms'] = atom_missing

    payload = {
        'n_points': int(np.count_nonzero(active)),
        'n_points_total': npoints,
        'n_alie_points': n_alie,
        'n_lea_points': n_lea,
        'density_min_au': float(np.min(d_active)),
        'density_max_au': float(np.max(d_active)),
        'density_mean_au': d_mean,
        'alie_min_ev': alie_min,
        'alie_max_ev': alie_max,
        'alie_mean_ev': alie_mean,
        'alie_sigma2_ev2': alie_s2,
        'alie_pi_ev': alie_pi,
        'lea_min_ev': lea_min,
        'lea_max_ev': lea_max,
        'lea_mean_ev': lea_mean,
        'lea_sigma2_ev2': lea_s2,
        'lea_pi_ev': lea_pi,
        'atoms': atoms_rows,
        'units': {
            'n_points': 'count',
            'n_points_total': 'count',
            'n_alie_points': 'count',
            'n_lea_points': 'count',
            'density_min_au': 'e/Bohr^3',
            'density_max_au': 'e/Bohr^3',
            'density_mean_au': 'e/Bohr^3',
            'alie_min_ev': 'eV',
            'alie_max_ev': 'eV',
            'alie_mean_ev': 'eV',
            'alie_sigma2_ev2': 'eV^2',
            'alie_pi_ev': 'eV',
            'lea_min_ev': 'eV',
            'lea_max_ev': 'eV',
            'lea_mean_ev': 'eV',
            'lea_sigma2_ev2': 'eV^2',
            'lea_pi_ev': 'eV',
        },
        'analysis_settings': {
            'protocol': REALSPACE_SURFACE_PROTOCOL,
            'surface_model': surface_model,
            'atom_assignment': atom_assignment,
            'weight_mode': weight_mode,
            'total_weight': w_total,
            'density_cutoff_au': cutoff,
            'variance_convention': 'weighted population variance (ddof=0)',
            'method': res.method,
            'backend': res.backend,
            'conformer_id': int(res.metadata['conformer_id']),
        },
        'invariance': (
            'invariant to rigid translation/rotation when surface points co-transform with the molecule, '
            'point permutation, MO phase, degenerate orbital subspace rotation, and atom renumbering'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))


def _fibonacci_directions(n_per_atom: int) -> np.ndarray:
    k = np.arange(n_per_atom, dtype=float)
    z = 1.0 - 2.0 * (k + 0.5) / n_per_atom
    angle = k * (np.pi * (3.0 - np.sqrt(5.0)))
    r = np.sqrt(1.0 - z * z)
    return np.column_stack((r * np.cos(angle), r * np.sin(angle), z))


def _voronoi_ray_limits(coords: np.ndarray, directions: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return ray origins, unit directions, and Euclidean Voronoi cell exit distances (Angstrom)."""
    natoms = len(coords)
    n_dir = len(directions)
    n_rays = natoms * n_dir
    owners = np.repeat(np.arange(natoms, dtype=int), n_dir)
    origins = coords[owners]
    ray_dirs = np.tile(directions, (natoms, 1))
    t_cell = np.full(n_rays, np.inf, dtype=float)

    if natoms > 1:
        for a in range(natoms):
            sl = slice(a * n_dir, (a + 1) * n_dir)
            diff = coords - coords[a]  # [natoms, 3]
            dist2 = np.sum(diff * diff, axis=1)  # [natoms]
            # Exclude self
            other = np.arange(natoms) != a
            diff_o = diff[other]  # [natoms-1, 3]
            dist2_o = dist2[other]  # [natoms-1]
            # dot[k, j] = directions[k] . (R_j - R_a)
            dots = directions @ diff_o.T  # [n_dir, natoms-1]
            with np.errstate(divide='ignore', invalid='ignore'):
                limits = np.where(dots > 1e-12, 0.5 * dist2_o[None, :] / dots, np.inf)
            t_cell[sl] = np.min(limits, axis=1)

    return owners, origins, ray_dirs, t_cell


def CalcDensityIsosurface(
    mol: Any,
    *,
    isovalueAU: float = 0.001,
    pointsPerAtom: int = 120,
    includeESP: bool = True,
    includeALIE: bool = True,
    occThreshold: float = 1e-7,
    leaWindowEV: float = 0.0,
    densityCutoffAU: float = 1e-18,
    includeKcalMol: bool = True,
    bisectionSteps: int = 18,
    options: XTBOptions | None = None,
) -> dict[str, Any]:
    """Compute exact electron-density isosurface (default ``0.001 a.u.``) in one xTB SCF.

    Uses a single native ``libxtb`` single-point calculation and evaluates the
    exact Cartesian GTO valence density ``rho(r)`` in memory to locate the
    ``isovalueAU`` isosurface on a Voronoi-partitioned Fibonacci ray mesh.

    Returns molecular volume (``volume_a3``, ``volume_bohr3``), surface area
    (``area_a2``, ``area_bohr2``), sphericity (``sphericity``), the exact
    isosurface grid (``grid``), and optional Politzer-Murray surface ESP
    (``esp_surface``) and ALIE/LEA surface descriptors (``realspace_surface``)
    evaluated on the exact isosurface points.
    """
    from rdkit import Chem

    _require_mol(mol)
    iso = _validate_positive_real(isovalueAU, 'isovalueAU')
    occ_thr = _validate_nonnegative_real(occThreshold, 'occThreshold')
    lea_win = _validate_nonnegative_real(leaWindowEV, 'leaWindowEV')
    cutoff = _validate_nonnegative_real(densityCutoffAU, 'densityCutoffAU')
    if (
        isinstance(pointsPerAtom, (bool, np.bool_))
        or not isinstance(pointsPerAtom, Integral)
        or not 12 <= pointsPerAtom <= 20000
    ):
        raise ValueError('pointsPerAtom must be an integer between 12 and 20000')
    if (
        isinstance(bisectionSteps, (bool, np.bool_))
        or not isinstance(bisectionSteps, Integral)
        or not 8 <= bisectionSteps <= 40
    ):
        raise ValueError('bisectionSteps must be an integer between 8 and 40')
    for flag_name, flag_val in (
        ('includeESP', includeESP),
        ('includeALIE', includeALIE),
        ('includeKcalMol', includeKcalMol),
    ):
        if not isinstance(flag_val, bool):
            raise ValueError(f'{flag_name} must be bool')

    snapshot = Chem.Mol(mol)
    base = options or XTBOptions()
    _xyz_from_mol(snapshot, base.confId)
    charge, unpaired = _effective_charge_uhf(snapshot, base.charge, base.numUnpairedElectrons)
    conf = snapshot.GetConformer(base.confId)
    conf_id = int(conf.GetId())
    coords = np.asarray(conf.GetPositions(), dtype=float)
    radii = _vdw_radii(snapshot)

    n_per_atom = int(pointsPerAtom)
    directions = _fibonacci_directions(n_per_atom)
    owners, origins, ray_dirs, t_cell = _voronoi_ray_limits(coords, directions)
    n_rays = len(owners)
    d_omega = 4.0 * math.pi / float(n_per_atom)

    selected = _resolve_options(
        replace(
            base,
            confId=conf_id,
            charge=charge,
            numUnpairedElectrons=unpaired,
            espPointsAngstrom=(),
            realSpacePointsAngstrom=(),
            realSpaceOccupationThreshold=occ_thr,
            realSpaceLEAWindowEV=lea_win,
        )
    )

    iso_state: dict[str, Any] = {}

    def _eval_rho_batch(library, environment, results, pts_angstrom: np.ndarray) -> np.ndarray:
        if len(pts_angstrom) == 0:
            return np.empty(0, dtype=float)
        pts_tuples = tuple(tuple(map(float, row)) for row in pts_angstrom)
        data = library.realspace_data(
            environment, results, pts_tuples, occ_threshold=occ_thr, lea_window_ev=lea_win
        )
        return np.asarray(data['realSpaceDensityAU'], dtype=float)

    def _hook(library, environment, results, natoms_cb, nao_cb):
        library.bind_realspace_api()

        # Initial lower and upper radial brackets (Angstrom).
        t_lo = np.minimum(0.05, 0.5 * t_cell)
        t_hi_init = np.maximum(3.5, 2.2 * radii[owners])
        t_hi = np.minimum(t_cell, t_hi_init)

        # Evaluate density at t_lo and t_hi for all rays.
        rho_lo = _eval_rho_batch(library, environment, results, origins + t_lo[:, None] * ray_dirs)
        rho_hi = _eval_rho_batch(library, environment, results, origins + t_hi[:, None] * ray_dirs)

        # Expand t_hi outward for exterior rays where density is still > isovalueAU.
        for _ in range(7):
            need_expand = (rho_hi > iso) & (t_hi < t_cell - 1e-12)
            if not np.any(need_expand):
                break
            idx = np.flatnonzero(need_expand)
            t_hi[idx] = np.minimum(t_cell[idx], t_hi[idx] * 1.6)
            pts_exp = origins[idx] + t_hi[idx, None] * ray_dirs[idx]
            rho_hi[idx] = _eval_rho_batch(library, environment, results, pts_exp)

        # Rays that hit the Voronoi wall with rho > iso are buried inside the molecule.
        buried = (rho_hi > iso) & np.isfinite(t_cell) & (t_hi >= t_cell - 1e-12)
        # Active exterior rays bracket the isosurface root: rho(t_lo) >= iso >= rho(t_hi).
        active = (~buried) & (rho_lo >= iso) & (rho_hi <= iso)

        t_star = np.zeros(n_rays, dtype=float)
        t_star[buried] = t_cell[buried]
        # Rare case: if rho(t_lo) < iso (very high isovalue), contract t_lo toward nucleus.
        contract = (~buried) & (rho_lo < iso)
        if np.any(contract):
            c_idx = np.flatnonzero(contract)
            t_lo[c_idx] = 1e-4
            pts_c = origins[c_idx] + t_lo[c_idx, None] * ray_dirs[c_idx]
            rho_lo[c_idx] = _eval_rho_batch(library, environment, results, pts_c)
            newly_active = rho_lo[c_idx] >= iso
            active[c_idx[newly_active]] = True

        active_idx = np.flatnonzero(active)
        if not len(active_idx):
            raise ValueError(
                f'could not locate electron-density isosurface at isovalueAU={iso:g}; '
                'try a smaller isovalueAU'
            )

        lo_a = t_lo[active_idx].copy()
        hi_a = t_hi[active_idx].copy()
        orig_a = origins[active_idx]
        dir_a = ray_dirs[active_idx]

        for _ in range(int(bisectionSteps)):
            mid_a = 0.5 * (lo_a + hi_a)
            pts_mid = orig_a + mid_a[:, None] * dir_a
            rho_mid = _eval_rho_batch(library, environment, results, pts_mid)
            above = rho_mid > iso
            lo_a = np.where(above, mid_a, lo_a)
            hi_a = np.where(above, hi_a, mid_a)

        r_surf = 0.5 * (lo_a + hi_a)
        t_star[active_idx] = r_surf
        surf_pts = orig_a + r_surf[:, None] * dir_a
        surf_owners = owners[active_idx]

        # Estimate outward surface normal n = -grad(rho) / ||grad(rho)|| via 6-point finite differences
        # for accurate surface area elements dA = r^2 / (n . d) dOmega.
        h_grad = 1e-3  # Angstrom
        n_surf = len(surf_pts)
        grad_pts = np.empty((6 * n_surf, 3), dtype=float)
        axes = np.eye(3, dtype=float) * h_grad
        for axis in range(3):
            grad_pts[2 * axis :: 6] = surf_pts + axes[axis]
            grad_pts[2 * axis + 1 :: 6] = surf_pts - axes[axis]
        rho_grad_vals = _eval_rho_batch(library, environment, results, grad_pts)
        gx = (rho_grad_vals[0::6] - rho_grad_vals[1::6]) / (2.0 * h_grad)
        gy = (rho_grad_vals[2::6] - rho_grad_vals[3::6]) / (2.0 * h_grad)
        gz = (rho_grad_vals[4::6] - rho_grad_vals[5::6]) / (2.0 * h_grad)
        grad_vec = np.column_stack((gx, gy, gz))
        grad_norm = np.linalg.norm(grad_vec, axis=1)
        outward_normal = np.where(
            grad_norm[:, None] > 1e-15,
            -grad_vec / np.maximum(grad_norm[:, None], 1e-15),
            dir_a,
        )
        cos_theta = np.sum(outward_normal * dir_a, axis=1)
        cos_clamped = np.clip(cos_theta, 0.25, 1.0)
        area_weights = (r_surf * r_surf) * d_omega / cos_clamped

        # Exact Voronoi-star-shaped volume integral: V = sum_all_rays (1/3) * r_eff^3 * dOmega.
        volume_a3 = float(np.sum((t_star ** 3) * (d_omega / 3.0)))
        area_a2 = float(np.sum(area_weights))

        surf_tuples = _normalize_esp_points(surf_pts)
        rs_data = library.realspace_data(
            environment, results, surf_tuples, occ_threshold=occ_thr, lea_window_ev=lea_win
        )
        surf_rho = np.asarray(rs_data['realSpaceDensityAU'], dtype=float)
        max_iso_err = float(np.max(np.abs(surf_rho - iso)))

        esp_data = library.esp_data(environment, results, surf_tuples) if includeESP else {}

        iso_state.update({
            'points': surf_tuples,
            'owners': tuple(int(v) for v in surf_owners),
            'area_weights': tuple(float(v) for v in area_weights),
            'shell_scales': tuple(float(r_surf[i] / radii[surf_owners[i]]) for i in range(n_surf)),
            'volume_a3': volume_a3,
            'area_a2': area_a2,
            'max_isovalue_error_au': max_iso_err,
            'n_rays_buried': int(np.count_nonzero(buried)),
        })
        out: dict[str, Any] = {}
        if includeALIE:
            out['realspace'] = rs_data
        if includeESP:
            out['esp'] = esp_data
        return out

    opts_dict = selected.ToDict()
    opts_dict.pop('backend', None)
    result = _run_xtb_native(
        snapshot,
        **opts_dict,
        _keep_esp_for_hook=includeESP,
        _post_singlepoint_hook=_hook,
    )

    volume_a3 = float(iso_state['volume_a3'])
    area_a2 = float(iso_state['area_a2'])
    volume_bohr3 = float(volume_a3 * (ANGSTROM_TO_BOHR ** 3))
    area_bohr2 = float(area_a2 * (ANGSTROM_TO_BOHR ** 2))
    sphericity = float((math.pi ** (1.0 / 3.0)) * ((6.0 * volume_a3) ** (2.0 / 3.0)) / area_a2)

    grid = ESPSurfaceGrid(
        pointsAngstrom=iso_state['points'],
        atomIndices=iso_state['owners'],
        areaWeightsAngstrom2=iso_state['area_weights'],
        shellScales=iso_state['shell_scales'],
        surfaceModel=SURFACE_MODEL_DENSITY_ISOSURFACE,
        settings={
            'conf_id': conf_id,
            'isovalue_au': iso,
            'points_per_atom': n_per_atom,
            'bisection_steps': int(bisectionSteps),
            'point_count': len(iso_state['points']),
            'total_area_a2': area_a2,
            'volume_a3': volume_a3,
            'sphericity': sphericity,
            'max_isovalue_error_au': iso_state['max_isovalue_error_au'],
            'surface_definition': (
                f'exact valence electron-density isosurface rho(r) = {iso:g} e/Bohr^3 '
                'via Voronoi-partitioned radial bisection'
            ),
            'area_weight_formula': 'r^2 * (4*pi/pointsPerAtom) / clamp(n_outward . d_ray, 0.25, 1.0)',
        },
        units={
            'points_angstrom': 'angstrom',
            'atom_indices': 'index',
            'area_weights_a2': 'angstrom^2',
            'shell_scales': 'dimensionless',
            'total_area_a2': 'angstrom^2',
            'volume_a3': 'angstrom^3',
        },
        invariance=(
            'deterministic and translation invariant; Fibonacci ray axes are fixed in the Cartesian '
            'frame so arbitrary rotations introduce finite-mesh variation that decreases with pointsPerAtom'
        ),
    )

    missing: dict[str, str] = {}
    if includeESP:
        esp_surface = CalcESPSurfaceDescriptors(
            snapshot, result, surfaceGrid=grid, includeKcalMol=includeKcalMol
        )
    else:
        esp_surface = None
        missing['esp_surface'] = 'not_requested'

    if includeALIE:
        realspace_surface = CalcRealSpaceSurfaceDescriptors(
            snapshot, result, surfaceGrid=grid, densityCutoffAU=cutoff
        )
    else:
        realspace_surface = None
        missing['realspace_surface'] = 'not_requested'

    payload = {
        'surface_model': SURFACE_MODEL_DENSITY_ISOSURFACE,
        'isovalue_au': iso,
        'volume_a3': volume_a3,
        'volume_bohr3': volume_bohr3,
        'area_a2': area_a2,
        'area_bohr2': area_bohr2,
        'sphericity': sphericity,
        'n_points': len(iso_state['points']),
        'n_rays_total': int(n_rays),
        'n_rays_buried': int(iso_state['n_rays_buried']),
        'max_isovalue_error_au': float(iso_state['max_isovalue_error_au']),
        'grid': grid.ToDict(),
        'esp_surface': esp_surface,
        'realspace_surface': realspace_surface,
        'method': result.method,
        'backend': result.backend,
        'units': {
            'isovalue_au': 'e/Bohr^3',
            'volume_a3': 'angstrom^3',
            'volume_bohr3': 'Bohr^3',
            'area_a2': 'angstrom^2',
            'area_bohr2': 'Bohr^2',
            'sphericity': 'dimensionless',
            'max_isovalue_error_au': 'e/Bohr^3',
        },
        'analysis_settings': {
            'protocol': ISOSURFACE_PROTOCOL,
            'isovalue_au': iso,
            'points_per_atom': n_per_atom,
            'bisection_steps': int(bisectionSteps),
            'include_esp': includeESP,
            'include_alie': includeALIE,
            'occupation_threshold': occ_thr,
            'lea_window_ev': lea_win,
            'density_cutoff_au': cutoff,
            'include_kcal_mol': includeKcalMol,
            'conformer_id': conf_id,
            'method': result.method,
        },
        'invariance': (
            'volume, surface area, sphericity, and surface statistics are translation invariant and '
            'rotation invariant up to finite Fibonacci angular mesh discretization'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))
