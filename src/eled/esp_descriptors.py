"""Politzer-Murray surface ESP descriptors from sampled xTB potentials.

Every analysis function in this module operates on an existing xTB result with
sampled ESP values (``XTBCalculator(...).Calculate(..., espPointsAngstrom=...)``
or ``rdXTBTools.RunMol(..., espPointsAngstrom=...)``). It never invokes xTB.

:func:`GenerateESPSurfaceGrid` constructs exposed Fibonacci spheres around
RDKit periodic-table van der Waals radii scaled by ``scales`` (default
``(1.4,)``). This is an explicit scaled vdW shell surface, **not** a
0.001 a.u. quantum-mechanical electron-density isosurface.

Units are explicit in field names: coordinates are Angstrom, area weights are
Angstrom^2, primary electrostatic potentials are Hartree/e (``*_au``), primary
variances are (Hartree/e)^2 (``*_au2``), and optional kcal/mol convenience
fields use ``HARTREE_TO_KCAL_MOL = 627.5094740631`` (``*_kcalmol``,
``*_kcalmol2``). Original RDKit atom ordering is preserved, and undefined
quantities are returned as ``None`` with explicit ``missing`` reasons rather
than encoded as zero.
"""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
import json
from numbers import Integral
from typing import Any

import numpy as np

from .electronic_descriptors import _as_result
from .esp import CalcESP, _real, _vector
from .esp_input import _normalize_esp_points
from .input import _xyz_from_mol
from .mapping import _validate_mapping

from .constants import HARTREE2_TO_KCAL_MOL2, HARTREE_TO_KCAL_MOL
ESP_SURFACE_PROTOCOL = 'eled-esp-surface-v1'
SURFACE_MODEL_VDW_FIBONACCI = 'rdkit_vdw_fibonacci_shell'

__all__ = [
    'HARTREE_TO_KCAL_MOL',
    'ESPSurfaceGrid',
    'GenerateESPSurfaceGrid',
    'CalcESPSurfaceDescriptors',
]


def _require_mol(mol: Any) -> Any:
    if mol is None or not (
        hasattr(mol, 'GetNumAtoms') and hasattr(mol, 'GetConformer') and hasattr(mol, 'GetAtoms')
    ):
        raise TypeError('mol must be an RDKit Mol')
    return mol


def _validate_conf_id(confId: Any) -> int:
    if isinstance(confId, (bool, np.bool_)) or not isinstance(confId, Integral) or confId < -1:
        raise ValueError('confId must be -1 or a nonnegative integer')
    return int(confId)


def _vdw_radii(mol: Any) -> np.ndarray:
    from rdkit import Chem

    table = Chem.GetPeriodicTable()
    radii = np.array([table.GetRvdw(atom.GetAtomicNum()) for atom in mol.GetAtoms()], dtype=float)
    if not radii.size or not np.isfinite(radii).all() or np.any(radii <= 0):
        raise ValueError('positive vdW radii are required for every atom')
    return radii


@dataclass(frozen=True)
class ESPSurfaceGrid:
    """Exposed RDKit vdW Fibonacci-shell grid with atom ownership and area weights.

    Coordinates are in Angstrom and area weights are in Angstrom^2. The surface
    uses scaled RDKit periodic-table vdW shells rather than a 0.001 a.u.
    electron-density isosurface.
    """

    pointsAngstrom: tuple[tuple[float, float, float], ...]
    atomIndices: tuple[int, ...]
    areaWeightsAngstrom2: tuple[float, ...]
    shellScales: tuple[float, ...]
    surfaceModel: str
    settings: Mapping[str, Any]
    units: Mapping[str, str]
    invariance: str

    def GetPoints(self) -> np.ndarray:
        """Return an owned NumPy copy of grid coordinates in Angstrom, ``[point, xyz]``."""
        return np.array(self.pointsAngstrom, dtype=float, copy=True)

    def GetAtomIndices(self) -> np.ndarray:
        """Return an owned NumPy copy of owning RDKit atom indices, ``[point]``."""
        return np.array(self.atomIndices, dtype=int, copy=True)

    def GetAreaWeights(self) -> np.ndarray:
        """Return an owned NumPy copy of per-point surface area weights in Angstrom^2."""
        return np.array(self.areaWeightsAngstrom2, dtype=float, copy=True)

    def GetShellScales(self) -> np.ndarray:
        """Return an owned NumPy copy of per-point vdW shell scale factors."""
        return np.array(self.shellScales, dtype=float, copy=True)

    def ToDict(self) -> dict[str, Any]:
        """Return a JSON-compatible dictionary with explicit unit-suffixed keys."""
        payload = {
            'points_angstrom': [list(row) for row in self.pointsAngstrom],
            'atom_indices': list(self.atomIndices),
            'area_weights_a2': list(self.areaWeightsAngstrom2),
            'shell_scales': list(self.shellScales),
            'surface_model': self.surfaceModel,
            'settings': dict(self.settings),
            'units': dict(self.units),
            'invariance': self.invariance,
        }
        return json.loads(json.dumps(payload, allow_nan=False))


def GenerateESPSurfaceGrid(
    mol: Any,
    confId: int = -1,
    *,
    scales: Iterable[float] = (1.4,),
    pointsPerAtom: int = 100,
) -> ESPSurfaceGrid:
    """Generate exposed RDKit vdW Fibonacci-shell points with atom ownership and area weights.

    Coordinates match :func:`eled.esp.GenerateESPGrid` for the same ``confId``,
    ``scales`` and ``pointsPerAtom``. Each exposed point on atom ``i`` at shell
    scale ``s`` carries the owning RDKit atom index ``i`` and area weight
    ``4 * pi * (s * r_vdw_i)^2 / pointsPerAtom`` in Angstrom^2.

    Note
    ----
    This grid samples scaled RDKit periodic-table vdW shells, **not** a
    0.001 a.u. electron-density isosurface.
    """
    _require_mol(mol)
    conf_id = _validate_conf_id(confId)
    _xyz_from_mol(mol, conf_id)
    if (
        isinstance(pointsPerAtom, (bool, np.bool_))
        or not isinstance(pointsPerAtom, Integral)
        or not 6 <= pointsPerAtom <= 100000
    ):
        raise ValueError('pointsPerAtom must be an integer between 6 and 100000')
    if isinstance(scales, (str, bytes)) or not isinstance(scales, Iterable):
        raise ValueError('scales must be an iterable of distinct values greater than 1')
    scale_values = tuple(_real(v, 'scale', True) for v in scales)
    if not scale_values or any(v <= 1.0 for v in scale_values) or len(set(scale_values)) != len(scale_values):
        raise ValueError('scales must be distinct values greater than 1')

    conf = mol.GetConformer(conf_id)
    resolved_conf_id = int(conf.GetId())
    coords = np.asarray(conf.GetPositions(), dtype=float)
    radii = _vdw_radii(mol)

    n_per_atom = int(pointsPerAtom)
    k = np.arange(n_per_atom, dtype=float)
    z = 1.0 - 2.0 * (k + 0.5) / n_per_atom
    angle = k * (np.pi * (3.0 - np.sqrt(5.0)))
    r = np.sqrt(1.0 - z * z)
    directions = np.column_stack((r * np.cos(angle), r * np.sin(angle), z))

    blocks: list[tuple[float, float, float]] = []
    owners: list[int] = []
    areas: list[float] = []
    point_scales: list[float] = []
    for scale in scale_values:
        scaled_radii = scale * radii
        for atom, center in enumerate(coords):
            shell_radius = float(scaled_radii[atom])
            points = center + shell_radius * directions
            distances = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)
            exposed = np.all(distances >= scaled_radii[None, :] - 1e-10, axis=1)
            exposed_points = points[exposed]
            if len(exposed_points):
                point_area = float(4.0 * np.pi * (shell_radius * shell_radius) / n_per_atom)
                for row in exposed_points:
                    blocks.append((float(row[0]), float(row[1]), float(row[2])))
                    owners.append(int(atom))
                    areas.append(point_area)
                    point_scales.append(float(scale))
    if not blocks:
        raise ValueError('ESP surface grid has no exposed points')

    normalized_points = _normalize_esp_points(blocks)
    settings = {
        'conf_id': resolved_conf_id,
        'scales': list(scale_values),
        'points_per_atom': n_per_atom,
        'point_count': len(normalized_points),
        'total_area_a2': float(sum(areas)),
        'vdw_radii_source': 'rdkit.Chem.PeriodicTable.GetRvdw',
        'surface_definition': (
            'exposed RDKit periodic-table vdW Fibonacci shells '
            '(not a 0.001 a.u. electron-density isosurface)'
        ),
        'area_weight_formula': '4 * pi * (scale * r_vdw_i)^2 / pointsPerAtom',
    }
    units = {
        'points_angstrom': 'angstrom',
        'atom_indices': 'index',
        'area_weights_a2': 'angstrom^2',
        'shell_scales': 'dimensionless',
        'total_area_a2': 'angstrom^2',
    }
    invariance = (
        'deterministic and translation invariant; Fibonacci sphere axes are fixed '
        'in the Cartesian frame, so arbitrary rotations introduce finite-grid '
        'discretization variation that decreases with pointsPerAtom unless the '
        'grid is co-rotated with the molecule'
    )
    return ESPSurfaceGrid(
        pointsAngstrom=normalized_points,
        atomIndices=tuple(owners),
        areaWeightsAngstrom2=tuple(areas),
        shellScales=tuple(point_scales),
        surfaceModel=SURFACE_MODEL_VDW_FIBONACCI,
        settings=settings,
        units=units,
        invariance=invariance,
    )


def _validate_atom_indices(values: Any, npoints: int, natoms: int) -> np.ndarray:
    if isinstance(values, (str, bytes)):
        raise ValueError('atomIndices must be a sequence of integer atom indices')
    if not isinstance(values, np.ndarray) and isinstance(values, Iterable):
        values = list(values)
        if any(isinstance(v, (bool, np.bool_)) for v in values):
            raise ValueError('atomIndices must contain integer RDKit atom indices (not bool or float)')
    raw = np.asarray(values)
    if raw.ndim != 1 or len(raw) != npoints:
        raise ValueError(f'atomIndices must be a 1D integer vector of length {npoints}')
    if raw.dtype.kind not in 'iu':
        raise ValueError('atomIndices must contain integer RDKit atom indices (not bool or float)')
    indices = raw.astype(int, copy=True)
    if np.any(indices < 0) or np.any(indices >= natoms):
        raise ValueError(f'atomIndices must lie in [0, {natoms})')
    return indices


def _validate_weights(values: Any, npoints: int) -> np.ndarray:
    if isinstance(values, (str, bytes)):
        raise ValueError('weights must be a real vector')
    if not isinstance(values, np.ndarray) and isinstance(values, Iterable):
        values = list(values)
        if any(isinstance(v, (bool, np.bool_)) for v in values):
            raise ValueError('weights must be real numbers (not bool)')
    return _vector(values, 'weights', npoints)


def _assign_nearest_vdw_atoms(points: np.ndarray, coords: np.ndarray, radii: np.ndarray) -> np.ndarray:
    distances = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)
    normalized = distances / radii[None, :]
    return np.argmin(normalized, axis=1).astype(int)


def _scale_optional(value: float | None, factor: float) -> float | None:
    return None if value is None else float(value * factor)


def _esp_units(include_kcal_mol: bool) -> dict[str, str]:
    units = {
        'n_points': 'count',
        'n_points_total': 'count',
        'n_pos_points': 'count',
        'n_neg_points': 'count',
        'v_max_au': 'hartree/e',
        'v_min_au': 'hartree/e',
        'v_mean_au': 'hartree/e',
        'v_pos_mean_au': 'hartree/e',
        'v_neg_mean_au': 'hartree/e',
        'v_pos_frac': 'dimensionless',
        'v_neg_frac': 'dimensionless',
        'pi_esp_au': 'hartree/e',
        'sigma2_pos_au2': '(hartree/e)^2',
        'sigma2_neg_au2': '(hartree/e)^2',
        'sigma2_tot_au2': '(hartree/e)^2',
        'nu': 'dimensionless',
        'nu_sigma2_au2': '(hartree/e)^2',
        'sigma2_au2': '(hartree/e)^2',
    }
    if include_kcal_mol:
        units.update({
            'v_max_kcalmol': 'kcal/mol',
            'v_min_kcalmol': 'kcal/mol',
            'v_mean_kcalmol': 'kcal/mol',
            'v_pos_mean_kcalmol': 'kcal/mol',
            'v_neg_mean_kcalmol': 'kcal/mol',
            'pi_esp_kcalmol': 'kcal/mol',
            'sigma2_pos_kcalmol2': '(kcal/mol)^2',
            'sigma2_neg_kcalmol2': '(kcal/mol)^2',
            'sigma2_tot_kcalmol2': '(kcal/mol)^2',
            'nu_sigma2_kcalmol2': '(kcal/mol)^2',
            'sigma2_kcalmol2': '(kcal/mol)^2',
        })
    return units


def CalcESPSurfaceDescriptors(
    mol: Any,
    result: Any,
    *,
    atomIndices: Any = None,
    weights: Any = None,
    surfaceGrid: ESPSurfaceGrid | None = None,
    includeKcalMol: bool = True,
) -> dict[str, Any]:
    """Compute molecular Politzer-Murray GIPF and per-atom local surface ESP descriptors.

    Validates ``mol`` structure/geometry provenance against ``result`` (both
    native :class:`eled.XTBResults` and compiled ``rdXTBTools.XTBResults``) and
    reads sampled potentials via :func:`eled.esp.CalcESP` and
    ``result.GetESPPoints()``.

    Parameters
    ----------
    mol:
        Indexed RDKit molecule matching the xTB calculation provenance.
    result:
        Native ``XTBResults`` or compiled ``rdXTBTools.XTBResults`` with
        sampled ESP values.
    atomIndices:
        Optional per-point RDKit atom ownership indices. When omitted (and
        ``surfaceGrid`` is not provided), each point is assigned to the atom
        minimizing vdW-normalized distance ``d(p, atom_i) / r_vdw_i``.
    weights:
        Optional nonnegative per-point weights with a positive sum (for
        example, surface area weights from :func:`GenerateESPSurfaceGrid`).
        Defaults to uniform weights ``1.0``.
    surfaceGrid:
        Optional :class:`ESPSurfaceGrid` returned by
        :func:`GenerateESPSurfaceGrid`. Supplies both ``atomIndices`` and
        ``areaWeightsAngstrom2`` after verifying grid coordinate agreement with
        ``result.GetESPPoints()``. Cannot be combined with explicit
        ``atomIndices`` or ``weights``.
    includeKcalMol:
        When ``True`` (default), include matching ``*_kcalmol`` and
        ``*_kcalmol2`` convenience fields alongside the primary Hartree/e
        (``*_au``, ``*_au2``) fields.
    """
    _require_mol(mol)
    if not isinstance(includeKcalMol, bool):
        raise ValueError('includeKcalMol must be bool')
    normalized_result = _as_result(result)
    _validate_mapping(mol, normalized_result)

    values = np.asarray(CalcESP(normalized_result), dtype=float)
    points = np.asarray(normalized_result.GetESPPoints(), dtype=float)
    if points.ndim != 2 or points.shape != (len(values), 3) or not np.isfinite(points).all():
        raise ValueError('ESP sampling points must be a finite [npoints, 3] array matching ESP values')

    natoms = int(mol.GetNumAtoms())
    npoints = int(len(values))

    if surfaceGrid is not None:
        if not isinstance(surfaceGrid, ESPSurfaceGrid):
            raise TypeError('surfaceGrid must be an ESPSurfaceGrid instance')
        if atomIndices is not None or weights is not None:
            raise ValueError('pass either surfaceGrid or explicit atomIndices/weights, not both')
        grid_points = surfaceGrid.GetPoints()
        if grid_points.shape != points.shape or not np.allclose(grid_points, points, atol=1e-10, rtol=0.0):
            raise ValueError('surfaceGrid coordinates do not match result ESP sampling points')
        assigned = _validate_atom_indices(surfaceGrid.atomIndices, npoints, natoms)
        w = _vector(surfaceGrid.areaWeightsAngstrom2, 'weights', npoints)
        atom_assignment = 'surface_grid'
        weight_mode = 'surface_grid_area_a2'
        surface_model = surfaceGrid.surfaceModel
        vdw_radii_source: str | None = surfaceGrid.settings.get('vdw_radii_source')
    else:
        if atomIndices is not None:
            assigned = _validate_atom_indices(atomIndices, npoints, natoms)
            atom_assignment = 'explicit'
            vdw_radii_source = None
        else:
            conf_id = int(normalized_result.metadata['conformer_id'])
            coords = np.asarray(mol.GetConformer(conf_id).GetPositions(), dtype=float)
            radii = _vdw_radii(mol)
            assigned = _assign_nearest_vdw_atoms(points, coords, radii)
            atom_assignment = 'nearest_vdw_normalized'
            vdw_radii_source = 'rdkit.Chem.PeriodicTable.GetRvdw'
        if weights is None:
            w = np.ones(npoints, dtype=float)
            weight_mode = 'uniform'
        else:
            w = _validate_weights(weights, npoints)
            weight_mode = 'explicit'
        surface_model = 'sampled_esp_points'

    if np.any(w < 0) or not np.any(w > 0):
        raise ValueError('weights must be nonnegative with a positive total')

    active = w > 0
    v_active = values[active]
    w_active = w[active]
    w_total = float(np.sum(w_active))

    v_max_au = float(np.max(v_active))
    v_min_au = float(np.min(v_active))
    v_mean_au = float(np.sum(w_active * v_active) / w_total)
    pi_esp_au = float(np.sum(w_active * np.abs(v_active - v_mean_au)) / w_total)

    pos_mask = v_active > 0.0
    neg_mask = v_active < 0.0
    n_pos = int(np.count_nonzero(pos_mask))
    n_neg = int(np.count_nonzero(neg_mask))
    w_pos = float(np.sum(w_active[pos_mask])) if n_pos else 0.0
    w_neg = float(np.sum(w_active[neg_mask])) if n_neg else 0.0
    v_pos_frac = float(w_pos / w_total)
    v_neg_frac = float(w_neg / w_total)

    mol_missing: dict[str, str] = {}
    if n_pos:
        v_pos = v_active[pos_mask]
        wp = w_active[pos_mask]
        v_pos_mean_au: float | None = float(np.sum(wp * v_pos) / w_pos)
        sigma2_pos_au2: float | None = float(np.sum(wp * (v_pos - v_pos_mean_au) ** 2) / w_pos)
    else:
        v_pos_mean_au = None
        sigma2_pos_au2 = None
        mol_missing['v_pos_mean_au'] = 'no_positive_surface_points'
        mol_missing['sigma2_pos_au2'] = 'no_positive_surface_points'

    if n_neg:
        v_neg = v_active[neg_mask]
        wn = w_active[neg_mask]
        v_neg_mean_au: float | None = float(np.sum(wn * v_neg) / w_neg)
        sigma2_neg_au2: float | None = float(np.sum(wn * (v_neg - v_neg_mean_au) ** 2) / w_neg)
    else:
        v_neg_mean_au = None
        sigma2_neg_au2 = None
        mol_missing['v_neg_mean_au'] = 'no_negative_surface_points'
        mol_missing['sigma2_neg_au2'] = 'no_negative_surface_points'

    if sigma2_pos_au2 is not None and sigma2_neg_au2 is not None:
        sigma2_tot_au2: float | None = float(sigma2_pos_au2 + sigma2_neg_au2)
        if sigma2_tot_au2 > 0.0:
            nu: float | None = float((sigma2_pos_au2 * sigma2_neg_au2) / (sigma2_tot_au2 * sigma2_tot_au2))
            nu_sigma2_au2: float | None = float((sigma2_pos_au2 * sigma2_neg_au2) / sigma2_tot_au2)
        else:
            nu = None
            nu_sigma2_au2 = None
            mol_missing['nu'] = 'zero_total_variance'
            mol_missing['nu_sigma2_au2'] = 'zero_total_variance'
    elif sigma2_pos_au2 is not None or sigma2_neg_au2 is not None:
        sigma2_tot_au2 = float(sigma2_pos_au2 if sigma2_pos_au2 is not None else sigma2_neg_au2)
        nu = None
        nu_sigma2_au2 = None
        mol_missing['nu'] = 'single_sign_surface'
        mol_missing['nu_sigma2_au2'] = 'single_sign_surface'
    else:
        sigma2_tot_au2 = None
        nu = None
        nu_sigma2_au2 = None
        mol_missing['sigma2_tot_au2'] = 'no_nonzero_surface_points'
        mol_missing['nu'] = 'no_nonzero_surface_points'
        mol_missing['nu_sigma2_au2'] = 'no_nonzero_surface_points'

    output: dict[str, Any] = {
        'n_points': int(len(v_active)),
        'n_points_total': npoints,
        'n_pos_points': n_pos,
        'n_neg_points': n_neg,
        'v_max_au': v_max_au,
        'v_min_au': v_min_au,
        'v_mean_au': v_mean_au,
        'v_pos_mean_au': v_pos_mean_au,
        'v_neg_mean_au': v_neg_mean_au,
        'v_pos_frac': v_pos_frac,
        'v_neg_frac': v_neg_frac,
        'pi_esp_au': pi_esp_au,
        'sigma2_pos_au2': sigma2_pos_au2,
        'sigma2_neg_au2': sigma2_neg_au2,
        'sigma2_tot_au2': sigma2_tot_au2,
        'nu': nu,
        'nu_sigma2_au2': nu_sigma2_au2,
    }

    if includeKcalMol:
        output.update({
            'v_max_kcalmol': float(v_max_au * HARTREE_TO_KCAL_MOL),
            'v_min_kcalmol': float(v_min_au * HARTREE_TO_KCAL_MOL),
            'v_mean_kcalmol': float(v_mean_au * HARTREE_TO_KCAL_MOL),
            'v_pos_mean_kcalmol': _scale_optional(v_pos_mean_au, HARTREE_TO_KCAL_MOL),
            'v_neg_mean_kcalmol': _scale_optional(v_neg_mean_au, HARTREE_TO_KCAL_MOL),
            'pi_esp_kcalmol': float(pi_esp_au * HARTREE_TO_KCAL_MOL),
            'sigma2_pos_kcalmol2': _scale_optional(sigma2_pos_au2, HARTREE2_TO_KCAL_MOL2),
            'sigma2_neg_kcalmol2': _scale_optional(sigma2_neg_au2, HARTREE2_TO_KCAL_MOL2),
            'sigma2_tot_kcalmol2': _scale_optional(sigma2_tot_au2, HARTREE2_TO_KCAL_MOL2),
            'nu_sigma2_kcalmol2': _scale_optional(nu_sigma2_au2, HARTREE2_TO_KCAL_MOL2),
        })
        for au_key, kcal_key in (
            ('v_pos_mean_au', 'v_pos_mean_kcalmol'),
            ('v_neg_mean_au', 'v_neg_mean_kcalmol'),
            ('sigma2_pos_au2', 'sigma2_pos_kcalmol2'),
            ('sigma2_neg_au2', 'sigma2_neg_kcalmol2'),
            ('sigma2_tot_au2', 'sigma2_tot_kcalmol2'),
            ('nu_sigma2_au2', 'nu_sigma2_kcalmol2'),
        ):
            if au_key in mol_missing:
                mol_missing[kcal_key] = mol_missing[au_key]

    atoms_rows: list[dict[str, Any]] = []
    atom_missing: dict[str, dict[str, str]] = {}
    atom_value_keys = (
        'v_min_au',
        'v_max_au',
        'v_mean_au',
        'v_pos_mean_au',
        'v_neg_mean_au',
        'sigma2_au2',
    )
    atom_kcal_pairs = (
        ('v_min_au', 'v_min_kcalmol', HARTREE_TO_KCAL_MOL),
        ('v_max_au', 'v_max_kcalmol', HARTREE_TO_KCAL_MOL),
        ('v_mean_au', 'v_mean_kcalmol', HARTREE_TO_KCAL_MOL),
        ('v_pos_mean_au', 'v_pos_mean_kcalmol', HARTREE_TO_KCAL_MOL),
        ('v_neg_mean_au', 'v_neg_mean_kcalmol', HARTREE_TO_KCAL_MOL),
        ('sigma2_au2', 'sigma2_kcalmol2', HARTREE2_TO_KCAL_MOL2),
    )

    for atom_idx in range(natoms):
        assigned_mask = assigned == atom_idx
        active_atom_mask = assigned_mask & active
        n_atom_active = int(np.count_nonzero(active_atom_mask))
        row_missing: dict[str, str] = {}
        if n_atom_active == 0:
            reason = (
                'zero_surface_weight_assigned'
                if np.any(assigned_mask)
                else 'no_surface_points_assigned'
            )
            row: dict[str, Any] = {'atom_index': atom_idx, 'n_points': 0}
            for key in atom_value_keys:
                row[key] = None
                row_missing[key] = reason
            if includeKcalMol:
                for _, kcal_key, _ in atom_kcal_pairs:
                    row[kcal_key] = None
                    row_missing[kcal_key] = reason
        else:
            va = values[active_atom_mask]
            wa = w[active_atom_mask]
            wa_sum = float(np.sum(wa))
            a_min = float(np.min(va))
            a_max = float(np.max(va))
            a_mean = float(np.sum(wa * va) / wa_sum)
            a_sigma2 = float(np.sum(wa * (va - a_mean) ** 2) / wa_sum)

            a_pos_mask = va > 0.0
            a_neg_mask = va < 0.0
            if np.any(a_pos_mask):
                a_pos_mean: float | None = float(
                    np.sum(wa[a_pos_mask] * va[a_pos_mask]) / np.sum(wa[a_pos_mask])
                )
            else:
                a_pos_mean = None
                row_missing['v_pos_mean_au'] = 'no_positive_surface_points'
            if np.any(a_neg_mask):
                a_neg_mean: float | None = float(
                    np.sum(wa[a_neg_mask] * va[a_neg_mask]) / np.sum(wa[a_neg_mask])
                )
            else:
                a_neg_mean = None
                row_missing['v_neg_mean_au'] = 'no_negative_surface_points'

            row = {
                'atom_index': atom_idx,
                'n_points': n_atom_active,
                'v_min_au': a_min,
                'v_max_au': a_max,
                'v_mean_au': a_mean,
                'v_pos_mean_au': a_pos_mean,
                'v_neg_mean_au': a_neg_mean,
                'sigma2_au2': a_sigma2,
            }
            if includeKcalMol:
                for au_key, kcal_key, factor in atom_kcal_pairs:
                    row[kcal_key] = _scale_optional(row[au_key], factor)
                    if au_key in row_missing:
                        row_missing[kcal_key] = row_missing[au_key]

        row['missing'] = row_missing
        if row_missing:
            atom_missing[str(atom_idx)] = dict(row_missing)
        atoms_rows.append(row)

    missing: dict[str, Any] = dict(mol_missing)
    if atom_missing:
        missing['atoms'] = atom_missing

    output['atoms'] = atoms_rows
    output['units'] = _esp_units(includeKcalMol)
    output['invariance'] = (
        'invariant to rigid translation/rotation when ESP points co-transform with the molecule, '
        'invariant to ESP point permutation, and invariant to RDKit atom renumbering '
        '(with remapped atom_index / atomIndices)'
    )
    output['analysis_settings'] = {
        'protocol': ESP_SURFACE_PROTOCOL,
        'surface_model': surface_model,
        'atom_assignment': atom_assignment,
        'vdw_radii_source': vdw_radii_source,
        'weight_mode': weight_mode,
        'total_weight': w_total,
        'sign_convention': (
            'positive: V > 0; negative: V < 0; zero-potential points are excluded '
            'from sign-partitioned means and variances'
        ),
        'variance_convention': (
            'weighted population variance (ddof=0); sigma2_tot = sigma2_pos + sigma2_neg, '
            'reducing to the present-sign variance on a single-sign surface while the '
            'absent-sign variance and nu remain None'
        ),
        'atom_variance_convention': 'weighted population variance (ddof=0) of V around atom v_mean',
        'include_kcal_mol': includeKcalMol,
        'hartree_to_kcal_mol': HARTREE_TO_KCAL_MOL,
        'point_count': npoints,
        'active_point_count': int(len(v_active)),
        'method': normalized_result.method,
        'backend': normalized_result.backend,
        'conformer_id': int(normalized_result.metadata['conformer_id']),
    }
    output['missing'] = missing
    return json.loads(json.dumps(output, allow_nan=False))
