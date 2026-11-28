"""Directed surface extrema: sigma-holes, pi-holes, and H-bond donor/acceptor sites.

Uses EleD's reviewed ``realspace_extension`` + ``esp_extension`` in **one**
native ``libxtb`` single-point calculation to locate exact ``isovalueAU``
(default ``0.001 e/Bohr^3``) electron-density isosurface intersections along
chemically directed vectors:

- **Sigma-holes**: along the 180-degree bond extension (and a 15-degree cone
  around it) of terminal halogen/chalcogen/pnictogen bonds ``A-X`` (default
  monovalent F/Cl/Br/I/At, or user-specified ``sigmaHoleBonds``).
- **Hydrogen-bond donors (HBD)**: along every polar ``X-H`` bond axis
  (``X in {N, O, F, S, Cl}``), yielding ``v_hbd_au`` / ``v_hbd_kcalmol`` and
  ``alie_hbd_ev`` at the exact ``0.001 a.u.`` isosurface.
- **Hydrogen-bond acceptors (HBA)**: minimum ESP (``v_hba_min_au`` /
  ``v_hba_min_kcalmol``) and minimum ALIE (``alie_hba_min_ev``) on each exposed
  heteroatom acceptor's (``N, O, F, S, P, Cl``) Voronoi ``0.001 a.u.`` surface.
- **Pi-holes**: along both ``+n`` and ``-n`` normal vectors from the centroid
  of every planar aromatic ring (or user-supplied ``rings``) to the exact
  ``0.001 a.u.`` isosurface.

Absent site classes return empty lists and ``None`` molecular summaries with
explicit ``missing`` reasons; undefined values are never encoded as zero.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from numbers import Integral
from typing import Any

import numpy as np

from .backends.native import _run_xtb_native
from .esp_descriptors import HARTREE_TO_KCAL_MOL, _require_mol, _vdw_radii
from .esp_input import _normalize_esp_points
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .realspace import (
    _fibonacci_directions,
    _validate_nonnegative_real,
    _validate_positive_real,
    _voronoi_ray_limits,
)
from .xtb import _resolve_options

DIRECTED_SURFACE_PROTOCOL = 'eled-directed-surface-extrema-v1'
_HALOGEN_Z = frozenset({9, 17, 35, 53, 85})
_HBD_DONOR_Z = frozenset({7, 8, 9, 16, 17})
_HBA_ACCEPTOR_Z = frozenset({7, 8, 9, 15, 16, 17})

__all__ = [
    'CalcDirectedSurfaceDescriptors',
]


def _orthonormal_basis(axis: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return two unit vectors orthogonal to unit vector ``axis``."""
    ref = np.array([1.0, 0.0, 0.0], dtype=float)
    if abs(float(np.dot(axis, ref))) > 0.9:
        ref = np.array([0.0, 1.0, 0.0], dtype=float)
    u = np.cross(axis, ref)
    u /= float(np.linalg.norm(u))
    v = np.cross(axis, u)
    v /= float(np.linalg.norm(v))
    return u, v


def _cone_directions(axis: np.ndarray, half_angle_deg: float = 15.0, n_ring: int = 6) -> tuple[np.ndarray, np.ndarray]:
    """Return unit directions [axis + ring at half_angle_deg] and their angles in degrees."""
    u, v = _orthonormal_basis(axis)
    theta = math.radians(half_angle_deg)
    ct, st = math.cos(theta), math.sin(theta)
    dirs = [axis]
    angles = [0.0]
    for k in range(n_ring):
        phi = 2.0 * math.pi * k / float(n_ring)
        d = ct * axis + st * (math.cos(phi) * u + math.sin(phi) * v)
        d /= float(np.linalg.norm(d))
        dirs.append(d)
        angles.append(float(half_angle_deg))
    return np.asarray(dirs, dtype=float), np.asarray(angles, dtype=float)


def _ring_plane_normal(coords: np.ndarray) -> tuple[np.ndarray, np.ndarray] | None:
    """Return (centroid, unit_normal) for a noncollinear ring via SVD."""
    centroid = np.mean(coords, axis=0)
    centered = coords - centroid
    _, svals, vh = np.linalg.svd(centered, full_matrices=False)
    if len(svals) < 3 or svals[1] < 1e-6:
        return None
    normal = vh[2]
    normal = normal / float(np.linalg.norm(normal))
    return centroid, normal


def _scale_opt(value: float | None, factor: float) -> float | None:
    return None if value is None else float(value * factor)


def CalcDirectedSurfaceDescriptors(
    mol: Any,
    *,
    isovalueAU: float = 0.001,
    includeSigmaHoles: bool = True,
    includeHydrogenBonds: bool = True,
    includePiHoles: bool = True,
    sigmaHoleBonds: Any = None,
    rings: Any = None,
    hbaPointsPerAtom: int = 60,
    bisectionSteps: int = 18,
    occThreshold: float = 1e-7,
    densityCutoffAU: float = 1e-18,
    includeKcalMol: bool = True,
    options: XTBOptions | None = None,
) -> dict[str, Any]:
    """Compute directed sigma-hole, pi-hole, and HBD/HBA extrema on the exact isosurface.

    Evaluates all directed rays in **one** native ``libxtb`` single-point
    calculation using EleD's ``realspace_extension`` and ``esp_extension``.
    """
    from rdkit import Chem

    _require_mol(mol)
    iso = _validate_positive_real(isovalueAU, 'isovalueAU')
    occ_thr = _validate_nonnegative_real(occThreshold, 'occThreshold')
    cutoff = _validate_nonnegative_real(densityCutoffAU, 'densityCutoffAU')
    for name, flag in (
        ('includeSigmaHoles', includeSigmaHoles),
        ('includeHydrogenBonds', includeHydrogenBonds),
        ('includePiHoles', includePiHoles),
        ('includeKcalMol', includeKcalMol),
    ):
        if not isinstance(flag, bool):
            raise ValueError(f'{name} must be bool')
    if (
        isinstance(hbaPointsPerAtom, (bool, np.bool_))
        or not isinstance(hbaPointsPerAtom, Integral)
        or not 12 <= hbaPointsPerAtom <= 2000
    ):
        raise ValueError('hbaPointsPerAtom must be an integer between 12 and 2000')
    if (
        isinstance(bisectionSteps, (bool, np.bool_))
        or not isinstance(bisectionSteps, Integral)
        or not 8 <= bisectionSteps <= 40
    ):
        raise ValueError('bisectionSteps must be an integer between 8 and 40')

    snapshot = Chem.Mol(mol)
    base = options or XTBOptions()
    _xyz_from_mol(snapshot, base.confId)
    charge, unpaired = _effective_charge_uhf(snapshot, base.charge, base.numUnpairedElectrons)
    conf = snapshot.GetConformer(base.confId)
    conf_id = int(conf.GetId())
    coords = np.asarray(conf.GetPositions(), dtype=float)
    radii = _vdw_radii(snapshot)
    natoms = len(coords)

    # 1. Identify sigma-hole bonds (anchor_idx -> target_idx).
    sigma_specs: list[tuple[int, int, int | None]] = []
    if includeSigmaHoles:
        if sigmaHoleBonds is not None:
            for item in sigmaHoleBonds:
                if len(item) != 2:
                    raise ValueError('each entry in sigmaHoleBonds must be a (anchor_idx, target_idx) pair')
                a_i, x_i = int(item[0]), int(item[1])
                if not (0 <= a_i < natoms and 0 <= x_i < natoms) or a_i == x_i:
                    raise ValueError('sigmaHoleBonds indices out of range')
                bond = snapshot.GetBondBetweenAtoms(a_i, x_i)
                sigma_specs.append((a_i, x_i, int(bond.GetIdx()) if bond is not None else None))
        else:
            for atom in snapshot.GetAtoms():
                if atom.GetAtomicNum() in _HALOGEN_Z and atom.GetDegree() == 1:
                    x_i = int(atom.GetIdx())
                    bond = atom.GetBonds()[0]
                    a_i = int(bond.GetOtherAtomIdx(x_i))
                    sigma_specs.append((a_i, x_i, int(bond.GetIdx())))

    # 2. Identify HBD bonds (donor_idx -> h_idx) and HBA atoms.
    hbd_specs: list[tuple[int, int, int]] = []
    hba_atoms: list[int] = []
    if includeHydrogenBonds:
        for bond in snapshot.GetBonds():
            a = bond.GetBeginAtom()
            b = bond.GetEndAtom()
            za, zb = a.GetAtomicNum(), b.GetAtomicNum()
            if za == 1 and zb in _HBD_DONOR_Z:
                hbd_specs.append((int(b.GetIdx()), int(a.GetIdx()), int(bond.GetIdx())))
            elif zb == 1 and za in _HBD_DONOR_Z:
                hbd_specs.append((int(a.GetIdx()), int(b.GetIdx()), int(bond.GetIdx())))
        for atom in snapshot.GetAtoms():
            if atom.GetAtomicNum() in _HBA_ACCEPTOR_Z:
                hba_atoms.append(int(atom.GetIdx()))

    # 3. Identify pi-hole rings.
    pi_specs: list[tuple[list[int], np.ndarray, np.ndarray]] = []
    if includePiHoles:
        raw_rings: list[list[int]] = []
        if rings is not None:
            for r in rings:
                idxs = [int(v) for v in r]
                if len(idxs) < 3 or any(not 0 <= v < natoms for v in idxs):
                    raise ValueError('each ring in rings must contain >= 3 valid atom indices')
                raw_rings.append(idxs)
        else:
            ring_info = snapshot.GetRingInfo()
            for r in ring_info.AtomRings():
                idxs = list(r)
                if 3 <= len(idxs) <= 8 and all(snapshot.GetAtomWithIdx(i).GetIsAromatic() for i in idxs):
                    raw_rings.append(idxs)
        for idxs in raw_rings:
            plane = _ring_plane_normal(coords[idxs])
            if plane is not None:
                pi_specs.append((idxs, plane[0], plane[1]))

    # Build unified ray batch: (origin, unit_dir, t_lo_init, t_hi_limit, r_vdw_hint).
    ray_origins: list[np.ndarray] = []
    ray_dirs: list[np.ndarray] = []
    ray_tcell: list[float] = []
    ray_r_hint: list[float] = []
    ray_tags: list[tuple[str, int, Any]] = []

    for s_idx, (a_i, x_i, _) in enumerate(sigma_specs):
        vec = coords[x_i] - coords[a_i]
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-6:
            raise ValueError(f'degenerate sigma-hole bond ({a_i}, {x_i})')
        axis = vec / norm
        dirs, angles = _cone_directions(axis, half_angle_deg=15.0, n_ring=6)
        for k in range(len(dirs)):
            ray_origins.append(coords[x_i])
            ray_dirs.append(dirs[k])
            ray_tcell.append(np.inf)
            ray_r_hint.append(float(radii[x_i]))
            ray_tags.append(('sigma', s_idx, float(angles[k])))

    for h_idx, (d_i, h_atom, _) in enumerate(hbd_specs):
        vec = coords[h_atom] - coords[d_i]
        norm = float(np.linalg.norm(vec))
        if norm <= 1e-6:
            raise ValueError(f'degenerate HBD bond ({d_i}, {h_atom})')
        axis = vec / norm
        ray_origins.append(coords[h_atom])
        ray_dirs.append(axis)
        ray_tcell.append(np.inf)
        ray_r_hint.append(float(radii[h_atom]))
        ray_tags.append(('hbd', h_idx, 0.0))

    if hba_atoms:
        fib = _fibonacci_directions(int(hbaPointsPerAtom))
        owners_all, origins_all, dirs_all, tcell_all = _voronoi_ray_limits(coords, fib)
        n_dir = len(fib)
        for a_pos, atom_idx in enumerate(hba_atoms):
            sl = slice(atom_idx * n_dir, (atom_idx + 1) * n_dir)
            for k in range(n_dir):
                g_i = atom_idx * n_dir + k
                ray_origins.append(origins_all[g_i])
                ray_dirs.append(dirs_all[g_i])
                ray_tcell.append(float(tcell_all[g_i]))
                ray_r_hint.append(float(radii[atom_idx]))
                ray_tags.append(('hba', a_pos, atom_idx))

    for p_idx, (_, centroid, normal) in enumerate(pi_specs):
        for sgn, face in ((1.0, 'plus'), (-1.0, 'minus')):
            ray_origins.append(centroid)
            ray_dirs.append(sgn * normal)
            ray_tcell.append(np.inf)
            ray_r_hint.append(2.0)
            ray_tags.append(('pi', p_idx, face))

    n_rays = len(ray_origins)
    selected = _resolve_options(
        replace(
            base,
            confId=conf_id,
            charge=charge,
            numUnpairedElectrons=unpaired,
            espPointsAngstrom=(),
            realSpacePointsAngstrom=(),
        )
    )

    surf_out: dict[str, Any] = {}

    if n_rays > 0:
        orig_arr = np.asarray(ray_origins, dtype=float)
        dir_arr = np.asarray(ray_dirs, dtype=float)
        tcell_arr = np.asarray(ray_tcell, dtype=float)
        hint_arr = np.asarray(ray_r_hint, dtype=float)

        def _eval_rho(library, environment, results, pts: np.ndarray) -> np.ndarray:
            if len(pts) == 0:
                return np.empty(0, dtype=float)
            d = library.realspace_data(
                environment, results, tuple(map(tuple, pts.tolist())), occ_threshold=occ_thr
            )
            return np.asarray(d['realSpaceDensityAU'], dtype=float)

        def _hook(library, environment, results, natoms_cb, nao_cb):
            library.bind_realspace_api()
            # Probe multiple inner radii so heavy-atom valence shells (which can have low density
            # right at r ~ 0.05 A due to absence of core electrons) and ring centroids are bracketed.
            t_lo = np.zeros(n_rays, dtype=float)
            rho_lo = np.zeros(n_rays, dtype=float)
            for probe_r in (0.0, 0.05, 0.35, 0.70, 1.05):
                t_cand = np.minimum(probe_r, 0.9 * tcell_arr)
                rho_cand = _eval_rho(library, environment, results, orig_arr + t_cand[:, None] * dir_arr)
                better = (rho_cand >= iso) & ((rho_lo < iso) | (t_cand > t_lo))
                t_lo = np.where(better, t_cand, t_lo)
                rho_lo = np.where(better, rho_cand, rho_lo)

            t_hi = np.minimum(tcell_arr, np.maximum(t_lo + 0.5, np.maximum(3.5, 2.2 * hint_arr)))
            rho_hi = _eval_rho(library, environment, results, orig_arr + t_hi[:, None] * dir_arr)

            for _ in range(7):
                need = (rho_hi > iso) & (t_hi < tcell_arr - 1e-12)
                if not np.any(need):
                    break
                ii = np.flatnonzero(need)
                t_hi[ii] = np.minimum(tcell_arr[ii], t_hi[ii] * 1.6)
                rho_hi[ii] = _eval_rho(
                    library, environment, results, orig_arr[ii] + t_hi[ii, None] * dir_arr[ii]
                )

            buried = (rho_hi > iso) & np.isfinite(tcell_arr) & (t_hi >= tcell_arr - 1e-12)
            active = (~buried) & (rho_lo >= iso) & (rho_hi <= iso)
            active_idx = np.flatnonzero(active)
            if len(active_idx) == 0:
                surf_out['active_idx'] = active_idx
                return {}

            lo_a = t_lo[active_idx].copy()
            hi_a = t_hi[active_idx].copy()
            oa = orig_arr[active_idx]
            da = dir_arr[active_idx]
            for _ in range(int(bisectionSteps)):
                mid = 0.5 * (lo_a + hi_a)
                rmid = _eval_rho(library, environment, results, oa + mid[:, None] * da)
                above = rmid > iso
                lo_a = np.where(above, mid, lo_a)
                hi_a = np.where(above, hi_a, mid)
            r_conv = 0.5 * (lo_a + hi_a)
            pts_conv = oa + r_conv[:, None] * da
            pts_tuples = _normalize_esp_points(pts_conv)
            rs = library.realspace_data(environment, results, pts_tuples, occ_threshold=occ_thr)
            esp = library.esp_data(environment, results, pts_tuples)
            surf_out.update({
                'active_idx': active_idx,
                'r_conv': r_conv,
                'pts_conv': pts_conv,
                'esp': np.asarray(esp['espValuesAU'], dtype=float),
                'density': np.asarray(rs['realSpaceDensityAU'], dtype=float),
                'alie': np.asarray(rs['realSpaceALIEEV'], dtype=float),
            })
            return {}

        opts_dict = selected.ToDict()
        opts_dict.pop('backend', None)
        result = _run_xtb_native(
            snapshot, **opts_dict, _keep_esp_for_hook=True, _post_singlepoint_hook=_hook
        )
        method_name, backend_name = result.method, result.backend
    else:
        method_name, backend_name = selected.method, 'native'

    # Map converged ray index -> slot in active_idx arrays.
    ray_to_slot: dict[int, int] = {}
    if 'active_idx' in surf_out:
        for slot, r_i in enumerate(surf_out['active_idx'].tolist()):
            ray_to_slot[r_i] = slot

    # Assemble sigma-holes.
    sigma_rows = []
    sigma_groups: dict[int, list[tuple[int, float]]] = {i: [] for i in range(len(sigma_specs))}
    hbd_rays: dict[int, int] = {}
    hba_groups: dict[int, list[int]] = {i: [] for i in range(len(hba_atoms))}
    pi_groups: dict[int, dict[str, int]] = {i: {} for i in range(len(pi_specs))}
    for r_i, (kind, idx_g, extra) in enumerate(ray_tags):
        if kind == 'sigma':
            sigma_groups[idx_g].append((r_i, float(extra)))
        elif kind == 'hbd':
            hbd_rays[idx_g] = r_i
        elif kind == 'hba':
            hba_groups[idx_g].append(r_i)
        elif kind == 'pi':
            pi_groups[idx_g][str(extra)] = r_i

    for s_idx, (a_i, x_i, b_idx) in enumerate(sigma_specs):
        entries = [(ray_to_slot[r_i], ang) for r_i, ang in sigma_groups[s_idx] if r_i in ray_to_slot]
        row_missing: dict[str, str] = {}
        if not entries:
            row_missing['v_sigma_max_au'] = 'isosurface_ray_buried_or_unbracketed'
            row_missing['v_sigma_axis_au'] = 'isosurface_ray_buried_or_unbracketed'
            row: dict[str, Any] = {
                'anchor_atom_index': a_i,
                'atom_index': x_i,
                'bond_index': b_idx,
                'v_sigma_axis_au': None,
                'v_sigma_max_au': None,
                'angle_from_axis_deg': None,
                'radius_a': None,
                'point_angstrom': None,
                'alie_ev': None,
            }
        else:
            axis_slots = [sl for sl, ang in entries if ang == 0.0]
            v_axis = float(surf_out['esp'][axis_slots[0]]) if axis_slots else None
            best_sl, best_ang = max(entries, key=lambda pair: float(surf_out['esp'][pair[0]]))
            v_max = float(surf_out['esp'][best_sl])
            d_best = float(surf_out['density'][best_sl])
            alie_v = float(surf_out['alie'][best_sl]) if d_best > cutoff else None
            if alie_v is None:
                row_missing['alie_ev'] = 'density_below_cutoff'
            if v_axis is None:
                row_missing['v_sigma_axis_au'] = 'axis_ray_buried'
            row = {
                'anchor_atom_index': a_i,
                'atom_index': x_i,
                'bond_index': b_idx,
                'v_sigma_axis_au': v_axis,
                'v_sigma_max_au': v_max,
                'angle_from_axis_deg': best_ang,
                'radius_a': float(surf_out['r_conv'][best_sl]),
                'point_angstrom': [float(v) for v in surf_out['pts_conv'][best_sl]],
                'alie_ev': alie_v,
            }
        if includeKcalMol:
            row['v_sigma_axis_kcalmol'] = _scale_opt(row['v_sigma_axis_au'], HARTREE_TO_KCAL_MOL)
            row['v_sigma_max_kcalmol'] = _scale_opt(row['v_sigma_max_au'], HARTREE_TO_KCAL_MOL)
        row['missing'] = row_missing
        sigma_rows.append(row)

    # Assemble HBD rows.
    hbd_rows = []
    for h_idx, (d_i, h_atom, b_idx) in enumerate(hbd_specs):
        r_i = hbd_rays.get(h_idx)
        row_missing = {}
        if r_i is None or r_i not in ray_to_slot:
            row_missing['v_hbd_au'] = 'isosurface_ray_buried_or_unbracketed'
            row = {
                'donor_atom_index': d_i,
                'hydrogen_atom_index': h_atom,
                'bond_index': b_idx,
                'v_hbd_au': None,
                'radius_a': None,
                'point_angstrom': None,
                'alie_ev': None,
            }
        else:
            sl = ray_to_slot[r_i]
            d_val = float(surf_out['density'][sl])
            alie_v = float(surf_out['alie'][sl]) if d_val > cutoff else None
            if alie_v is None:
                row_missing['alie_ev'] = 'density_below_cutoff'
            row = {
                'donor_atom_index': d_i,
                'hydrogen_atom_index': h_atom,
                'bond_index': b_idx,
                'v_hbd_au': float(surf_out['esp'][sl]),
                'radius_a': float(surf_out['r_conv'][sl]),
                'point_angstrom': [float(v) for v in surf_out['pts_conv'][sl]],
                'alie_ev': alie_v,
            }
        if includeKcalMol:
            row['v_hbd_kcalmol'] = _scale_opt(row['v_hbd_au'], HARTREE_TO_KCAL_MOL)
        row['missing'] = row_missing
        hbd_rows.append(row)

    # Assemble HBA rows.
    hba_rows = []
    for a_pos, atom_idx in enumerate(hba_atoms):
        slots = [ray_to_slot[r_i] for r_i in hba_groups[a_pos] if r_i in ray_to_slot]
        row_missing = {}
        if not slots:
            row_missing['v_hba_min_au'] = 'no_exposed_surface_points'
            row_missing['alie_min_ev'] = 'no_exposed_surface_points'
            row = {
                'atom_index': atom_idx,
                'n_points': 0,
                'v_hba_min_au': None,
                'point_angstrom': None,
                'alie_min_ev': None,
            }
        else:
            esp_a = surf_out['esp'][slots]
            best_local = int(np.argmin(esp_a))
            best_sl = slots[best_local]
            d_a = surf_out['density'][slots]
            alie_a = surf_out['alie'][slots]
            valid_a = d_a > cutoff
            alie_min = float(np.min(alie_a[valid_a])) if np.any(valid_a) else None
            if alie_min is None:
                row_missing['alie_min_ev'] = 'density_below_cutoff'
            row = {
                'atom_index': atom_idx,
                'n_points': len(slots),
                'v_hba_min_au': float(esp_a[best_local]),
                'point_angstrom': [float(v) for v in surf_out['pts_conv'][best_sl]],
                'alie_min_ev': alie_min,
            }
        if includeKcalMol:
            row['v_hba_min_kcalmol'] = _scale_opt(row['v_hba_min_au'], HARTREE_TO_KCAL_MOL)
        row['missing'] = row_missing
        hba_rows.append(row)

    # Assemble pi-hole rows.
    pi_rows = []
    for p_idx, (idxs, _, _) in enumerate(pi_specs):
        faces = pi_groups[p_idx]
        sl_p = ray_to_slot.get(faces.get('plus', -1))
        sl_m = ray_to_slot.get(faces.get('minus', -1))
        v_p = float(surf_out['esp'][sl_p]) if sl_p is not None else None
        v_m = float(surf_out['esp'][sl_m]) if sl_m is not None else None
        present = [v for v in (v_p, v_m) if v is not None]
        row_missing = {}
        if v_p is None:
            row_missing['v_pi_plus_au'] = 'isosurface_ray_unbracketed'
        if v_m is None:
            row_missing['v_pi_minus_au'] = 'isosurface_ray_unbracketed'
        v_max = max(present) if present else None
        v_min = min(present) if present else None
        if v_max is None:
            row_missing['v_pi_max_au'] = 'isosurface_ray_unbracketed'
            row_missing['v_pi_min_au'] = 'isosurface_ray_unbracketed'
        row = {
            'ring_index': p_idx,
            'atom_indices': idxs,
            'v_pi_plus_au': v_p,
            'v_pi_minus_au': v_m,
            'v_pi_max_au': v_max,
            'v_pi_min_au': v_min,
        }
        if includeKcalMol:
            row['v_pi_plus_kcalmol'] = _scale_opt(v_p, HARTREE_TO_KCAL_MOL)
            row['v_pi_minus_kcalmol'] = _scale_opt(v_m, HARTREE_TO_KCAL_MOL)
            row['v_pi_max_kcalmol'] = _scale_opt(v_max, HARTREE_TO_KCAL_MOL)
            row['v_pi_min_kcalmol'] = _scale_opt(v_min, HARTREE_TO_KCAL_MOL)
        row['missing'] = row_missing
        pi_rows.append(row)

    # Molecular-level extrema summaries with strict None + missing semantics.
    missing: dict[str, str] = {}
    sigma_vals = [r['v_sigma_max_au'] for r in sigma_rows if r['v_sigma_max_au'] is not None]
    hbd_vals = [r['v_hbd_au'] for r in hbd_rows if r['v_hbd_au'] is not None]
    hba_vals = [r['v_hba_min_au'] for r in hba_rows if r['v_hba_min_au'] is not None]
    pi_vals = [r['v_pi_max_au'] for r in pi_rows if r['v_pi_max_au'] is not None]

    v_sigma_max = max(sigma_vals) if sigma_vals else None
    v_hbd_max = max(hbd_vals) if hbd_vals else None
    v_hba_min = min(hba_vals) if hba_vals else None
    v_pi_max = max(pi_vals) if pi_vals else None

    if v_sigma_max is None:
        missing['v_sigma_max_au'] = 'not_requested' if not includeSigmaHoles else 'no_sigma_hole_sites'
    if v_hbd_max is None:
        missing['v_hbd_max_au'] = 'not_requested' if not includeHydrogenBonds else 'no_hbd_sites'
    if v_hba_min is None:
        missing['v_hba_min_au'] = 'not_requested' if not includeHydrogenBonds else 'no_hba_sites'
    if v_pi_max is None:
        missing['v_pi_max_au'] = 'not_requested' if not includePiHoles else 'no_pi_hole_rings'

    payload: dict[str, Any] = {
        'isovalue_au': iso,
        'v_sigma_max_au': v_sigma_max,
        'v_hbd_max_au': v_hbd_max,
        'v_hba_min_au': v_hba_min,
        'v_pi_max_au': v_pi_max,
        'sigma_holes': sigma_rows,
        'hbd_sites': hbd_rows,
        'hba_sites': hba_rows,
        'pi_holes': pi_rows,
        'method': method_name,
        'backend': backend_name,
        'units': {
            'isovalue_au': 'e/Bohr^3',
            'v_sigma_max_au': 'Hartree/e',
            'v_hbd_max_au': 'Hartree/e',
            'v_hba_min_au': 'Hartree/e',
            'v_pi_max_au': 'Hartree/e',
            'radius_a': 'angstrom',
            'angle_from_axis_deg': 'degree',
            'alie_ev': 'eV',
            'alie_min_ev': 'eV',
            **(
                {
                    'v_sigma_max_kcalmol': 'kcal/(mol*e)',
                    'v_hbd_max_kcalmol': 'kcal/(mol*e)',
                    'v_hba_min_kcalmol': 'kcal/(mol*e)',
                    'v_pi_max_kcalmol': 'kcal/(mol*e)',
                }
                if includeKcalMol
                else {}
            ),
        },
        'analysis_settings': {
            'protocol': DIRECTED_SURFACE_PROTOCOL,
            'isovalue_au': iso,
            'include_sigma_holes': includeSigmaHoles,
            'include_hydrogen_bonds': includeHydrogenBonds,
            'include_pi_holes': includePiHoles,
            'hba_points_per_atom': int(hbaPointsPerAtom),
            'bisection_steps': int(bisectionSteps),
            'occupation_threshold': occ_thr,
            'density_cutoff_au': cutoff,
            'include_kcal_mol': includeKcalMol,
            'conformer_id': conf_id,
        },
        'invariance': (
            'v_sigma_axis, v_hbd and v_pi_max/min are strictly rigid-motion invariant; '
            'v_hba_min and 15-degree sigma-cone max are translation invariant and rotation invariant '
            'up to finite angular mesh discretization'
        ),
        'missing': missing,
    }
    if includeKcalMol:
        payload['v_sigma_max_kcalmol'] = _scale_opt(v_sigma_max, HARTREE_TO_KCAL_MOL)
        payload['v_hbd_max_kcalmol'] = _scale_opt(v_hbd_max, HARTREE_TO_KCAL_MOL)
        payload['v_hba_min_kcalmol'] = _scale_opt(v_hba_min, HARTREE_TO_KCAL_MOL)
        payload['v_pi_max_kcalmol'] = _scale_opt(v_pi_max, HARTREE_TO_KCAL_MOL)
        for au_k, kcal_k in (
            ('v_sigma_max_au', 'v_sigma_max_kcalmol'),
            ('v_hbd_max_au', 'v_hbd_max_kcalmol'),
            ('v_hba_min_au', 'v_hba_min_kcalmol'),
            ('v_pi_max_au', 'v_pi_max_kcalmol'),
        ):
            if au_k in missing:
                missing[kcal_k] = missing[au_k]
    return json.loads(json.dumps(payload, allow_nan=False))
