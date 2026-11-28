"""Conformer-specific steric, plane and ordered-ring descriptors.

All indices are original RDKit indices. Coordinates and radii are Angstrom.
Functions analyze the supplied conformer and return owned JSON-compatible data.
"""
from __future__ import annotations

from collections.abc import Mapping
from numbers import Integral, Real

import numpy as np
from rdkit import Chem, rdBase

from .electronic_descriptors import _indices
from .input import _input_metadata, _xyz_from_mol
from .reports import DESCRIPTOR_VERSION, SCHEMA_VERSION
from .selections import NormalizeSelections


HOMA_PARAMETERS = {
    (6, 6): (1.388, 257.7),
    (6, 7): (1.334, 93.52),
    (6, 8): (1.265, 157.38),
    (6, 16): (1.677, 94.09),
}
HOMA_PARAMETER_SOURCE = 'Krygowski 1993, J. Chem. Inf. Comput. Sci. 33, 70-78'


def _positive(value, name, *, integer=False):
    kind = Integral if integer else Real
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, kind) or not np.isfinite(value) or value <= 0:
        raise ValueError(f'{name} must be a positive {"integer" if integer else "real"}')
    return int(value) if integer else float(value)


def _geometry(mol, confId):
    if isinstance(confId, (bool, np.bool_)) or not isinstance(confId, Integral) or confId < -1:
        raise ValueError('confId must be -1 or a nonnegative integer')
    xyz = _xyz_from_mol(mol, int(confId))
    conf = mol.GetConformer(int(confId))
    metadata = _input_metadata(
        mol,
        int(confId),
        Chem.GetFormalCharge(mol),
        sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()),
        xyz,
    )
    return np.asarray(conf.GetPositions()), metadata


def _report(kind, metadata, selection, settings, values, units, missing=None):
    return dict(
        schema_version=SCHEMA_VERSION,
        descriptor_version=DESCRIPTOR_VERSION,
        report_type=kind,
        protocol_version='1.0',
        selection=selection,
        provenance=metadata,
        analysis_settings=settings,
        values=values,
        units=units,
        missing=missing or {},
    )


def AnalyzeLocalGeometry(
    mol,
    *,
    selections=None,
    confId=-1,
    sterimol=None,
    buriedVolume=None,
    planes=None,
    homa=None,
    puckering=None,
    pyramidalization=None,
    coneAngle=None,
):
    """Run selected local geometry descriptors using one canonical selection set."""
    normalized = NormalizeSelections(mol, selections)
    out = {
        'schema_version': SCHEMA_VERSION,
        'descriptor_version': DESCRIPTOR_VERSION,
        'report_type': 'local_geometry',
        'protocol_version': '1.0',
        'selections': normalized.ToDict(),
        'reports': {},
        'errors': {},
    }

    def run(name, function, *args, **kwargs):
        try:
            out['reports'][name] = function(mol, *args, confId=confId, **kwargs)
        except Exception as exc:
            out['errors'][name] = {'type': type(exc).__name__, 'message': str(exc)}

    for name, options in (sterimol or {}).items():
        if name not in normalized.atoms:
            raise ValueError(f'sterimol request {name!r} has no atom selection')
        options = dict(options)
        axis = options.pop('axis')
        run(f'sterimol:{name}', CalcSterimol, normalized.atoms[name], axis=axis, **options)
    for name, options in (buriedVolume or {}).items():
        if name not in normalized.atoms:
            raise ValueError(f'buried-volume request {name!r} has no atom selection')
        options = dict(options)
        center = options.pop('center')
        run(f'buried_volume:{name}', CalcBuriedVolume, normalized.atoms[name], center=center, **options)
    for name, options in (planes or {}).items():
        if name not in normalized.atoms:
            raise ValueError(f'plane request {name!r} has no atom selection')
        options = dict(options)
        other_name = options.pop('otherAtomSelection', None)
        other = None if other_name is None else normalized.atoms[other_name]
        run(f'plane:{name}', CalcPlaneDescriptors, normalized.atoms[name], otherAtomIndices=other, **options)
    for name, options in (homa or {}).items():
        if name not in normalized.rings:
            raise ValueError(f'HOMA request {name!r} has no ring selection')
        run(f'homa:{name}', CalcHOMA, normalized.rings[name], **dict(options))
    for name, options in (puckering or {}).items():
        if name not in normalized.rings:
            raise ValueError(f'puckering request {name!r} has no ring selection')
        run(f'puckering:{name}', CalcRingPuckering, normalized.rings[name], **dict(options))
    for name, options in (pyramidalization or {}).items():
        options = dict(options)
        atom_idx = options.pop('atomIdx')
        run(f'pyramidalization:{name}', CalcPyramidalization, atom_idx, **options)
    for name, options in (coneAngle or {}).items():
        if name not in normalized.atoms:
            raise ValueError(f'cone-angle request {name!r} has no atom selection')
        options = dict(options)
        center = options.pop('center')
        run(f'cone_angle:{name}', CalcConeAngle, normalized.atoms[name], center=center, **options)
    return out


def RunLocalGeometryBatch(
    molecules,
    selections,
    *,
    confIds=None,
    returnErrors=True,
    sterimol=None,
    buriedVolume=None,
    planes=None,
    homa=None,
    puckering=None,
    pyramidalization=None,
    coneAngle=None,
):
    """Run ``AnalyzeLocalGeometry`` over molecules with stable ordered rows."""
    molecules = list(molecules)
    if confIds is None:
        confIds = [-1] * len(molecules)
    else:
        confIds = list(confIds)
    if len(confIds) != len(molecules):
        raise ValueError('confIds must match molecules')
    rows = []
    for mol, conf_id in zip(molecules, confIds):
        try:
            rows.append(
                AnalyzeLocalGeometry(
                    mol,
                    selections=selections,
                    confId=conf_id,
                    sterimol=sterimol,
                    buriedVolume=buriedVolume,
                    planes=planes,
                    homa=homa,
                    puckering=puckering,
                    pyramidalization=pyramidalization,
                    coneAngle=coneAngle,
                )
            )
        except Exception as exc:
            if not returnErrors:
                raise
            rows.append({
                'schema_version': SCHEMA_VERSION,
                'descriptor_version': DESCRIPTOR_VERSION,
                'report_type': 'local_geometry',
                'protocol_version': '1.0',
                'selections': None,
                'reports': {},
                'errors': {'molecule': {'type': type(exc).__name__, 'message': str(exc)}},
            })
    return rows


def _radii(mol, radiiAngstrom, radiusScale):
    scale = _positive(radiusScale, 'radiusScale')
    if radiiAngstrom is None:
        values = [Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()]
        source = f'RDKit {rdBase.rdkitVersion} van der Waals radii'
    else:
        values = radiiAngstrom
        source = 'user_supplied_per_atom'
    raw = np.asarray(values)
    if raw.dtype.kind not in 'iuf' or raw.shape != (mol.GetNumAtoms(),):
        raise ValueError('radiiAngstrom must contain one real radius per original atom')
    values = np.asarray(raw, dtype=float) * scale
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError('atomic radii must be finite and positive')
    return values, dict(radius_source=source, radius_scale=scale, effective_radii_angstrom=values.tolist())


def CalcSterimol(
    mol,
    atomIndices,
    *,
    axis,
    confId=-1,
    radiiAngstrom=None,
    radiusScale=1.0,
    lCorrectionAngstrom=0.4,
    angularSamples=7200,
):
    """Sterimol L/B1/B5 from atom spheres around directed attachment axis."""
    coordinates, metadata = _geometry(mol, confId)
    selected = _indices(atomIndices, mol.GetNumAtoms(), 'atomIndices')
    directed = _indices(axis, mol.GetNumAtoms(), 'axis')
    if len(directed) != 2 or directed[0] in selected or directed[1] not in selected:
        raise ValueError('axis must be (excluded base, selected attached atom)')
    if mol.GetBondBetweenAtoms(int(directed[0]), int(directed[1])) is None:
        raise ValueError('axis atoms must share a bond')
    count = _positive(angularSamples, 'angularSamples', integer=True)
    if not 16 <= count <= 1000000:
        raise ValueError('angularSamples must be between 16 and 1000000')
    if (isinstance(lCorrectionAngstrom, (bool, np.bool_))
            or not isinstance(lCorrectionAngstrom, Real)
            or not np.isfinite(lCorrectionAngstrom)):
        raise ValueError('lCorrectionAngstrom must be finite and real')
    radii, settings = _radii(mol, radiiAngstrom, radiusScale)
    origin = coordinates[directed[0]]
    direction = coordinates[directed[1]] - origin
    norm = np.linalg.norm(direction)
    if norm < 1e-12:
        raise ValueError('axis atoms must have distinct coordinates')
    direction /= norm
    relative = coordinates[selected] - origin
    axial = relative @ direction
    transverse = relative - axial[:, None] * direction
    distances = np.linalg.norm(transverse, axis=1)
    if distances.max() > 1e-12:
        first = transverse[np.argmax(distances)] / distances.max()
    else:
        basis = np.eye(3)[np.argmin(np.abs(direction))]
        first = np.cross(direction, basis)
        first /= np.linalg.norm(first)
    second = np.cross(direction, first)
    x, y = transverse @ first, transverse @ second
    b1 = float('inf')
    for start in range(0, count, 256):
        angles = np.arange(start, min(start + 256, count)) * (2 * np.pi / count)
        support = x[:, None] * np.cos(angles) + y[:, None] * np.sin(angles) + radii[selected, None]
        b1 = min(b1, float(support.max(axis=0).min()))
    settings.update(
        angular_samples=count,
        l_correction_angstrom=float(lCorrectionAngstrom),
        b1_error_bound_angstrom=float(distances.max() * np.pi / count),
    )
    return _report(
        'sterimol',
        metadata,
        dict(atom_indices=selected.tolist(), axis=directed.tolist()),
        settings,
        dict(
            l_a=float(np.max(axial + radii[selected]) + lCorrectionAngstrom),
            b1_a=b1,
            b5_a=float(np.max(distances + radii[selected])),
        ),
        dict(l_a='angstrom', b1_a='angstrom', b5_a='angstrom'),
    )


def CalcBuriedVolume(
    mol,
    atomIndices,
    *,
    center,
    confId=-1,
    sphereRadiusAngstrom=3.5,
    spacingAngstrom=0.1,
    radiiAngstrom=None,
    radiusScale=1.0,
    zAxisAtom=None,
    xzPlaneAtom=None,
):
    """Percent of a center sphere occupied by the union of selected atom spheres.

    When ``zAxisAtom`` is provided (atom index defining the +z axis from
    ``center`` toward ``zAxisAtom``; optional ``xzPlaneAtom`` fixes the +x
    half-plane), also computes SambVca-style four-quadrant buried volumes
    (``vbur_q1_percent`` .. ``vbur_q4_percent`` in the xy plane perpendicular
    to z) and hemisphere breakdown (``vbur_z_pos_percent``, ``vbur_z_neg_percent``).
    """
    coordinates, metadata = _geometry(mol, confId)
    selected = _indices(atomIndices, mol.GetNumAtoms(), 'atomIndices')
    central = int(_indices([center], mol.GetNumAtoms(), 'center')[0])
    if central in selected:
        raise ValueError('exclude the center atom from atomIndices')
    radius = _positive(sphereRadiusAngstrom, 'sphereRadiusAngstrom')
    spacing = _positive(spacingAngstrom, 'spacingAngstrom')
    cells = np.ceil(2 * radius / spacing)
    if not 4 <= cells <= 256:
        raise ValueError('choose spacing giving 4..256 grid cells per sphere diameter')
    n = int(cells)
    radii, settings = _radii(mol, radiiAngstrom, radiusScale)

    frame_axes = None
    z_idx = None
    xz_idx = None
    if zAxisAtom is not None:
        z_idx = int(_indices([zAxisAtom], mol.GetNumAtoms(), 'zAxisAtom')[0])
        if z_idx == central:
            raise ValueError('zAxisAtom must differ from center')
        z_vec = coordinates[z_idx] - coordinates[central]
        z_norm = float(np.linalg.norm(z_vec))
        if z_norm < 1e-10:
            raise ValueError('zAxisAtom must not coincide with center')
        z_hat = z_vec / z_norm
        if xzPlaneAtom is not None:
            xz_idx = int(_indices([xzPlaneAtom], mol.GetNumAtoms(), 'xzPlaneAtom')[0])
            if xz_idx in (central, z_idx):
                raise ValueError('xzPlaneAtom must differ from center and zAxisAtom')
            ref = coordinates[xz_idx] - coordinates[central]
            x_vec = ref - float(np.dot(ref, z_hat)) * z_hat
            x_norm = float(np.linalg.norm(x_vec))
            if x_norm < 1e-10:
                raise ValueError('xzPlaneAtom must not be collinear with center-zAxisAtom')
            x_hat = x_vec / x_norm
        else:
            basis = np.eye(3)[int(np.argmin(np.abs(z_hat)))]
            x_vec = basis - float(np.dot(basis, z_hat)) * z_hat
            x_hat = x_vec / float(np.linalg.norm(x_vec))
        y_hat = np.cross(z_hat, x_hat)
        frame_axes = np.vstack((x_hat, y_hat, z_hat))  # rows: x_hat, y_hat, z_hat
    elif xzPlaneAtom is not None:
        raise ValueError('xzPlaneAtom requires zAxisAtom')

    step = 2 * radius / n
    grid = (np.arange(n) + 0.5) * step - radius
    y, z = np.meshgrid(grid, grid, indexing='ij')
    yz = np.column_stack((y.ravel(), z.ravel()))
    centers = coordinates[selected] - coordinates[central]
    active = np.linalg.norm(centers, axis=1) < radius + radii[selected]
    centers, atom_radii = centers[active], radii[selected][active]
    total, occupied = 0, 0
    q_tot = np.zeros(4, dtype=float)
    q_occ = np.zeros(4, dtype=float)
    z_tot = np.zeros(2, dtype=float)
    z_occ = np.zeros(2, dtype=float)
    for x in grid:
        plane = yz[np.sum(yz ** 2, axis=1) + x * x <= radius * radius]
        points = np.column_stack((np.full(len(plane), x), plane))
        inside = np.zeros(len(points), dtype=bool)
        for position, atom_radius in zip(centers, atom_radii):
            inside |= np.sum((points - position) ** 2, axis=1) <= atom_radius ** 2
        total += len(points)
        occupied += int(inside.sum())
        if frame_axes is not None and len(points):
            local = points @ frame_axes.T
            lx, ly, lz = local[:, 0], local[:, 1], local[:, 2]
            occ_f = inside.astype(float)
            # Boundary points exactly on a coordinate plane split 0.5/0.5 so
            # odd and even grids keep exact quadrant/hemisphere symmetry.
            wx_p = np.where(lx > 1e-12, 1.0, np.where(lx < -1e-12, 0.0, 0.5))
            wx_n = 1.0 - wx_p
            wy_p = np.where(ly > 1e-12, 1.0, np.where(ly < -1e-12, 0.0, 0.5))
            wy_n = 1.0 - wy_p
            wz_p = np.where(lz > 1e-12, 1.0, np.where(lz < -1e-12, 0.0, 0.5))
            wz_n = 1.0 - wz_p
            weights = (wx_p * wy_p, wx_n * wy_p, wx_n * wy_n, wx_p * wy_n)
            for qi, qw in enumerate(weights):
                q_tot[qi] += float(qw.sum())
                q_occ[qi] += float((qw * occ_f).sum())
            z_tot[0] += float(wz_p.sum())
            z_occ[0] += float((wz_p * occ_f).sum())
            z_tot[1] += float(wz_n.sum())
            z_occ[1] += float((wz_n * occ_f).sum())

    settings.update(
        sphere_radius_angstrom=radius,
        requested_spacing_angstrom=spacing,
        effective_spacing_angstrom=step,
        grid_points=total,
        occupied_grid_points=occupied,
        quadrature='centered_cartesian_midpoint',
    )
    selection = dict(atom_indices=selected.tolist(), center=central)
    values = dict(vbur_percent=100.0 * occupied / total)
    units = dict(vbur_percent='percent')
    if frame_axes is not None:
        selection['z_axis_atom'] = z_idx
        if xz_idx is not None:
            selection['xz_plane_atom'] = xz_idx
        for qi in range(4):
            key = f'vbur_q{qi + 1}_percent'
            values[key] = float(100.0 * q_occ[qi] / q_tot[qi]) if q_tot[qi] > 0 else None
            units[key] = 'percent'
        values['vbur_z_pos_percent'] = float(100.0 * z_occ[0] / z_tot[0]) if z_tot[0] > 0 else None
        values['vbur_z_neg_percent'] = float(100.0 * z_occ[1] / z_tot[1]) if z_tot[1] > 0 else None
        units['vbur_z_pos_percent'] = 'percent'
        units['vbur_z_neg_percent'] = 'percent'
    return _report(
        'buried_volume',
        metadata,
        selection,
        settings,
        values,
        units,
    )


def _plane(points):
    if len(points) < 3:
        raise ValueError('a plane requires at least three atoms')
    center = points.mean(axis=0)
    shifted = points - center
    eigenvalues, vectors = np.linalg.eigh(shifted.T @ shifted)
    scale = max(float(eigenvalues[-1]), 1e-24)
    if eigenvalues[1] <= scale * 1e-12 or eigenvalues[1] - eigenvalues[0] <= scale * 1e-12:
        return center, None, None
    normal = vectors[:, 0]
    return center, normal, shifted @ normal


def CalcPlaneDescriptors(mol, atomIndices, *, otherAtomIndices=None, confId=-1):
    """Unweighted least-squares plane distances and acute angle between planes."""
    coordinates, metadata = _geometry(mol, confId)
    selected = _indices(atomIndices, mol.GetNumAtoms(), 'atomIndices')
    centroid, normal, distances = _plane(coordinates[selected])
    values = dict(plane_rms_a=None, plane_max_a=None)
    units = dict(plane_rms_a='angstrom', plane_max_a='angstrom')
    missing = {}
    if normal is None:
        missing.update(plane_rms_a='plane_normal_not_unique', plane_max_a='plane_normal_not_unique')
    else:
        values.update(
            plane_rms_a=float(np.sqrt(np.mean(distances ** 2))),
            plane_max_a=float(np.max(np.abs(distances))),
        )
    selection = dict(atom_indices=selected.tolist())
    if otherAtomIndices is not None:
        other = _indices(otherAtomIndices, mol.GetNumAtoms(), 'otherAtomIndices')
        other_center, other_normal, _ = _plane(coordinates[other])
        selection['other_atom_indices'] = other.tolist()
        values.update(centroid_distance_a=float(np.linalg.norm(centroid - other_center)), plane_angle_deg=None)
        units.update(centroid_distance_a='angstrom', plane_angle_deg='degree')
        if normal is None or other_normal is None:
            missing['plane_angle_deg'] = 'plane_normal_not_unique'
        else:
            values['plane_angle_deg'] = float(np.degrees(np.arccos(np.clip(abs(normal @ other_normal), 0, 1))))
    return _report('plane_geometry', metadata, selection, dict(weighting='uniform'), values, units, missing)


def _ring(mol, ring):
    indices = _indices(ring, mol.GetNumAtoms(), 'ringAtomIndices')
    if len(indices) < 3:
        raise ValueError('a ring requires at least three distinct atoms')
    for first, second in zip(indices, np.roll(indices, -1)):
        if mol.GetBondBetweenAtoms(int(first), int(second)) is None:
            raise ValueError('ringAtomIndices must form an ordered bonded closed path')
    return indices


def CalcHOMA(mol, ringAtomIndices, *, confId=-1, parameters=None, parameterSource=None):
    """HOMA=1-mean(alpha_AB*(R_AB-Ropt_AB)^2) on an ordered ring."""
    coordinates, metadata = _geometry(mol, confId)
    indices = _ring(mol, ringAtomIndices)
    if parameters is None:
        parameters = HOMA_PARAMETERS
        parameterSource = HOMA_PARAMETER_SOURCE
    elif not isinstance(parameterSource, str) or not parameterSource.strip():
        raise ValueError('custom HOMA parameters require parameterSource')
    if not isinstance(parameters, Mapping):
        raise ValueError('parameters must map element pairs to (optimal_length, alpha)')
    table = {}
    for pair, row in parameters.items():
        if (not isinstance(pair, tuple) or len(pair) != 2
                or any(isinstance(z, (bool, np.bool_)) or not isinstance(z, Integral)
                       or not 1 <= z <= 118 for z in pair)):
            raise ValueError('parameter keys must be pairs of atomic numbers')
        key = tuple(sorted(pair))
        if key in table or not isinstance(row, (tuple, list, np.ndarray)) or np.shape(row) != (2,):
            raise ValueError('duplicate element pair or invalid HOMA parameter row')
        table[key] = (_positive(row[0], 'optimal_length'), _positive(row[1], 'alpha'))
    penalties, bonds, unsupported = [], [], []
    for first, second in zip(indices, np.roll(indices, -1)):
        length = float(np.linalg.norm(coordinates[first] - coordinates[second]))
        if length <= 1e-12:
            raise ValueError('ring bond length must be positive')
        pair = tuple(
            sorted((mol.GetAtomWithIdx(int(first)).GetAtomicNum(), mol.GetAtomWithIdx(int(second)).GetAtomicNum()))
        )
        bonds.append(dict(begin_atom_index=int(first), end_atom_index=int(second), length_a=length))
        if pair not in table:
            unsupported.append(list(pair))
        else:
            optimum, alpha = table[pair]
            penalties.append(alpha * (length - optimum) ** 2)
    settings = dict(
        parameter_source=parameterSource,
        parameters=[
            dict(atomic_numbers=list(pair), optimal_length_a=r, alpha_inv_a2=a)
            for pair, (r, a) in sorted(table.items())
        ],
    )
    missing = dict(homa='unsupported_element_pair') if unsupported else {}
    values = dict(homa=None if unsupported else float(1 - np.mean(penalties)))
    report = _report(
        'homa',
        metadata,
        dict(ring_atom_indices=indices.tolist()),
        settings,
        values,
        dict(homa='dimensionless'),
        missing,
    )
    report.update(bonds=bonds, bond_units=dict(length_a='angstrom'), unsupported_element_pairs=unsupported)
    return report


def CalcRingPuckering(mol, ringAtomIndices, *, confId=-1):
    """Cremer-Pople amplitudes for ordered rings of at least four atoms."""
    coordinates, metadata = _geometry(mol, confId)
    indices = _ring(mol, ringAtomIndices)
    n = len(indices)
    if n < 4:
        raise ValueError('Cremer-Pople analysis requires at least four ring atoms')
    points = coordinates[indices] - coordinates[indices].mean(axis=0)
    angle = 2 * np.pi * np.arange(n) / n
    sine, cosine = np.sin(angle) @ points, np.cos(angle) @ points
    normal = np.cross(sine, cosine)
    norm = np.linalg.norm(normal)
    selection = dict(ring_atom_indices=indices.tolist())
    settings = dict(definition='Cremer-Pople 1975', normal='sine_cross_cosine', phase_zero_tolerance_a=1e-10)
    if norm <= 1e-12 * max(float(np.linalg.norm(sine) * np.linalg.norm(cosine)), 1e-24):
        values, units = dict(q_total_a=None), dict(q_total_a='angstrom')
        if n % 2 == 0:
            values['q_alternating_a'] = None
            units['q_alternating_a'] = 'angstrom'
        if n == 6:
            values['theta_deg'] = None
            units['theta_deg'] = 'degree'
        missing = dict.fromkeys(values, 'ring_normal_undefined')
        modes = []
        for mode in range(2, (n - 1) // 2 + 1):
            modes.append(dict(mode=mode, q_a=None, phi_deg=None))
            missing[f'mode_{mode}_q_a'] = 'ring_normal_undefined'
            missing[f'mode_{mode}_phi_deg'] = 'ring_normal_undefined'
        report = _report('ring_puckering', metadata, selection, settings, values, units, missing)
        report.update(modes=modes, mode_units=dict(q_a='angstrom', phi_deg='degree'))
        return report
    z = points @ (normal / norm)
    modes = []
    missing = {}
    for mode in range(2, (n - 1) // 2 + 1):
        c = float(np.sqrt(2 / n) * (z @ np.cos(mode * angle)))
        s = float(-np.sqrt(2 / n) * (z @ np.sin(mode * angle)))
        q = float(np.hypot(c, s))
        phase = float(np.degrees(np.arctan2(s, c)) % 360) if q > 1e-10 else None
        modes.append(dict(mode=mode, q_a=q, phi_deg=phase))
        if phase is None:
            missing[f'mode_{mode}_phi_deg'] = 'zero_mode_amplitude'
    alternating = float(z @ ((-1.0) ** np.arange(n)) / np.sqrt(n)) if n % 2 == 0 else None
    values = dict(q_total_a=float(np.linalg.norm(z)))
    units = dict(q_total_a='angstrom')
    if alternating is not None:
        values['q_alternating_a'] = alternating
        units['q_alternating_a'] = 'angstrom'
    if n == 6:
        values['theta_deg'] = (
            float(np.degrees(np.arctan2(modes[0]['q_a'], alternating)))
            if values['q_total_a'] > 1e-10
            else None
        )
        units['theta_deg'] = 'degree'
        if values['theta_deg'] is None:
            missing['theta_deg'] = 'planar_ring'
    report = _report('ring_puckering', metadata, selection, settings, values, units, missing)
    report.update(modes=modes, mode_units=dict(q_a='angstrom', phi_deg='degree'))
    return report


def CalcPyramidalization(mol, atomIdx, *, neighborIndices=None, confId=-1):
    """Haddon POAV1 pyramidalization angle, P=sin(theta_p), angle sum, and height.

    Analyzes a 3-coordinate center atom ``atomIdx`` using its three RDKit-bonded
    neighbors (or explicit ``neighborIndices`` of length 3). Computes:

    - ``pyramidalization_angle_deg`` (``theta_p = theta_sigma_pi - 90 deg``):
      ``0`` deg for trigonal planar, ``19.4712`` deg for ideal tetrahedral.
    - ``pyramidalization_p`` (``sin(theta_p)`` in ``[0, 1]``): ``0`` planar,
      ``1/3`` ideal tetrahedral.
    - ``poav_angle_deg`` (``theta_sigma_pi = 90 + theta_p``): ``90`` deg planar,
      ``109.4712`` deg ideal tetrahedral.
    - ``angle_sum_deg``: sum of the three bond angles around ``atomIdx``
      (``360`` deg planar, ``328.41`` deg ideal tetrahedral).
    - ``pyramid_height_a``: perpendicular distance of ``atomIdx`` from the plane
      of the three neighbor coordinates in Angstrom.

    Degenerate geometries (collinear unit bond tips or collinear neighbors)
    return ``None`` with explicit ``missing`` reasons rather than zero.
    """
    coordinates, metadata = _geometry(mol, confId)
    natoms = mol.GetNumAtoms()
    center = int(_indices([atomIdx], natoms, 'atomIdx')[0])
    if neighborIndices is None:
        neighbors = [int(n.GetIdx()) for n in mol.GetAtomWithIdx(center).GetNeighbors()]
        if len(neighbors) != 3:
            raise ValueError(
                f'CalcPyramidalization requires exactly 3 bonded neighbors for atom {center}, '
                f'found {len(neighbors)} (or pass neighborIndices explicitly)'
            )
    else:
        neighbors = _indices(neighborIndices, natoms, 'neighborIndices').tolist()
        if len(neighbors) != 3 or center in neighbors:
            raise ValueError('neighborIndices must contain 3 distinct atoms different from atomIdx')

    vecs = coordinates[neighbors] - coordinates[center]
    norms = np.linalg.norm(vecs, axis=1)
    if np.any(norms < 1e-10):
        raise ValueError('neighbor atoms must have distinct coordinates from atomIdx')
    unit = vecs / norms[:, None]

    # Sum of the three bond angles around the center atom.
    a01 = float(np.degrees(np.arccos(np.clip(float(np.dot(unit[0], unit[1])), -1.0, 1.0))))
    a12 = float(np.degrees(np.arccos(np.clip(float(np.dot(unit[1], unit[2])), -1.0, 1.0))))
    a20 = float(np.degrees(np.arccos(np.clip(float(np.dot(unit[2], unit[0])), -1.0, 1.0))))
    angle_sum = a01 + a12 + a20

    values = dict(
        pyramidalization_p=None,
        pyramidalization_angle_deg=None,
        poav_angle_deg=None,
        angle_sum_deg=float(angle_sum),
        pyramid_height_a=None,
    )
    units = dict(
        pyramidalization_p='dimensionless',
        pyramidalization_angle_deg='degree',
        poav_angle_deg='degree',
        angle_sum_deg='degree',
        pyramid_height_a='angstrom',
    )
    missing = {}

    # Haddon POAV1 normal to the plane of the three unit bond tips.
    tip_normal = np.cross(unit[1] - unit[0], unit[2] - unit[0])
    tip_norm = float(np.linalg.norm(tip_normal))
    if tip_norm <= 1e-12:
        for key in ('pyramidalization_p', 'pyramidalization_angle_deg', 'poav_angle_deg'):
            missing[key] = 'degenerate_unit_bond_geometry'
    else:
        n_hat = tip_normal / tip_norm
        p_sin = float(np.clip(abs(float(np.dot(n_hat, unit[0]))), 0.0, 1.0))
        theta_p = float(np.degrees(np.arcsin(p_sin)))
        values['pyramidalization_p'] = p_sin
        values['pyramidalization_angle_deg'] = theta_p
        values['poav_angle_deg'] = float(90.0 + theta_p)

    # Height of center above the plane through the three neighbor coordinates.
    p_normal = np.cross(coordinates[neighbors[1]] - coordinates[neighbors[0]],
                        coordinates[neighbors[2]] - coordinates[neighbors[0]])
    p_norm = float(np.linalg.norm(p_normal))
    if p_norm <= 1e-12:
        missing['pyramid_height_a'] = 'collinear_neighbors'
    else:
        values['pyramid_height_a'] = float(
            abs(float(np.dot(coordinates[center] - coordinates[neighbors[0]], p_normal / p_norm)))
        )

    selection = dict(atom_index=center, neighbor_indices=neighbors)
    settings = dict(
        definition='Haddon POAV1 pyramidalization (theta_p = theta_sigma_pi - 90 deg, P = sin(theta_p))',
    )
    return _report('pyramidalization', metadata, selection, settings, values, units, missing)


def _solve_exact_min_cone(dirs: np.ndarray, alphas: np.ndarray):
    """Exact minimum enclosing spherical cone of caps {(dirs[i], alphas[i])} on S^2.

    Implements the Bilbrey-Kazez-Locklin-Allen (2013) exact 1/2/3-cap solver:
    every minimal enclosing cone on the sphere is tangent to 1, 2, or 3 caps.
    """
    n = len(dirs)
    # Prune caps strictly enclosed by another cap to reduce combinatorics.
    dots = np.clip(dirs @ dirs.T, -1.0, 1.0)
    gammas = np.arccos(dots)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        for j in range(n):
            if i == j or not keep[j]:
                continue
            if gammas[i, j] + alphas[j] <= alphas[i] + 1e-12:
                if abs(alphas[i] - alphas[j]) <= 1e-12 and gammas[i, j] <= 1e-12 and j < i:
                    continue
                keep[j] = False
    active_idx = np.flatnonzero(keep)
    v = dirs[active_idx]
    a = alphas[active_idx]
    m = len(active_idx)
    tol = 1e-9

    def _enclosing_half_angle(axis: np.ndarray) -> float:
        return float(np.max(np.arccos(np.clip(v @ axis, -1.0, 1.0)) + a))

    best_theta = float('inf')
    best_axis = None
    best_support = ()

    # 1-cap candidates
    for i in range(m):
        cand_axis = v[i]
        theta = _enclosing_half_angle(cand_axis)
        if theta < best_theta - 1e-12:
            best_theta = theta
            best_axis = cand_axis
            best_support = (int(active_idx[i]),)

    # 2-cap candidates
    pair_theta = np.zeros((m, m), dtype=float)
    for i in range(m):
        for j in range(i + 1, m):
            gij = float(np.arccos(np.clip(float(np.dot(v[i], v[j])), -1.0, 1.0)))
            th = 0.5 * (gij + a[i] + a[j])
            pair_theta[i, j] = pair_theta[j, i] = th
            if th >= best_theta + tol:
                continue
            sin_g = np.sin(gij)
            if sin_g <= 1e-12:
                continue
            phi_i = th - a[i]
            if phi_i < -tol or phi_i > gij + tol:
                continue
            cand_axis = (np.sin(gij - phi_i) * v[i] + np.sin(phi_i) * v[j]) / sin_g
            cand_norm = float(np.linalg.norm(cand_axis))
            if cand_norm <= 1e-12:
                continue
            cand_axis /= cand_norm
            max_reach = _enclosing_half_angle(cand_axis)
            if max_reach <= th + 1e-7 and max_reach < best_theta - 1e-12:
                best_theta = max_reach
                best_axis = cand_axis
                best_support = (int(active_idx[i]), int(active_idx[j]))

    # 3-cap candidates (batched 3x3 linear solve for speed and exactness)
    if m >= 3:
        triples = []
        for i in range(m):
            for j in range(i + 1, m):
                if pair_theta[i, j] >= best_theta + 1e-7:
                    continue
                for k in range(j + 1, m):
                    if pair_theta[i, k] >= best_theta + 1e-7 or pair_theta[j, k] >= best_theta + 1e-7:
                        continue
                    triples.append((i, j, k))
        if triples:
            tri = np.asarray(triples, dtype=int)
            cos_a = np.cos(a)
            sin_a = np.sin(a)
            batch = 4096
            for start in range(0, len(tri), batch):
                sub = tri[start:start + batch]
                V = v[sub]  # (B, 3, 3)
                dets = np.linalg.det(V)
                good = np.abs(dets) > 1e-11
                if not np.any(good):
                    continue
                sub = sub[good]
                V = V[good]
                rhs = np.stack((cos_a[sub], sin_a[sub]), axis=-1)  # (B, 3, 2)
                sol = np.linalg.solve(V, rhs)  # (B, 3, 2)
                p = sol[:, :, 0]
                q = sol[:, :, 1]
                A = np.sum(p * p, axis=1) - 1.0
                B = 2.0 * np.sum(p * q, axis=1)
                C = np.sum(q * q, axis=1) - 1.0
                u = A - C
                vv = B
                w = -(A + C)
                R = np.hypot(u, vv)
                valid = (R > 1e-14) & (np.abs(w) <= R + 1e-10)
                if not np.any(valid):
                    continue
                sub = sub[valid]
                p = p[valid]
                q = q[valid]
                u = u[valid]
                vv = vv[valid]
                w = w[valid]
                R = R[valid]
                base_ang = np.arctan2(vv, u)
                delta_ang = np.arccos(np.clip(w / R, -1.0, 1.0))
                max_a_sub = np.max(a[sub], axis=1)
                for sign in (-1.0, 1.0):
                    phi = (base_ang + sign * delta_ang) % (2.0 * np.pi)
                    thetas = 0.5 * phi
                    cand_mask = (thetas >= max_a_sub - 1e-9) & (thetas < best_theta + 1e-7)
                    if not np.any(cand_mask):
                        continue
                    for idx in np.flatnonzero(cand_mask):
                        th = float(thetas[idx])
                        if th >= best_theta + 1e-7:
                            continue
                        cand_axis = np.cos(th) * p[idx] + np.sin(th) * q[idx]
                        cn = float(np.linalg.norm(cand_axis))
                        if abs(cn - 1.0) > 1e-5:
                            continue
                        cand_axis /= cn
                        max_reach = _enclosing_half_angle(cand_axis)
                        if max_reach <= th + 1e-6 and max_reach < best_theta - 1e-12:
                            best_theta = max_reach
                            best_axis = cand_axis
                            i0, j0, k0 = sub[idx]
                            best_support = (int(active_idx[i0]), int(active_idx[j0]), int(active_idx[k0]))

    # Recompute exact tangent atom set among all supplied caps within tolerance.
    reaches = np.arccos(np.clip(dirs @ best_axis, -1.0, 1.0)) + alphas
    final_theta = float(np.max(reaches))
    tangent = [int(i) for i in np.flatnonzero(reaches >= final_theta - 1e-6)]
    if not tangent:
        tangent = list(best_support)
    return final_theta, best_axis, tangent


def CalcConeAngle(
    mol,
    atomIndices,
    *,
    center,
    axisAtom=None,
    confId=-1,
    radiiAngstrom=None,
    radiusScale=1.0,
):
    """Exact minimum enclosing ligand cone angle (or Tolman axial cone angle).

    Given apex ``center`` (e.g. a metal atom) and ligand ``atomIndices``,
    computes the enclosing cone around the atomic van der Waals spheres:

    - When ``axisAtom is None``: computes the exact Bilbrey-Allen (2013)
      minimum enclosing cone over all unit axes on the sphere (strictly
      rigid-motion and atom-permutation invariant).
    - When ``axisAtom`` is provided: computes the Tolman-style axial cone angle
      around the fixed axis from ``center`` toward ``axisAtom``.

    If any selected atom sphere encloses ``center`` (``r_i >= d_i``), returns
    ``cone_angle_deg=None`` with explicit ``missing`` reason
    ``'apex_inside_atom_sphere'`` rather than returning a clipped angle.
    """
    coordinates, metadata = _geometry(mol, confId)
    natoms = mol.GetNumAtoms()
    selected = _indices(atomIndices, natoms, 'atomIndices')
    central = int(_indices([center], natoms, 'center')[0])
    if central in selected:
        raise ValueError('exclude the center atom from atomIndices')
    radii, settings = _radii(mol, radiiAngstrom, radiusScale)

    axis_idx = None
    fixed_axis = None
    if axisAtom is not None:
        axis_idx = int(_indices([axisAtom], natoms, 'axisAtom')[0])
        if axis_idx == central:
            raise ValueError('axisAtom must differ from center')
        axis_vec = coordinates[axis_idx] - coordinates[central]
        axis_norm = float(np.linalg.norm(axis_vec))
        if axis_norm < 1e-10:
            raise ValueError('axisAtom must have distinct coordinates from center')
        fixed_axis = axis_vec / axis_norm

    vecs = coordinates[selected] - coordinates[central]
    dists = np.linalg.norm(vecs, axis=1)
    if np.any(dists < 1e-10):
        raise ValueError('selected atoms must have distinct coordinates from center')

    sel_radii = radii[selected]
    overlap = sel_radii >= dists
    selection = dict(atom_indices=selected.tolist(), center=central)
    if axis_idx is not None:
        selection['axis_atom'] = axis_idx
    settings.update(
        mode='fixed_axis' if fixed_axis is not None else 'exact_minimum_enclosing',
        definition=(
            'Tolman axial cone angle around center->axisAtom'
            if fixed_axis is not None
            else 'Bilbrey-Allen (2013) exact minimum enclosing cone over spherical caps'
        ),
    )

    if np.any(overlap):
        overlapping_atoms = [int(selected[i]) for i in np.flatnonzero(overlap)]
        values = dict(cone_angle_deg=None, half_angle_deg=None)
        units = dict(cone_angle_deg='degree', half_angle_deg='degree')
        missing = dict(
            cone_angle_deg='apex_inside_atom_sphere',
            half_angle_deg='apex_inside_atom_sphere',
        )
        report = _report('cone_angle', metadata, selection, settings, values, units, missing)
        report.update(cone_axis=None, tangent_atom_indices=[], overlapping_atom_indices=overlapping_atoms)
        return report

    dirs = vecs / dists[:, None]
    alphas = np.arcsin(sel_radii / dists)
    if fixed_axis is not None:
        reaches = np.arccos(np.clip(dirs @ fixed_axis, -1.0, 1.0)) + alphas
        half_rad = float(np.max(reaches))
        cone_axis = fixed_axis
        tangent_local = [int(i) for i in np.flatnonzero(reaches >= half_rad - 1e-6)]
    else:
        half_rad, cone_axis, tangent_local = _solve_exact_min_cone(dirs, alphas)

    half_deg = min(180.0, float(np.degrees(half_rad)))
    cone_deg = min(360.0, 2.0 * half_deg)
    tangent_atoms = [int(selected[i]) for i in tangent_local]
    values = dict(cone_angle_deg=cone_deg, half_angle_deg=half_deg)
    units = dict(cone_angle_deg='degree', half_angle_deg='degree')
    report = _report('cone_angle', metadata, selection, settings, values, units)
    report.update(
        cone_axis=[float(x) for x in cone_axis],
        tangent_atom_indices=tangent_atoms,
        overlapping_atom_indices=[],
    )
    return report
