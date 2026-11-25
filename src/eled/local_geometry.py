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
):
    """Percent of a center sphere occupied by the union of selected atom spheres."""
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
    step = 2 * radius / n
    grid = (np.arange(n) + 0.5) * step - radius
    y, z = np.meshgrid(grid, grid, indexing='ij')
    yz = np.column_stack((y.ravel(), z.ravel()))
    centers = coordinates[selected] - coordinates[central]
    active = np.linalg.norm(centers, axis=1) < radius + radii[selected]
    centers, atom_radii = centers[active], radii[selected][active]
    total, occupied = 0, 0
    for x in grid:
        plane = yz[np.sum(yz ** 2, axis=1) + x * x <= radius * radius]
        points = np.column_stack((np.full(len(plane), x), plane))
        inside = np.zeros(len(points), dtype=bool)
        for position, atom_radius in zip(centers, atom_radii):
            inside |= np.sum((points - position) ** 2, axis=1) <= atom_radius ** 2
        total += len(points)
        occupied += int(inside.sum())
    settings.update(
        sphere_radius_angstrom=radius,
        requested_spacing_angstrom=spacing,
        effective_spacing_angstrom=step,
        grid_points=total,
        occupied_grid_points=occupied,
        quadrature='centered_cartesian_midpoint',
    )
    return _report(
        'buried_volume',
        metadata,
        dict(atom_indices=selected.tolist(), center=central),
        settings,
        dict(vbur_percent=100.0 * occupied / total),
        dict(vbur_percent='percent'),
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
