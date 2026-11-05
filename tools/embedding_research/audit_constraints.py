"""Replay the pinned Embedder.cpp geometric predicates in Angstrom coordinates.

Signed volumes are scalar triple products (six times tetrahedron volume).
The center containment tests use unnormalized face normals, as upstream does.
"""
import itertools

import numpy as np


def center_faces(positions, atoms, tolerance):
    center, *neighbors = atoms
    if center == neighbors[-1]:
        return []
    faces = []
    for i in range(4):
        a, b, c, opposite = neighbors[i:] + neighbors[:i]
        normal = np.cross(positions[b]-positions[a], positions[c]-positions[a])
        d1 = float(normal @ (positions[opposite]-positions[a]))
        d2 = float(normal @ (positions[center]-positions[a]))
        faces.append({'face': [a, b, c], 'opposite': opposite,
                      'opposite_signed_volume': d1, 'center_signed_volume': d2,
                      'ok': abs(d1) >= tolerance and abs(d2) >= tolerance and (d1 < 0) == (d2 < 0)})
    return faces


def audit(positions, constraints):
    p = np.asarray(positions, dtype=float)[:, :3]
    if not np.isfinite(p).all():
        raise ValueError('nonfinite coordinates cannot be audited')
    chiral_atoms = set()
    for c in constraints['chiral_sets']:
        if c['atoms'][0] != c['atoms'][-1]:
            chiral_atoms.update(c['atoms'])
    pairs = []
    for i, j, lower, upper in constraints['bounds']:
        i, j = int(i), int(j)
        distance = float(np.linalg.norm(p[i]-p[j]))
        excess = max(lower-distance, distance-upper, 0.)
        if excess > 0:
            pairs.append({'atoms': [i, j], 'distance': distance, 'lower': lower, 'upper': upper,
                          'excess_angstrom': excess, 'in_final_chiral_check': i in chiral_atoms and j in chiral_atoms,
                          'exceeds_final_tolerance': excess > .1*upper})
    chiral = []
    for c in constraints['chiral_sets']:
        center, a, b, d, anchor = c['atoms']
        volume = float((p[a]-p[anchor]) @ np.cross(p[b]-p[anchor], p[d]-p[anchor]))
        lower, upper = c['lower'], c['upper']
        ok = not ((lower > 0 and volume < lower and (volume/lower < .8 or np.signbit(volume) != np.signbit(lower)))
                  or (upper < 0 and volume > upper and (volume/upper < .8 or np.signbit(volume) != np.signbit(upper))))
        faces = center_faces(p, c['atoms'], .1)
        chiral.append({'center': center, 'volume': volume, 'volume_ok': bool(ok),
                       'sign_ok': bool(volume > 0 if lower > 0 else volume < 0),
                       'center_inside': all(f['ok'] for f in faces), 'faces': faces})
    tetrahedral = []
    for c in constraints['tetrahedral_sets']:
        center, *neighbors = c['atoms']
        v = p[center]-p[neighbors]
        lengths = np.linalg.norm(v, axis=1)
        if (lengths == 0).any():
            raise ValueError('degenerate tetrahedral center')
        v /= lengths[:, None]
        volumes = [abs(float(v[i] @ np.cross(v[j], v[k]))) for i, j, k in itertools.combinations(range(4), 3)]
        threshold = .125 if c['flags'] & 1 else .5
        faces = center_faces(p, c['atoms'], .3)
        tetrahedral.append({'center': center, 'normalized_volumes': volumes, 'threshold': threshold,
                            'volume_ok': min(volumes) >= threshold,
                            'center_inside': all(f['ok'] for f in faces), 'faces': faces})
    failures = {
        'FINAL_CHIRAL_BOUNDS': any(v['in_final_chiral_check'] and v['exceeds_final_tolerance'] for v in pairs),
        'FINAL_CENTER_IN_VOLUME': any(not c['center_inside'] for c in chiral),
        'CHIRAL_VOLUME': any(not c['volume_ok'] for c in chiral),
        'CHECK_TETRAHEDRAL_CENTERS': any(not c['volume_ok'] or not c['center_inside'] for c in tetrahedral),
    }
    return {'distance_violations': pairs, 'chiral': chiral, 'tetrahedral': tetrahedral,
            'failed_predicates': [name for name, failed in failures.items() if failed]}
