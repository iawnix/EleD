"""Experimental diagnostics for libxtb 6.7.1 and tblite 0.8 s/p/d bases.

No engine is imported here. These explicit conventions are not a universal
AO conversion API. See docs/orbital-alignment.md for source references.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


# Target tblite components: s; (y,z,x); (xy,yz,z2,xz,x2-y2).
# libxtb d(z2) has the opposite phase. Cartesian cross terms already include
# sqrt(3) in libxtb's contracted functions; there is no extra scaling here.
_COMPONENTS = {
    0: ((0,), (1,), ("s",)),
    1: ((1, 2, 0), (1, 1, 1), ("py", "pz", "px")),
    2: ((2, 4, 1, 3, 0), (1, 1, -1, 1, 1),
        ("dxy", "dyz", "dz2", "dxz", "dx2-y2")),
}


def _indices(values, name):
    array = np.asarray(values)
    if array.ndim != 1 or array.dtype.kind not in "iu" or not len(array):
        raise ValueError(f"{name} must be a nonempty integer vector")
    return array.astype(np.int64)


def validate_tblite_basis(orbital_map, shell_map, angular_momenta, natoms, nao):
    """Require complete, contiguous atom/shell/AO maps, with zero-based indices."""
    orbital_map = _indices(orbital_map, "orbital_map")
    shell_map = _indices(shell_map, "shell_map")
    angular = _indices(angular_momenta, "angular_momenta")
    if len(shell_map) != len(angular) or np.any(angular < 0) or np.any(angular > 4):
        raise ValueError("invalid shell angular momenta")
    if (np.any(shell_map < 0) or np.any(shell_map >= natoms)
            or np.any(np.diff(shell_map) < 0)
            or not np.array_equal(np.unique(shell_map), np.arange(natoms))):
        raise ValueError("shell_map must cover atoms in input order")
    expected = np.repeat(np.arange(len(angular)), 2 * angular + 1)
    if len(orbital_map) != nao or not np.array_equal(orbital_map, expected):
        raise ValueError("orbital_map does not match contiguous spherical shells")
    return orbital_map, shell_map, angular


def libxtb_shells(numbers, parameter_file):
    """Read per-atom shells independently from the parameter file used by libxtb."""
    parameters = Path(parameter_file).read_text()
    elements = {}
    for match in re.finditer(r"(?m)^\$Z=\s*(\d+)[^\n]*\n(.*?)(?=^\$end)",
                             parameters, flags=re.S | re.M):
        ao = re.search(r"(?m)^\s*ao\s*=\s*(\S+)", match[2])
        if ao is None:
            raise ValueError(f"missing ao entry for element {match[1]}")
        tokens = re.findall(r"(\d+)([spdfg])", ao[1])
        if not tokens or "".join(n + l for n, l in tokens) != ao[1]:
            raise ValueError(f"unsupported ao entry: {ao[1]}")
        elements[int(match[1])] = [(int(n), "spdfg".index(l)) for n, l in tokens]
    shells = []
    for atom, number in enumerate(numbers):
        if int(number) not in elements:
            raise ValueError(f"missing parameter shells for atomic number {number}")
        for local, (n, l) in enumerate(elements[int(number)]):
            shells.append({"atom": atom, "shell_on_atom": local, "n": n, "l": l})
    return shells


def libxtb_to_tblite_rows(numbers, parameter_file, orbital_map, shell_map,
                          angular_momenta, nao):
    """Return source row indices, row phases and target AO labels.

    Apply as ``C_tblite_basis = C_libxtb[rows, :] * phases[:, None]``.
    Shell topology must agree with independently read libxtb parameters.
    Unsupported angular momenta fail explicitly; no numerical row fitting.
    """
    _, atom_map, angular = validate_tblite_basis(
        orbital_map, shell_map, angular_momenta, len(numbers), nao)
    shells = libxtb_shells(numbers, parameter_file)
    if (len(shells) != len(angular)
            or [(s["atom"], s["l"]) for s in shells] != list(zip(atom_map, angular))):
        raise ValueError("libxtb parameter shells differ from tblite basis maps")
    rows, phases, labels = [], [], []
    offset = 0
    for shell_index, shell in enumerate(shells):
        l = shell["l"]
        if l not in _COMPONENTS:
            raise ValueError(f"AO alignment is only verified for s/p/d, got l={l}")
        permutation, signs, components = _COMPONENTS[l]
        for row, phase, component in zip(permutation, signs, components):
            rows.append(offset + row)
            phases.append(phase)
            labels.append({**shell, "shell": shell_index, "component": component})
        offset += 2 * l + 1
    if offset != nao:
        raise ValueError("libxtb and tblite AO counts differ")
    return np.asarray(rows), np.asarray(phases), labels


def _array(values, shape, name):
    array = np.asarray(values, dtype=float)
    if array.shape != shape or not np.isfinite(array).all():
        raise ValueError(f"{name} must be finite with shape {shape}")
    return array


def compare_orbitals(reference, candidate, overlap, reference_energies,
                     candidate_energies, reference_occupations, candidate_occupations,
                     *, degeneracy_ev=1e-3, occupation_tolerance=1e-6):
    """Compare real MO matrices in the SAME AO basis using the reference overlap.

    Sorted orbitals are grouped only if the full energy span in BOTH engines
    is <= degeneracy_ev, and occupations are compatible. Report sign-only
    errors separately from S-metric Procrustes errors within those groups.
    Group rotations are diagnostics, not a redefinition of eigenvectors.
    """
    reference = np.asarray(reference, dtype=float)
    if reference.ndim != 2 or reference.shape[0] != reference.shape[1] or not reference.size:
        raise ValueError("a full square AO/MO matrix is required")
    n = reference.shape[0]
    a = _array(reference, (n, n), "reference coefficients")
    b = _array(candidate, (n, n), "candidate coefficients")
    s = _array(overlap, (n, n), "overlap")
    ea = _array(reference_energies, (n,), "reference energies")
    eb = _array(candidate_energies, (n,), "candidate energies")
    oa = _array(reference_occupations, (n,), "reference occupations")
    ob = _array(candidate_occupations, (n,), "candidate occupations")
    if (not np.isfinite(degeneracy_ev) or degeneracy_ev < 0
            or not np.isfinite(occupation_tolerance) or occupation_tolerance < 0):
        raise ValueError("grouping tolerances must be finite and nonnegative")
    if np.any(np.diff(ea) < -1e-10) or np.any(np.diff(eb) < -1e-10):
        raise ValueError("orbital energies must be sorted")
    if not np.allclose(s, s.T, atol=1e-12, rtol=0):
        raise ValueError("overlap must be symmetric")
    # W.T @ W = S, so W @ C represents the same orbitals in an orthonormal basis.
    w = np.linalg.cholesky(s).T
    wa, wb = w @ a, w @ b
    if np.linalg.matrix_rank(wa) != n or np.linalg.matrix_rank(wb) != n:
        raise ValueError("coefficient matrix is rank deficient")
    column_signs = np.where(np.sum(wa * wb, axis=0) < 0, -1.0, 1.0)
    aligned = b.copy()
    groups = []
    start = 0
    while start < n:
        end = start + 1
        while end < n:
            occ = np.r_[oa[start:end + 1], ob[start:end + 1]]
            if (ea[end] - ea[start] > degeneracy_ev
                    or eb[end] - eb[start] > degeneracy_ev
                    or np.ptp(occ) > occupation_tolerance):
                break
            end += 1
        x, y = wa[:, start:end], wb[:, start:end]
        # QR normalization makes principal angles independent of norm defects,
        # which are reported separately and must also be checked by callers.
        qx, _ = np.linalg.qr(x)
        qy, _ = np.linalg.qr(y)
        singular = np.linalg.svd(qx.T @ qy, compute_uv=False)
        u, _, vh = np.linalg.svd(x.T @ y)
        aligned[:, start:end] = b[:, start:end] @ vh.T @ u.T
        # Stable projector distance, including angles too small for acos.
        projector_error = float(np.linalg.norm(qy - qx @ (qx.T @ qy), ord=2))
        groups.append({"start": start, "stop": end, "size": end - start,
                       "reference_span_ev": float(ea[end - 1] - ea[start]),
                       "candidate_span_ev": float(eb[end - 1] - eb[start]),
                       "minimum_principal_cosine": float(singular.min()),
                       "projector_distance": projector_error})
        start = end
    density_a, density_b = (wa * oa) @ wa.T, (wb * ob) @ wb.T
    return {
        "reference_orthonormality_error": float(np.max(abs(wa.T @ wa - np.eye(n)))),
        "candidate_orthonormality_error": float(np.max(abs(wb.T @ wb - np.eye(n)))),
        "sign_aligned_max_abs": float(np.max(abs(a - b * column_signs))),
        "subspace_aligned_max_abs": float(np.max(abs(a - aligned))),
        "max_projector_distance": max(g["projector_distance"] for g in groups),
        "density_orthonormal_max_abs": float(np.max(abs(density_a - density_b))),
        "density_trace_difference": float(abs(np.trace(density_a) - np.trace(density_b))),
        "degeneracy_ev": degeneracy_ev,
        "occupation_tolerance": occupation_tolerance,
        "groups": groups,
    }
