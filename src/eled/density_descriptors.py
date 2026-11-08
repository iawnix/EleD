"""Descriptors from the retained xTB density and AO basis.

Every function in this module is an analysis of an existing result.  It never
calls xTB.  ``XTBCalculator(...).calculate(..., keep_overlap=True)`` is required because orbital
coefficients without the same calculation's overlap matrix are not a density
in a defined non-orthogonal basis.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .electronic_descriptors import _as_result, _finite, _orbitals
from .errors import DescriptorUnavailableError

_VALENCE = {1: 1, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 14: 4, 15: 5,
            16: 6, 17: 7, 32: 4, 33: 5, 34: 6, 35: 7, 53: 7}


def _basis(result):
    result = _as_result(result)
    fields = (result.overlap_matrix, result.orbital_coefficients,
              result.ao_atom_indices, result.ao_shell_indices,
              result.shell_atom_indices, result.shell_angular_momenta)
    if not all(len(v) for v in fields):
        raise DescriptorUnavailableError(
            "density descriptors require same-calculation overlap, coefficients and shell maps; "
            "use XTBCalculator(...).calculate(..., keep_overlap=True)")
    c = _finite(result.orbital_coefficients, 'orbital coefficients', 2)
    s = _finite(result.overlap_matrix, 'overlap', 2)
    if c.shape != s.shape or s.shape[0] != len(result.ao_atom_indices):
        raise ValueError('density basis arrays have inconsistent dimensions')
    if not np.allclose(s, s.T, atol=1e-10, rtol=0):
        raise ValueError('overlap must be symmetric')
    eig, u = np.linalg.eigh(s)
    if eig.min() <= max(1., eig.max()) * 1e-10:
        raise ValueError('overlap must be positive definite')
    energies, occ, _, _ = _orbitals(result)
    if c.shape != (len(energies), len(energies)):
        raise ValueError('coefficient matrix does not match orbital arrays')
    metric = c.T @ s @ c
    err = float(np.max(np.abs(metric - np.eye(len(energies)))))
    if err > 1e-5:
        raise ValueError(f'coefficients are not S-orthonormal: error={err:g}')
    atom = np.asarray(result.ao_atom_indices, dtype=int)
    shell = np.asarray(result.ao_shell_indices, dtype=int)
    shell_atom = np.asarray(result.shell_atom_indices, dtype=int)
    shell_l = np.asarray(result.shell_angular_momenta, dtype=int)
    if len(shell) != len(atom) or np.any(shell < 0) or np.any(shell >= len(shell_atom)):
        raise ValueError('AO shell map is invalid')
    if not np.array_equal(shell_atom[shell], atom):
        raise ValueError('AO shell and atom maps disagree')
    if len(shell_l) != len(shell_atom):
        raise ValueError('shell angular momentum map is invalid')
    p = (c * occ[None, :]) @ c.T
    ps = p @ s
    lowdin_metric = (u * np.sqrt(eig)) @ u.T @ p @ (u * np.sqrt(eig)) @ u.T
    return result, c, s, p, ps, lowdin_metric, atom, shell, shell_atom, shell_l, err


def CalcDensityMatrix(result):
    """Return the spin-summed density ``P=C diag(n) C.T`` and its checks."""
    result, c, s, p, _, _, _, _, _, _, _ = _basis(result)
    electron_count = float(np.sum(result.orbital_occupations))
    trace = float(np.trace(p @ s))
    return {'density_matrix': p.tolist(), 'trace_PS_e': trace,
            'electron_count_e': electron_count,
            'electron_conservation_error_e': trace - electron_count,
            'units': {'density_matrix': 'basis-dependent (trace with S gives e)', 'trace_PS_e': 'e'},
            'invariance': 'invariant to AO and orbital phase changes when maps and S are transformed together'}


def _rows(result, mull, low, shell_atom, shell_l, natoms):
    out = []
    for a in range(natoms):
        row = {'atom_index': a}
        for label, values in (('mulliken', mull), ('lowdin', low)):
            for l, name in ((0, 's'), (1, 'p'), (2, 'd')):
                row[f'{label}_{name}_electron_count_e'] = float(values[a, l])
            p = values[a, 1]; s = values[a, 0]
            row[f'{label}_p_to_s_ratio'] = float(p / s) if abs(s) > 1e-14 else None
        if result.method == 'GFN0-xTB':
            row['mulliken_charge_reconstruction_e'] = None
            row['lowdin_charge_reconstruction_e'] = None
        else:
            numbers = result.metadata.get('atomic_numbers', ())
            valence = _VALENCE.get(int(numbers[a]), None) if a < len(numbers) else None
            row['mulliken_charge_reconstruction_e'] = float(valence - np.sum(mull[a])) if valence is not None else None
            row['lowdin_charge_reconstruction_e'] = float(valence - np.sum(low[a])) if valence is not None else None
        # Backward-friendly alias: the Löwdin reconstruction is the default
        # positive-definite population convention in this package.
        row['charge_reconstruction_e'] = row['lowdin_charge_reconstruction_e']
        out.append(row)
    return out


def CalcShellPopulations(result):
    """Return atom/shell angular-momentum populations and electron checks."""
    result, _, s, p, ps, lowdin, atom, shell, shell_atom, shell_l, _ = _basis(result)
    natoms = len(result.partial_charges)
    if not natoms:
        natoms = int(atom.max()) + 1
    if np.any(shell_l < 0) or np.any(shell_l > 2):
        # Keep unsupported shells visible without assigning them to d.
        supported = {0, 1, 2}
    else:
        supported = {0, 1, 2}
    mull = np.zeros((natoms, 3)); low = np.zeros((natoms, 3))
    shell_rows = []
    shell_mull_total = np.zeros(natoms); shell_low_total = np.zeros(natoms)
    for sh, (a, l) in enumerate(zip(shell_atom, shell_l)):
        ao = np.flatnonzero(shell == sh)
        m = float(np.trace(ps[np.ix_(ao, ao)]))
        q = float(np.trace(lowdin[np.ix_(ao, ao)]))
        shell_rows.append({'shell_index': int(sh), 'atom_index': int(a),
                           'angular_momentum': int(l), 'mulliken_electron_count_e': m,
                           'lowdin_electron_count_e': q})
        shell_mull_total[a] += m; shell_low_total[a] += q
        if l in supported and l < 3:
            mull[a, l] += m; low[a, l] += q
    trace = float(np.trace(p @ s)); electrons = float(np.sum(result.orbital_occupations))
    atom_rows = _rows(result, mull, low, shell_atom, shell_l, natoms)
    return {'atoms': atom_rows,
            'shells': shell_rows, 'trace_PS_e': trace,
            'electron_count_e': electrons,
            'electron_conservation_error_e': trace - electrons,
            'shell_atom_conservation_error_e': {
                'mulliken': float(np.max(np.abs(shell_mull_total - np.asarray(mull).sum(axis=1)))),
                'lowdin': float(np.max(np.abs(shell_low_total - np.asarray(low).sum(axis=1)))),
            },
            'charge_convention': ('GFN0 shell populations are density populations; they are not '
                                  'the GFN0 EEQ partial charges.' if result.method == 'GFN0-xTB' else
                                  'charge reconstruction is reference valence electrons minus density population'),
            'units': {'electron_count': 'e', 'p_to_s_ratio': 'dimensionless'}}


def _bond_atom_pairs(mol):
    return [(b.GetIdx(), b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]


def CalcMayerBondDescriptors(mol, result):
    """Return Mayer and Wiberg orders for RDKit bonds.

    xTB exposes a spin-summed density here.  Mayer order is therefore withheld
    for open-shell results instead of silently presenting a closed-shell proxy.
    """
    result, _, s, p, ps, _, atom, _, _, _, _ = _basis(result)
    wbo = np.asarray(result.wiberg_bond_orders, dtype=float)
    if wbo.shape != (mol.GetNumAtoms(), mol.GetNumAtoms()):
        raise ValueError('Wiberg matrix does not match molecule')
    open_shell = bool(result.number_of_unpaired_electrons)
    rows = []
    for index, a, b in _bond_atom_pairs(mol):
        ia, ib = np.flatnonzero(atom == a), np.flatnonzero(atom == b)
        value = None if open_shell else float(np.sum(ps[np.ix_(ia, ib)] * ps[np.ix_(ib, ia)].T))
        rows.append({'bond_index': index, 'begin_atom_index': a, 'end_atom_index': b,
                     'mayer_bond_order': value, 'wiberg_bond_order': float(wbo[a, b]),
                     'mayer_missing_reason': 'spin_density_unavailable_for_open_shell' if open_shell else None})
    return {'bonds': rows, 'units': {'mayer_bond_order': 'dimensionless', 'wiberg_bond_order': 'dimensionless'},
            'invariance': 'invariant to AO permutation and orbital phase'}


def _frontier_group(energies, occ, occupied, virtual, last, degeneracy_ev, occupation_tolerance):
    indices = occupied if last else virtual
    if not len(indices): return []
    anchor = int(indices[-1] if last else indices[0]); selected = [anchor]
    candidates = indices[-2::-1] if last else indices[1:]
    for i in candidates:
        if (abs(energies[i] - energies[anchor]) > degeneracy_ev or
                abs(occ[i] - occ[anchor]) > occupation_tolerance or abs(int(i)-selected[-1]) != 1): break
        selected.append(int(i))
    return sorted(selected)


def CalcFrontierBondDescriptors(mol, result, *, degeneracy_ev=1e-3,
                                   occupation_threshold=1e-7, occupation_tolerance=1e-6):
    """Mulliken overlap populations of the rotation-invariant HOMO/LUMO groups."""
    result, c, s, _, _, _, atom, _, _, _, _ = _basis(result)
    energies, occ, occupied, virtual = _orbitals(result, occupation_threshold)
    rows = []
    for label, indices in (('homo', _frontier_group(energies, occ, occupied, virtual, True, degeneracy_ev, occupation_tolerance)),
                           ('lumo', _frontier_group(energies, occ, occupied, virtual, False, degeneracy_ev, occupation_tolerance))):
        if not indices: continue
        # The projector, rather than individual MO coefficients, is invariant
        # under phase changes and rotations in a degenerate group.
        pg = c[:, indices] @ c[:, indices].T / len(indices)
        for bond, a, b in _bond_atom_pairs(mol):
            ia, ib = np.flatnonzero(atom == a), np.flatnonzero(atom == b)
            overlap = float(np.sum(pg[np.ix_(ia, ib)] * s[np.ix_(ia, ib)]) * 2.0)
            rows.append({'frontier_group': label, 'orbital_indices': indices,
                         'bond_index': bond, 'begin_atom_index': a, 'end_atom_index': b,
                         'mulliken_overlap_population': overlap})
    return {'bonds': rows, 'units': {'mulliken_overlap_population': 'e'},
            'invariance': 'degenerate groups use a projector; orbital phase and group rotation invariant',
            'missing': {} if rows else {'frontier_bonds': 'no occupied or virtual frontier group'}}


def _normal_for_atom(mol, atom_index, conf_id=-1):
    atom = mol.GetAtomWithIdx(atom_index); pos = np.asarray(mol.GetConformer(conf_id).GetPositions())
    ring = next((list(r) for r in mol.GetRingInfo().AtomRings() if atom_index in r and len(r) >= 3), None)
    if ring is not None:
        points = pos[ring] - pos[ring].mean(axis=0)
        _, _, vh = np.linalg.svd(points, full_matrices=False)
        n = vh[-1]
    else:
        vectors = [pos[n.GetIdx()] - pos[atom_index] for n in atom.GetNeighbors()]
        candidates = [np.cross(vectors[i], vectors[j]) for i in range(len(vectors)) for j in range(i)]
        n = max(candidates, key=lambda x: np.linalg.norm(x), default=None)
    if n is None or np.linalg.norm(n) < 1e-8: return None
    return n / np.linalg.norm(n)


def CalcPiDescriptors(mol, result, *, scheme='lowdin'):
    """Project p shells onto local normals for sp2/aromatic atoms and bonds."""
    if scheme not in ('mulliken', 'lowdin'): raise ValueError("scheme must be 'mulliken' or 'lowdin'")
    from .mapping import _validate_mapping
    probe = _as_result(result)
    if probe.metadata.get('ordered_structure_sha256'):
        _validate_mapping(mol, probe)
    elif len(probe.metadata.get('atomic_numbers', ())) != mol.GetNumAtoms():
        raise ValueError('molecule and result atom counts do not match')
    result, _, s, p, ps, lowdin, atom, shell, shell_atom, shell_l, _ = _basis(result)
    conf_id = result.metadata.get('conformer_id', -1)
    normals = {}; missing = {}
    for a in range(mol.GetNumAtoms()):
        if not (mol.GetAtomWithIdx(a).GetIsAromatic() or mol.GetAtomWithIdx(a).GetHybridization().name == 'SP2'):
            missing[str(a)] = 'atom_is_not_sp2_or_aromatic'; continue
        n = _normal_for_atom(mol, a, conf_id)
        if n is None: missing[str(a)] = 'local_plane_normal_undefined'
        else: normals[a] = n
    root = (lambda: None)
    values = []
    for a, n in normals.items():
        ao = np.flatnonzero((atom == a) & np.isin(shell_l[shell], [1]))
        if len(ao) != 3: missing[str(a)] = 'p_shell_mapping_not_three_components'; continue
        v = np.zeros(len(atom)); v[ao] = n
        d = p if scheme == 'mulliken' else lowdin
        electron = float(v @ d @ v)
        valence = _VALENCE.get(int(result.metadata.get('atomic_numbers', [])[a]))
        values.append({'atom_index': a, 'pi_electron_count_e': electron,
                       'pi_charge_e': float(valence - electron) if valence is not None else None,
                       'normal': n.tolist()})
    bonds = []
    for bond, a, b in _bond_atom_pairs(mol):
        if a not in normals or b not in normals: continue
        pa = np.flatnonzero((atom == a) & np.isin(shell_l[shell], [1])); pb = np.flatnonzero((atom == b) & np.isin(shell_l[shell], [1]))
        if len(pa) != 3 or len(pb) != 3: continue
        va = np.zeros(len(atom)); vb = np.zeros(len(atom)); va[pa] = normals[a]; vb[pb] = normals[b]
        # Reversing either local normal changes the sign of the projected
        # orbital, but cannot change a bond order.  Report its magnitude.
        bond_order = float(abs(2.0 * (va @ p @ vb) * (va @ s @ vb)))
        bonds.append({'bond_index': bond, 'begin_atom_index': a, 'end_atom_index': b, 'pi_bond_order': bond_order})
    return {'atoms': values, 'bonds': bonds, 'missing': missing,
            'units': {'pi_electron_count_e': 'e', 'pi_charge_e': 'e', 'pi_bond_order': 'dimensionless'},
            'invariance': 'rigid-motion invariant; normal sign does not change atom populations or bond order'}


def CalcDensityDescriptors(mol, result, *, include_frontier_bonds=True, include_pi=True):
    """Aggregate all zero-cost descriptors available from one retained result."""
    out = {'density': CalcDensityMatrix(result), 'shell_populations': CalcShellPopulations(result),
           'mayer_bonds': CalcMayerBondDescriptors(mol, result)}
    if include_frontier_bonds: out['frontier_bonds'] = CalcFrontierBondDescriptors(mol, result)
    if include_pi: out['pi'] = CalcPiDescriptors(mol, result)
    return out

