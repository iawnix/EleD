"""Descriptors from the retained xTB density and AO basis.

Every function in this module is an analysis of an existing result. It never
calls xTB. ``XTBCalculator(...).Calculate(..., keepOverlap=True)`` is required
because orbital coefficients without the same calculation's overlap matrix are
not a density in a defined non-orthogonal basis.
"""
from __future__ import annotations

from collections.abc import Mapping
from itertools import permutations
import numpy as np

from .electronic_descriptors import (
    _as_result,
    _finite,
    _indices,
    _nonnegative,
    _orbital_basis,
    _orbitals,
)
from .errors import DescriptorUnavailableError
from .mapping import _validate_mapping
from .populations import _frontier_group
from .selections import NormalizeSelections

_VALENCE = {1: 1, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7, 14: 4, 15: 5,
            16: 6, 17: 7, 32: 4, 33: 5, 34: 6, 35: 7, 53: 7}


def _basis(result):
    result = _as_result(result)
    if not all(len(v) for v in (result.aoShellMap, result.shellAtomMap, result.shellAngularMomenta)):
        raise DescriptorUnavailableError(
            "density descriptors require same-calculation overlap, coefficients and shell maps; "
            "use XTBCalculator(...).Calculate(..., keepOverlap=True)"
        )
    basis = _orbital_basis(result)
    atom = basis.ao_atom
    def integer_map(values, name):
        array = np.asarray(values)
        if array.ndim != 1 or array.dtype.kind not in 'iu':
            raise ValueError(f'{name} must contain integer indices')
        return array.astype(int)

    shell = integer_map(result.aoShellMap, 'AO shell map')
    shell_atom = integer_map(result.shellAtomMap, 'shell atom map')
    shell_l = integer_map(result.shellAngularMomenta, 'shell angular momentum map')
    if (len(shell) != len(atom) or len(shell_l) != len(shell_atom)
            or np.any(shell < 0) or np.any(shell >= len(shell_atom))
            or not np.array_equal(shell_atom[shell], atom)
            or not np.array_equal(np.bincount(shell, minlength=len(shell_l)), 2 * shell_l + 1)):
        raise ValueError('AO shell map is inconsistent with the atom and angular-momentum maps')
    if np.any(basis.occupations > 2 + 1e-8):
        raise ValueError('spin-summed occupations must not exceed two')
    return (result, basis.coefficients, basis.overlap, basis.density, basis.density_overlap,
            basis.lowdin_density, atom, shell, shell_atom, shell_l, basis.orthonormality_error)


def CalcDensityMatrix(result):
    """Return the spin-summed density ``P=C diag(n) C.T`` and its checks."""
    result, c, s, p, _, _, _, _, _, _, _ = _basis(result)
    electron_count = float(
        result.numElectrons if result.numElectrons is not None else np.sum(result.orbitalOccupancies)
    )
    trace = float(np.trace(p @ s))
    return {
        'density_matrix': p.tolist(),
        'tr_ps_e': trace,
        'n_e': electron_count,
        'n_err_e': trace - electron_count,
        'units': {
            'density_matrix': 'basis-dependent (trace with S gives e)',
            'tr_ps_e': 'e',
            'n_e': 'e',
            'n_err_e': 'e',
        },
        'invariance': 'P transforms with AO basis; trace(P S) is invariant; MO phase leaves P unchanged',
    }


def _rows(result, mull, low, natoms):
    out = []
    for a in range(natoms):
        row = {'atom_index': a}
        for label, values in (('mul', mull), ('low', low)):
            for l, name in ((0, 's'), (1, 'p'), (2, 'd')):
                row[f'{label}_{name}_e'] = float(values[a, l])
            p = values[a, 1]
            s = values[a, 0]
            row[f'{label}_ps'] = float(p / s) if abs(s) > 1e-14 else None
        if result.method == 'GFN0-xTB':
            row['mul_qrec_e'] = None
            row['low_qrec_e'] = None
        else:
            numbers = result.metadata.get('atomic_numbers', ())
            valence = _VALENCE.get(int(numbers[a]), None) if a < len(numbers) else None
            row['mul_qrec_e'] = float(valence - np.sum(mull[a])) if valence is not None else None
            row['low_qrec_e'] = float(valence - np.sum(low[a])) if valence is not None else None
        row['missing'] = {}
        for key, value in row.items():
            if value is None:
                row['missing'][key] = (
                    'zero_s_population'
                    if key.endswith('_ps')
                    else 'GFN0_density_is_not_EEQ_charge'
                    if result.method == 'GFN0-xTB'
                    else 'reference_valence_unavailable'
                )
        out.append(row)
    return out


def CalcShellPopulations(result):
    """Return atom/shell angular-momentum populations and electron checks."""
    result, _, s, p, ps, lowdin, atom, shell, shell_atom, shell_l, _ = _basis(result)
    natoms = len(result.atomicCharges)
    if not natoms:
        natoms = int(atom.max()) + 1
    if np.any(shell_l > 2):
        raise DescriptorUnavailableError('shell populations support s, p and d shells only')
    mull = np.zeros((natoms, 3))
    low = np.zeros((natoms, 3))
    shell_rows = []
    shell_mull_total = np.zeros(natoms)
    shell_low_total = np.zeros(natoms)
    for sh, (a, l) in enumerate(zip(shell_atom, shell_l)):
        ao = np.flatnonzero(shell == sh)
        m = float(np.trace(ps[np.ix_(ao, ao)]))
        q = float(np.trace(lowdin[np.ix_(ao, ao)]))
        shell_rows.append({
            'shell_index': int(sh),
            'atom_index': int(a),
            'angular_momentum': int(l),
            'mul_n_e': m,
            'low_n_e': q,
        })
        shell_mull_total[a] += m
        shell_low_total[a] += q
        mull[a, l] += m
        low[a, l] += q
    trace = float(np.trace(p @ s))
    electrons = float(
        result.numElectrons if result.numElectrons is not None else np.sum(result.orbitalOccupancies)
    )
    atom_rows = _rows(result, mull, low, natoms)
    return {
        'atoms': atom_rows,
        'shells': shell_rows,
        'tr_ps_e': trace,
        'n_e': electrons,
        'n_err_e': trace - electrons,
        'shell_err_e': {
            'mulliken': float(np.max(np.abs(shell_mull_total - np.asarray(mull).sum(axis=1)))),
            'lowdin': float(np.max(np.abs(shell_low_total - np.asarray(low).sum(axis=1)))),
        },
        'charge_convention': (
            'GFN0 shell populations are density populations; they are not '
            'the GFN0 EEQ partial charges.'
            if result.method == 'GFN0-xTB'
            else 'charge reconstruction is reference valence electrons minus density population'
        ),
        'units': {
            **{
                f'{label}_{field}_e': 'e'
                for label in ('mul', 'low')
                for field in ('s', 'p', 'd', 'n', 'qrec')
            },
            'mul_ps': 'dimensionless',
            'low_ps': 'dimensionless',
            'tr_ps_e': 'e',
            'n_e': 'e',
            'n_err_e': 'e',
            'shell_err_e': 'e',
        },
    }


def _bond_atom_pairs(mol):
    return [(b.GetIdx(), b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]


def CalcFrontierBondDescriptors(
    mol,
    result,
    *,
    degeneracyEV=1e-3,
    occupationTolerance=1e-6,
):
    """Mulliken overlap populations of the rotation-invariant HOMO/LUMO groups."""
    result = _as_result(result)
    _validate_mapping(mol, result)
    _nonnegative(degeneracyEV, "degeneracyEV")
    _nonnegative(occupationTolerance, "occupationTolerance")
    result, c, s, _, _, _, atom, _, _, _, _ = _basis(result)
    energies, occ, occupied, virtual = _orbitals(result)
    rows = []
    missing = {}
    for label, indices in (
        ('homo', _frontier_group(energies, occ, occupied, virtual, True, degeneracyEV, occupationTolerance)),
        ('lumo', _frontier_group(energies, occ, occupied, virtual, False, degeneracyEV, occupationTolerance)),
    ):
        if not indices:
            missing[label] = 'frontier_group_unavailable'
        pg = c[:, indices] @ c[:, indices].T / len(indices) if indices else None
        for bond, a, b in _bond_atom_pairs(mol):
            ia, ib = np.flatnonzero(atom == a), np.flatnonzero(atom == b)
            overlap = float(np.sum(pg[np.ix_(ia, ib)] * s[np.ix_(ia, ib)]) * 2.0) if pg is not None else None
            rows.append({
                'frontier_group': label,
                'orb_idx': indices,
                'bond_index': bond,
                'begin_atom_index': a,
                'end_atom_index': b,
                'mul_ovlp_e': overlap,
                'missing': {'mul_ovlp_e': missing[label]} if not indices else {},
            })
    return {
        'bonds': rows,
        'units': {'mul_ovlp_e': 'e'},
        'invariance': 'degenerate groups use a projector; orbital phase and group rotation invariant',
        'missing': missing,
    }


def _normal_for_atom(mol, atomIndex, confId=-1, planeTolAngstrom=0.1):
    """Fit the atom's ring planes (union for fused rings), or its neighbor plane."""
    positions = np.asarray(mol.GetConformer(confId).GetPositions())
    rings = [r for r in mol.GetRingInfo().AtomRings() if atomIndex in r]
    indices = (
        sorted({i for r in rings for i in r})
        if rings
        else [atomIndex, *[a.GetIdx() for a in mol.GetAtomWithIdx(atomIndex).GetNeighbors()]]
    )
    if len(indices) < 3:
        return None, 'insufficient_plane_neighbors'
    points = positions[indices] - positions[indices].mean(axis=0)
    _, singular, vh = np.linalg.svd(points, full_matrices=False)
    if singular[1] <= max(singular[0], 1.0) * 1e-8:
        return None, 'local_plane_degenerate'
    normal = vh[-1]
    if np.max(np.abs(points @ normal)) > planeTolAngstrom:
        return None, 'local_plane_nonplanar'
    return normal, None


def _pi_options(scheme, planeTolAngstrom, piRef, *, natoms=None):
    """Validate and snapshot pi settings, using JSON-stable reference keys."""
    if scheme not in ('mulliken', 'lowdin'):
        raise ValueError("pi scheme must be 'mulliken' or 'lowdin'")
    _nonnegative(planeTolAngstrom, 'planeTolAngstrom')
    if piRef is not None and not isinstance(piRef, Mapping):
        raise ValueError('piRef must map original atom indices to reference populations')
    references = {}
    for key, value in (piRef or {}).items():
        index = int(key) if isinstance(key, str) and key.isascii() and key.isdecimal() else key
        if (isinstance(index, bool) or not isinstance(index, (int, np.integer))
                or index < 0 or (isinstance(key, str) and str(index) != key)
                or (natoms is not None and index >= natoms)):
            raise ValueError('piRef keys must be original atom indices within the molecule')
        if str(index) in references:
            raise ValueError('piRef contains duplicate atom indices')
        _nonnegative(value, 'piRef population')
        if value > 2:
            raise ValueError('piRef population must be at most two electrons')
        references[str(index)] = float(value)
    return {'scheme': scheme, 'plane_tol_a': float(planeTolAngstrom), 'pi_ref': references}


def CalcPiDescriptors(mol, result, *, scheme='lowdin', planeTolAngstrom=0.1, piRef=None):
    """Local p projections; populations in e, closed-shell projected bond indices."""
    settings = _pi_options(scheme, planeTolAngstrom, piRef, natoms=mol.GetNumAtoms())
    result = _as_result(result)
    _validate_mapping(mol, result)
    result, _, _, _, ps, lowdin, atom, shell, _, shell_l, _ = _basis(result)
    references = {int(key): value for key, value in settings['pi_ref'].items()}
    directions = None
    if len(result.aoPDirections):
        directions = _finite(result.aoPDirections, 'AO p directions', 2)
        if directions.shape != (len(atom), 3):
            raise ValueError('AO p directions must have shape [AO, 3]')
        for sh, l in enumerate(shell_l):
            indices = np.flatnonzero(shell == sh)
            if l == 1:
                if len(indices) != 3 or not np.allclose(directions[indices] @ directions[indices].T, np.eye(3), atol=1e-8):
                    raise ValueError('p component directions must be orthonormal')
            elif np.any(directions[indices]):
                raise ValueError('non-p AO directions must be zero')
    conf_id = result.metadata['conformer_id']
    vectors, atoms, bonds = {}, [], []
    missing = {'atoms': {}, 'bonds': {}}
    population_metric = (ps + ps.T) / 2 if scheme == 'mulliken' else lowdin
    for a in mol.GetAtoms():
        i = a.GetIdx()
        row = {'atom_index': i, 'pi_n_e': None, 'pi_q_e': None, 'normal': None, 'missing': {}}
        reason = None
        if not (a.GetIsAromatic() or a.GetHybridization().name == 'SP2'):
            reason = 'atom_is_not_sp2_or_aromatic'
        normal = None
        if reason is None:
            normal, reason = _normal_for_atom(mol, i, conf_id, planeTolAngstrom)
        indices = np.flatnonzero((atom == i) & (shell_l[shell] == 1))
        if reason is None and len(indices) != 3:
            reason = 'requires_one_complete_p_shell'
        if reason is None and directions is None:
            reason = 'p_component_directions_unavailable'
        if reason is not None:
            row['missing'] = {'pi_n_e': reason, 'pi_q_e': reason}
        else:
            v = np.zeros(len(atom))
            v[indices] = directions[indices] @ normal
            vectors[i] = v
            population = float(v @ population_metric @ v)
            row.update(pi_n_e=population, normal=normal.tolist())
            reference = references.get(i)
            if reference is None and a.GetAtomicNum() == 6 and a.GetFormalCharge() == 0 and a.GetNumRadicalElectrons() == 0:
                reference = 1.0
            if reference is None:
                row['missing']['pi_q_e'] = 'pi_reference_unspecified'
            else:
                row['pi_q_e'] = float(reference - population)
        if row['missing']:
            missing['atoms'][str(i)] = dict(row['missing'])
        atoms.append(row)
    for index, a, b in _bond_atom_pairs(mol):
        reason = None
        if a not in vectors or b not in vectors:
            reason = 'endpoint_pi_projection_unavailable'
        elif result.numUnpairedElectrons is None:
            reason = 'spin_state_unknown'
        elif result.numUnpairedElectrons:
            reason = 'spin_density_unavailable_for_open_shell'
        value = None
        if reason is None:
            va, vb = vectors[a], vectors[b]
            value = float((va @ ps @ vb) * (vb @ ps @ va) if scheme == 'mulliken' else (va @ lowdin @ vb) ** 2)
        row = {
            'bond_index': index,
            'begin_atom_index': a,
            'end_atom_index': b,
            'pi_bo': value,
            'missing': {'pi_bo': reason} if reason else {},
        }
        bonds.append(row)
        if reason:
            missing['bonds'][str(index)] = dict(row['missing'])
    return {
        'atoms': atoms,
        'bonds': bonds,
        'missing': missing,
        **settings,
        'charge_reference': '1 e for neutral nonradical sp2 carbon; otherwise explicit piRef',
        'units': {'pi_n_e': 'e', 'pi_q_e': 'e', 'pi_bo': 'dimensionless'},
        'invariance': 'AO permutation/sign changes with component directions; rigid motion; local normal sign',
    }


def CalcMCI(mol, result, ringAtomIndices, *, scheme='mulliken'):
    """Calculate the closed-shell Giambiagi-Bultinck multicenter index (MCI) and I_ring.

    Uses AO-block matrix multiplication and trace on ``M = P @ S`` (Mulliken) or
    ``M = S^(1/2) @ P @ S^(1/2)`` (Löwdin), which is strictly invariant to rigid
    3D rotations/translations, AO permutations/sign flips, MO phase/degenerate
    occupied rotations, ring starting atom and ring direction.

    Open-shell or unknown-spin results return ``None`` with explicit ``missing``
    reasons instead of silently applying the closed-shell spin-summed identity.
    """
    if scheme not in {'mulliken', 'lowdin'}:
        raise ValueError("scheme must be 'mulliken' or 'lowdin'")
    result = _as_result(result)
    _validate_mapping(mol, result)
    result, _, _, _, ps, lowdin, atom, _, _, _, _ = _basis(result)
    natoms = result.metadata.get('atom_count', len(result.atomicCharges))
    indices = _indices(ringAtomIndices, natoms, 'ringAtomIndices')
    if len(indices) < 3:
        raise ValueError('ringAtomIndices must contain at least three unique atom indices')
    if len(indices) > 9:
        raise ValueError('ringAtomIndices supports rings of at most nine atoms')
    for first, second in zip(indices, np.roll(indices, -1)):
        if mol.GetBondBetweenAtoms(int(first), int(second)) is None:
            raise ValueError('ringAtomIndices must form an ordered bonded closed path')
    ring = tuple(map(int, indices))
    definition = (
        'Closed-shell Giambiagi-Bultinck multicenter index: '
        'I_ring = (1/4) * (Tr(M_{A1,A2}...M_{An,A1}) + Tr(M_{A1,An}...M_{A2,A1})), '
        'MCI = (1/4) * sum_{fixed-first perms} Tr(M_{A1,A2}...M_{An,A1}) '
        'with M = P*S (mulliken) or S^(1/2)*P*S^(1/2) (lowdin)'
    )
    unavailable = (
        'spin_state_unknown'
        if result.numUnpairedElectrons is None
        else 'spin_density_unavailable_for_open_shell'
        if result.numUnpairedElectrons
        else None
    )
    if unavailable is not None:
        return {
            'ring_atom_indices': list(ring),
            'mci': None,
            'i_ring': None,
            'scheme': scheme,
            'units': {'mci': 'dimensionless', 'i_ring': 'dimensionless'},
            'definition': definition,
            'missing': {'mci': unavailable, 'i_ring': unavailable},
        }
    metric = ps if scheme == 'mulliken' else lowdin
    ao_by_atom = {a: np.flatnonzero(atom == a) for a in ring}
    blocks = {
        (first, second): metric[np.ix_(ao_by_atom[first], ao_by_atom[second])]
        for first in ring
        for second in ring
        if first != second
    }

    def _cycle_trace(order):
        mat = blocks[(order[0], order[1])]
        for first, second in zip(order[1:-1], order[2:]):
            mat = mat @ blocks[(first, second)]
        return float(np.einsum('ij,ji->', mat, blocks[(order[-1], order[0])]))

    forward_trace = _cycle_trace(ring)
    reverse_order = (ring[0],) + ring[:0:-1]
    reverse_trace = forward_trace if scheme == 'lowdin' else _cycle_trace(reverse_order)
    i_ring = 0.25 * (forward_trace + reverse_trace)
    total_trace = sum(_cycle_trace((ring[0],) + tail) for tail in permutations(ring[1:]))
    mci = 0.25 * total_trace
    return {
        'ring_atom_indices': list(ring),
        'mci': float(mci),
        'i_ring': float(i_ring),
        'scheme': scheme,
        'units': {'mci': 'dimensionless', 'i_ring': 'dimensionless'},
        'definition': definition,
        'missing': {},
    }


def CalcMCIDescriptors(mol, result, rings, *, scheme='mulliken'):
    """Calculate MCI for a named mapping of ordered ring selections."""
    if hasattr(rings, 'rings'):
        rings = NormalizeSelections(mol, rings).rings
    else:
        rings = NormalizeSelections(mol, {'rings': rings}).rings
    if not rings:
        raise ValueError('rings must be a nonempty mapping of names to ordered atom indices')
    return {name: CalcMCI(mol, result, indices, scheme=scheme) for name, indices in rings.items()}


def CalcDensityDescriptors(
    mol,
    result,
    *,
    includeFrontierBonds=True,
    includePi=True,
    piScheme='lowdin',
    piPlaneTolAngstrom=0.1,
    piRef=None,
    mciRings=None,
    mciScheme='mulliken',
):
    """Aggregate all zero-cost descriptors available from one retained result."""
    pi_options = _pi_options(piScheme, piPlaneTolAngstrom, piRef, natoms=mol.GetNumAtoms())
    out = {
        'density': CalcDensityMatrix(result),
        'shell_populations': CalcShellPopulations(result),
    }
    if includeFrontierBonds:
        out['frontier_bonds'] = CalcFrontierBondDescriptors(mol, result)
    if includePi:
        out['pi'] = CalcPiDescriptors(
            mol,
            result,
            scheme=pi_options['scheme'],
            planeTolAngstrom=pi_options['plane_tol_a'],
            piRef=pi_options['pi_ref'],
        )
    if mciRings is not None:
        out['mci'] = CalcMCIDescriptors(mol, result, mciRings, scheme=mciScheme)
    return out
