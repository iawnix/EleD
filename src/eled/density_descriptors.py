"""Descriptors from the retained xTB density and AO basis.

Every function in this module is an analysis of an existing result.  It never
calls xTB.  ``XTBCalculator(...).calculate(..., keep_overlap=True)`` is required because orbital
coefficients without the same calculation's overlap matrix are not a density
in a defined non-orthogonal basis.
"""
from __future__ import annotations

from collections.abc import Mapping
import numpy as np

from .electronic_descriptors import _as_result, _finite, _orbitals, _nonnegative
from .errors import DescriptorUnavailableError
from .mapping import _validate_mapping
from .populations import _frontier_group

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
    natoms = result.metadata.get('atom_count', len(result.partial_charges))
    def integer_map(values, name):
        array = np.asarray(values)
        if array.ndim != 1 or array.dtype.kind not in 'iu':
            raise ValueError(f'{name} must contain integer indices')
        return array.astype(int)
    atom = integer_map(result.ao_atom_indices, 'AO atom map')
    shell = integer_map(result.ao_shell_indices, 'AO shell map')
    shell_atom = integer_map(result.shell_atom_indices, 'shell atom map')
    shell_l = integer_map(result.shell_angular_momenta, 'shell angular momentum map')
    if (len(shell) != len(atom) or np.any(shell < 0) or np.any(shell >= len(shell_atom))):
        raise ValueError('AO shell map is invalid')
    if (np.any(atom < 0) or np.any(atom >= natoms) or np.any(shell_atom < 0)
            or np.any(shell_atom >= natoms) or not np.array_equal(np.unique(atom),np.arange(natoms))):
        raise ValueError('AO atom map must cover the calculation atoms')
    if not np.array_equal(shell_atom[shell], atom):
        raise ValueError('AO shell and atom maps disagree')
    if len(shell_l) != len(shell_atom) or np.any(shell_l < 0):
        raise ValueError('shell angular momentum map is invalid')
    if not np.array_equal(np.bincount(shell, minlength=len(shell_l)), 2*shell_l+1):
        raise ValueError('shell size must equal 2*l+1')
    if np.any(occ > 2 + 1e-8):
        raise ValueError('spin-summed occupations must not exceed two')
    p = (c * occ[None, :]) @ c.T
    ps = p @ s
    lowdin_metric = (u * np.sqrt(eig)) @ u.T @ p @ (u * np.sqrt(eig)) @ u.T
    return result, c, s, p, ps, lowdin_metric, atom, shell, shell_atom, shell_l, err


def CalcDensityMatrix(result):
    """Return the spin-summed density ``P=C diag(n) C.T`` and its checks."""
    result, c, s, p, _, _, _, _, _, _, _ = _basis(result)
    electron_count = float(result.number_of_electrons if result.number_of_electrons is not None else np.sum(result.orbital_occupations))
    trace = float(np.trace(p @ s))
    return {'density_matrix': p.tolist(), 'tr_ps_e': trace,
            'n_e': electron_count,
            'n_err_e': trace - electron_count,
            'units': {'density_matrix': 'basis-dependent (trace with S gives e)', 'tr_ps_e': 'e', 'n_e': 'e', 'n_err_e': 'e'},
            'invariance': 'P transforms with AO basis; trace(P S) is invariant; MO phase leaves P unchanged'}


def _rows(result, mull, low, natoms):
    out = []
    for a in range(natoms):
        row = {'atom_index': a}
        for label, values in (('mul', mull), ('low', low)):
            for l, name in ((0, 's'), (1, 'p'), (2, 'd')):
                row[f'{label}_{name}_e'] = float(values[a, l])
            p = values[a, 1]; s = values[a, 0]
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
                row['missing'][key] = ('zero_s_population' if key.endswith('_ps') else
                    'GFN0_density_is_not_EEQ_charge' if result.method == 'GFN0-xTB' else
                    'reference_valence_unavailable')
        out.append(row)
    return out


def CalcShellPopulations(result):
    """Return atom/shell angular-momentum populations and electron checks."""
    result, _, s, p, ps, lowdin, atom, shell, shell_atom, shell_l, _ = _basis(result)
    natoms = len(result.partial_charges)
    if not natoms:
        natoms = int(atom.max()) + 1
    if np.any(shell_l > 2):
        raise DescriptorUnavailableError('shell populations support s, p and d shells only')
    mull = np.zeros((natoms, 3)); low = np.zeros((natoms, 3))
    shell_rows = []
    shell_mull_total = np.zeros(natoms); shell_low_total = np.zeros(natoms)
    for sh, (a, l) in enumerate(zip(shell_atom, shell_l)):
        ao = np.flatnonzero(shell == sh)
        m = float(np.trace(ps[np.ix_(ao, ao)]))
        q = float(np.trace(lowdin[np.ix_(ao, ao)]))
        shell_rows.append({'shell_index': int(sh), 'atom_index': int(a),
                           'angular_momentum': int(l), 'mul_n_e': m,
                           'low_n_e': q})
        shell_mull_total[a] += m; shell_low_total[a] += q
        mull[a, l] += m; low[a, l] += q
    trace = float(np.trace(p @ s)); electrons = float(result.number_of_electrons if result.number_of_electrons is not None else np.sum(result.orbital_occupations))
    atom_rows = _rows(result, mull, low, natoms)
    return {'atoms': atom_rows,
            'shells': shell_rows, 'tr_ps_e': trace,
            'n_e': electrons,
            'n_err_e': trace - electrons,
            'shell_err_e': {
                'mulliken': float(np.max(np.abs(shell_mull_total - np.asarray(mull).sum(axis=1)))),
                'lowdin': float(np.max(np.abs(shell_low_total - np.asarray(low).sum(axis=1)))),
            },
            'charge_convention': ('GFN0 shell populations are density populations; they are not '
                                  'the GFN0 EEQ partial charges.' if result.method == 'GFN0-xTB' else
                                  'charge reconstruction is reference valence electrons minus density population'),
            'units': {**{f'{label}_{field}_e': 'e' for label in ('mul','low') for field in ('s','p','d','n','qrec')},
                      'mul_ps': 'dimensionless', 'low_ps': 'dimensionless',
                      'tr_ps_e': 'e', 'n_e': 'e', 'n_err_e': 'e', 'shell_err_e': 'e'}}


def _bond_atom_pairs(mol):
    return [(b.GetIdx(), b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]


def CalcMayerBondDescriptors(mol, result):
    """Return Mayer and Wiberg orders for RDKit bonds.

    xTB exposes a spin-summed density here.  Mayer order is therefore withheld
    for open-shell results instead of silently presenting a closed-shell proxy.
    """
    result = _as_result(result)
    _validate_mapping(mol, result)
    result, _, s, p, ps, _, atom, _, _, _, _ = _basis(result)
    wbo = np.asarray(result.wiberg_bond_orders, dtype=float)
    if wbo.shape != (mol.GetNumAtoms(), mol.GetNumAtoms()):
        raise ValueError('Wiberg matrix does not match molecule')
    unavailable = ('spin_state_unknown' if result.number_of_unpaired_electrons is None else
                   'spin_density_unavailable_for_open_shell' if result.number_of_unpaired_electrons else None)
    rows = []
    for index, a, b in _bond_atom_pairs(mol):
        ia, ib = np.flatnonzero(atom == a), np.flatnonzero(atom == b)
        value = None if unavailable else float(np.sum(ps[np.ix_(ia, ib)] * ps[np.ix_(ib, ia)].T))
        rows.append({'bond_index': index, 'begin_atom_index': a, 'end_atom_index': b,
                     'mbo': value, 'wbo': float(wbo[a, b]),
                     'mayer_missing_reason': unavailable})
    return {'bonds': rows, 'units': {'mbo': 'dimensionless', 'wbo': 'dimensionless'},
            'invariance': 'invariant to AO permutation and orbital phase'}


def CalcFrontierBondDescriptors(mol, result, *, degeneracy_ev=1e-3,
                                   occupation_threshold=1e-7, occupation_tolerance=1e-6):
    """Mulliken overlap populations of the rotation-invariant HOMO/LUMO groups."""
    result = _as_result(result)
    _validate_mapping(mol, result)
    _nonnegative(degeneracy_ev, "degeneracy_ev")
    _nonnegative(occupation_tolerance, "occupation_tolerance")
    result, c, s, _, _, _, atom, _, _, _, _ = _basis(result)
    energies, occ, occupied, virtual = _orbitals(result, occupation_threshold)
    rows = []
    missing = {}
    for label, indices in (('homo', _frontier_group(energies, occ, occupied, virtual, True, degeneracy_ev, occupation_tolerance)),
                           ('lumo', _frontier_group(energies, occ, occupied, virtual, False, degeneracy_ev, occupation_tolerance))):
        if not indices:
            missing[label] = 'frontier_group_unavailable'
        # The projector, rather than individual MO coefficients, is invariant
        # under phase changes and rotations in a degenerate group.
        pg = c[:, indices] @ c[:, indices].T / len(indices) if indices else None
        for bond, a, b in _bond_atom_pairs(mol):
            ia, ib = np.flatnonzero(atom == a), np.flatnonzero(atom == b)
            overlap = float(np.sum(pg[np.ix_(ia, ib)] * s[np.ix_(ia, ib)]) * 2.0) if pg is not None else None
            rows.append({'frontier_group': label, 'orb_idx': indices,
                         'bond_index': bond, 'begin_atom_index': a, 'end_atom_index': b,
                         'mul_ovlp_e': overlap,
                         'missing': {'mul_ovlp_e': missing[label]} if not indices else {}})
    return {'bonds': rows, 'units': {'mul_ovlp_e': 'e'},
            'invariance': 'degenerate groups use a projector; orbital phase and group rotation invariant',
            'missing': missing}


def _normal_for_atom(mol, atom_index, conf_id=-1, plane_tol_a=0.1):
    """Fit the atom's ring planes (union for fused rings), or its neighbor plane."""
    positions = np.asarray(mol.GetConformer(conf_id).GetPositions())
    rings = [r for r in mol.GetRingInfo().AtomRings() if atom_index in r]
    indices = (sorted({i for r in rings for i in r}) if rings else
               [atom_index, *[a.GetIdx() for a in mol.GetAtomWithIdx(atom_index).GetNeighbors()]])
    if len(indices) < 3:
        return None, 'insufficient_plane_neighbors'
    points = positions[indices] - positions[indices].mean(axis=0)
    _, singular, vh = np.linalg.svd(points, full_matrices=False)
    if singular[1] <= max(singular[0], 1.) * 1e-8:
        return None, 'local_plane_degenerate'
    normal = vh[-1]
    if np.max(np.abs(points @ normal)) > plane_tol_a:
        return None, 'local_plane_nonplanar'
    return normal, None


def _pi_options(scheme, plane_tol_a, pi_ref, *, natoms=None):
    """Validate and snapshot pi settings, using JSON-stable reference keys."""
    if scheme not in ('mulliken', 'lowdin'):
        raise ValueError("pi scheme must be 'mulliken' or 'lowdin'")
    _nonnegative(plane_tol_a, 'pi_plane_tol_a')
    if pi_ref is not None and not isinstance(pi_ref, Mapping):
        raise ValueError('pi_ref must map original atom indices to reference populations')
    references = {}
    for key, value in (pi_ref or {}).items():
        # Decimal string keys allow analysis_settings to be reused after JSON export.
        index = int(key) if isinstance(key, str) and key.isascii() and key.isdecimal() else key
        if (isinstance(index, bool) or not isinstance(index, (int, np.integer))
                or index < 0 or (isinstance(key, str) and str(index) != key)
                or (natoms is not None and index >= natoms)):
            raise ValueError('pi_ref keys must be original atom indices within the molecule')
        if str(index) in references:
            raise ValueError('pi_ref contains duplicate atom indices')
        _nonnegative(value, 'pi_ref population')
        if value > 2: raise ValueError('pi_ref population must be at most two electrons')
        references[str(index)] = float(value)
    return {'scheme': scheme, 'plane_tol_a': float(plane_tol_a), 'pi_ref': references}


def CalcPiDescriptors(mol, result, *, scheme='lowdin', plane_tol_a=0.1, pi_ref=None):
    """Local p projections; populations in e, closed-shell projected bond indices.

    Mulliken uses sym(P S), Löwdin uses D=S**(1/2) P S**(1/2).
    The pi bond index is (vA.T P S vB)(vB.T P S vA) for Mulliken,
    or (vA.T D vB)**2 for Löwdin. Spin-resolved bond indices are unavailable.
    Neutral, nonradical sp2 carbon has a reference pi population of 1 e;
    other reference populations must be supplied explicitly in pi_ref.
    """
    settings = _pi_options(scheme, plane_tol_a, pi_ref, natoms=mol.GetNumAtoms())
    result = _as_result(result)
    _validate_mapping(mol, result)
    result, _, _, _, ps, lowdin, atom, shell, _, shell_l, _ = _basis(result)
    references = {int(key): value for key, value in settings['pi_ref'].items()}
    directions = None
    if len(result.ao_p_directions):
        directions = _finite(result.ao_p_directions, 'AO p directions', 2)
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
            normal, reason = _normal_for_atom(mol, i, conf_id, plane_tol_a)
        indices = np.flatnonzero((atom == i) & (shell_l[shell] == 1))
        if reason is None and len(indices) != 3:
            reason = 'requires_one_complete_p_shell'
        if reason is None and directions is None:
            reason = 'p_component_directions_unavailable'
        if reason is not None:
            row['missing'] = {'pi_n_e': reason, 'pi_q_e': reason}
        else:
            v = np.zeros(len(atom)); v[indices] = directions[indices] @ normal
            vectors[i] = v
            population = float(v @ population_metric @ v)
            row.update(pi_n_e=population, normal=normal.tolist())
            reference = references.get(i)
            if reference is None and a.GetAtomicNum() == 6 and a.GetFormalCharge() == 0 and a.GetNumRadicalElectrons() == 0:
                reference = 1.
            if reference is None:
                row['missing']['pi_q_e'] = 'pi_reference_unspecified'
            else:
                row['pi_q_e'] = float(reference - population)
        if row['missing']: missing['atoms'][str(i)] = dict(row['missing'])
        atoms.append(row)
    for index, a, b in _bond_atom_pairs(mol):
        reason = None
        if a not in vectors or b not in vectors:
            reason = 'endpoint_pi_projection_unavailable'
        elif result.number_of_unpaired_electrons is None:
            reason = 'spin_state_unknown'
        elif result.number_of_unpaired_electrons:
            reason = 'spin_density_unavailable_for_open_shell'
        value = None
        if reason is None:
            va, vb = vectors[a], vectors[b]
            value = float((va @ ps @ vb) * (vb @ ps @ va) if scheme == 'mulliken' else (va @ lowdin @ vb)**2)
        row = {'bond_index': index, 'begin_atom_index': a, 'end_atom_index': b,
               'pi_bo': value, 'missing': {'pi_bo': reason} if reason else {}}
        bonds.append(row)
        if reason: missing['bonds'][str(index)] = dict(row['missing'])
    return {'atoms': atoms, 'bonds': bonds, 'missing': missing, **settings,
            'charge_reference': '1 e for neutral nonradical sp2 carbon; otherwise explicit pi_ref',
            'units': {'pi_n_e': 'e', 'pi_q_e': 'e', 'pi_bo': 'dimensionless'},
            'invariance': 'AO permutation/sign changes with component directions; rigid motion; local normal sign'}


def CalcDensityDescriptors(mol, result, *, include_frontier_bonds=True, include_pi=True,
                           pi_scheme='lowdin', pi_plane_tol_a=0.1, pi_ref=None):
    """Aggregate all zero-cost descriptors available from one retained result."""
    pi_options = _pi_options(pi_scheme, pi_plane_tol_a, pi_ref, natoms=mol.GetNumAtoms())
    out = {'density': CalcDensityMatrix(result), 'shell_populations': CalcShellPopulations(result),
           'mayer_bonds': CalcMayerBondDescriptors(mol, result)}
    if include_frontier_bonds: out['frontier_bonds'] = CalcFrontierBondDescriptors(mol, result)
    if include_pi: out['pi'] = CalcPiDescriptors(mol, result, **pi_options)
    return out

