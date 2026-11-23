"""Same-basis orbital population analysis in a nonorthogonal AO basis.

These are ground-state MO descriptors, not transition densities or excited
states. Native coefficients alone do not supply the required overlap metric.
"""
from __future__ import annotations

import numpy as np

from .electronic_descriptors import _orbital_basis_arrays, _has_population_basis, _charges, _as_result, _fragments, _indices, _nonnegative, _orbitals
from .errors import DescriptorUnavailableError


def CalcOrbitalPopulations(result, *, orb_idx=None, scheme='lowdin',
                             orthonormality_tolerance=1e-6):
    """Return atom-by-orbital fractions and occupation-weighted populations.

    Lowdin: sum |S^(1/2) C|^2 on each atom (nonnegative).
    Mulliken: sum C * (S C) on each atom (signed; may be negative).
    Columns are divided by their measured C.T S C norm after validation;
    orbital populations sum to one, electron populations include occupations.
    Requires overlap and AO mapping from the SAME calculation as C. No engine
    is invoked here, and matrices from another backend are never substituted.
    """
    result = _as_result(result)
    if scheme not in {'lowdin', 'mulliken'}:
        raise ValueError("scheme must be 'lowdin' or 'mulliken'")
    _nonnegative(orthonormality_tolerance, 'orthonormality_tolerance')
    if not _has_population_basis(result):
        raise DescriptorUnavailableError(
            'orbital populations require same-calculation overlap and AO-to-atom mapping; '
            'use XTBCalculator(...).calculate(..., keep_overlap=True).')
    energies, occupations, _, _ = _orbitals(result)
    n = len(energies)
    natoms = len(_charges(result))
    coefficients, overlap, mapping = _orbital_basis_arrays(result, n, natoms)
    eigenvalues, eigenvectors = np.linalg.eigh(overlap)
    if eigenvalues.min() <= max(1., eigenvalues.max()) * 1e-10:
        raise ValueError('overlap must be positive definite and well conditioned')
    metric = coefficients.T @ overlap @ coefficients
    error = float(np.max(np.abs(metric - np.eye(n))))
    if error > orthonormality_tolerance or np.any(np.diag(metric) <= 0):
        raise ValueError(f'coefficients are not orthonormal in the supplied overlap: error={error:g}')
    selected = np.arange(n) if orb_idx is None else _indices(orb_idx, n, 'orb_idx')
    c = coefficients[:, selected]
    if scheme == 'lowdin':
        root = (eigenvectors * np.sqrt(eigenvalues)) @ eigenvectors.T
        ao = (root @ c) ** 2
    else:
        ao = c * (overlap @ c)
    ao /= np.diag(metric)[selected]
    atom = np.zeros((natoms, len(selected)))
    np.add.at(atom, mapping, ao)
    return {
        'scheme': scheme, 'axes': ['atom', 'orbital'],
        'orb_idx': selected.tolist(),
        'orbital_occupations': occupations[selected].tolist(),
        'atom_pop': atom.tolist(),
        'mo_n_e': (atom @ occupations[selected]).tolist(),
        'orthonormality_error': error,
        'normalization': 'unit norm per MO; electron populations include occupations',
    }


def CalcFrontierDescriptors(mol, result, *, fragments=None,
                              occupation_threshold=1e-7, degeneracy_ev=1e-3,
                              occupation_tolerance=1e-6):
    """Lowdin HOMO/LUMO group descriptors, with atom-centroid spatial proxies.

    Frontier groups include contiguous near-degenerate orbitals with matching
    occupations in the same occupied/virtual class. Each group is equally
    averaged and normalized to one, making the result invariant to rotations
    within that subspace. IPR=sum(p_atom**2); participation ratio=1/IPR.
    Entropy uses natural logarithms. Centroids weight atom positions by these
    group populations; their separation is NOT an excited-state electron-hole
    distance. Atom overlap is sum(sqrt(p_H*p_L)), not an orbital overlap integral.
    """
    from .mapping import _validate_mapping
    result = _as_result(result)
    _validate_mapping(mol, result)
    _nonnegative(degeneracy_ev, 'degeneracy_ev')
    _nonnegative(occupation_tolerance, 'occupation_tolerance')
    energies, occupations, occupied, virtual = _orbitals(result, occupation_threshold)
    populations = CalcOrbitalPopulations(result, scheme='lowdin')
    atom = np.asarray(populations['atom_pop'])
    groups = _fragments(fragments, mol.GetNumAtoms()) if fragments is not None else {}
    positions = np.asarray(mol.GetConformer(result.metadata['conformer_id']).GetPositions())

    def frontier(indices, last):
        if not len(indices):
            return None
        selected = _frontier_group(energies, occupations, occupied, virtual,
                                   last, degeneracy_ev, occupation_tolerance)
        weights = atom[:, selected].mean(axis=1)
        ipr = float(weights @ weights)
        positive = weights[weights > 0]
        return {
            'orb_idx': selected, 'n_orb': len(selected),
            'emin_ev': float(energies[selected].min()),
            'emax_ev': float(energies[selected].max()),
            'atom_pop': weights.tolist(),
            'frag_pop': {name: float(weights[indices].sum()) for name, indices in groups.items()},
            'ipr': ipr, 'pr': 1. / ipr,
            'ent_nat': float(-(positive * np.log(positive)).sum()),
            'cent_a': (weights @ positions).tolist(),
        }

    homo, lumo = frontier(occupied, True), frontier(virtual, False)
    overlap = separation = None
    if homo is not None and lumo is not None:
        overlap = float(np.sqrt(np.asarray(homo['atom_pop']) * lumo['atom_pop']).sum())
        separation = float(np.linalg.norm(np.asarray(homo['cent_a']) - lumo['cent_a']))
    assigned = set(np.concatenate(list(groups.values())).tolist()) if groups else set()
    return {
        'scheme': 'lowdin', 'occupation_threshold': float(occupation_threshold),
        'degeneracy_ev': float(degeneracy_ev), 'occupation_tolerance': float(occupation_tolerance),
        'homo': homo, 'lumo': lumo,
        'hl_ovlp': overlap,
        'hl_dist_a': separation,
        'unassigned_atom_indices': sorted(set(range(mol.GetNumAtoms())) - assigned),
        'orthonormality_error': populations['orthonormality_error'],
        'interpretation': 'ground-state frontier subspace populations and atom-centroid proxies; not excitation or transition-density descriptors',
    }


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
