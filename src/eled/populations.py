"""Same-basis orbital population analysis in a nonorthogonal AO basis.

These are ground-state MO descriptors, not transition densities or excited
states. Native coefficients alone do not supply the required overlap metric.
"""
from __future__ import annotations

import numpy as np

from .electronic_descriptors import _as_result, _finite, _fragments, _indices, _nonnegative, _orbitals
from .errors import DescriptorUnavailableError
from .result import XTBResult


def CalcOrbitalPopulations(result, *, orbital_indices=None, scheme='lowdin',
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
    if (not isinstance(result, XTBResult) or not len(result.overlap_matrix)
            or not len(result.ao_atom_indices)):
        raise DescriptorUnavailableError(
            'orbital populations require same-calculation overlap and AO-to-atom mapping; '
            'use XTBCalculator(...).calculate(..., keep_overlap=True) with the verified native basis extension. '
            'Official libraries support basic descriptors without this extension.')
    energies, occupations, _, _ = _orbitals(result)
    n = len(energies)
    coefficients = _finite(result.orbital_coefficients, 'orbital coefficients', 2)
    overlap = _finite(result.overlap_matrix, 'overlap', 2)
    charges = _finite(result.partial_charges, 'partial charges', 1)
    natoms = len(charges)
    if result.metadata.get('atom_count', natoms) != natoms:
        raise ValueError('partial charges do not match atom count')
    if not n or not natoms or coefficients.shape != (n, n) or overlap.shape != (n, n):
        raise ValueError('coefficients and overlap must match the full square AO/MO basis')
    mapping = np.asarray(result.ao_atom_indices)
    if (mapping.shape != (n,) or mapping.dtype.kind not in 'iu'
            or np.any(mapping < 0) or np.any(mapping >= natoms)
            or not np.array_equal(np.unique(mapping), np.arange(natoms))):
        raise ValueError('AO-to-atom mapping must cover every atom with valid integer indices')
    if not np.allclose(overlap, overlap.T, atol=1e-10, rtol=0):
        raise ValueError('overlap must be symmetric')
    eigenvalues, eigenvectors = np.linalg.eigh(overlap)
    if eigenvalues.min() <= max(1., eigenvalues.max()) * 1e-10:
        raise ValueError('overlap must be positive definite and well conditioned')
    metric = coefficients.T @ overlap @ coefficients
    error = float(np.max(np.abs(metric - np.eye(n))))
    if error > orthonormality_tolerance or np.any(np.diag(metric) <= 0):
        raise ValueError(f'coefficients are not orthonormal in the supplied overlap: error={error:g}')
    selected = np.arange(n) if orbital_indices is None else _indices(orbital_indices, n, 'orbital_indices')
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
        'orbital_indices': selected.tolist(),
        'orbital_occupations': occupations[selected].tolist(),
        'atom_populations': atom.tolist(),
        'electron_populations_e': (atom @ occupations[selected]).tolist(),
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
    atom = np.asarray(populations['atom_populations'])
    groups = _fragments(fragments, mol.GetNumAtoms()) if fragments is not None else {}
    positions = np.asarray(mol.GetConformer(result.metadata['conformer_id']).GetPositions())

    def frontier(indices, last):
        if not len(indices):
            return None
        anchor = int(indices[-1] if last else indices[0])
        # Anchor is the occupied group's high-energy edge or virtual low edge.
        selected = [anchor]
        candidates = indices[-2::-1] if last else indices[1:]
        for index in candidates:
            if (abs(energies[index] - energies[anchor]) > degeneracy_ev
                    or abs(occupations[index] - occupations[anchor]) > occupation_tolerance
                    or abs(int(index) - selected[-1]) != 1):
                break
            selected.append(int(index))
        selected.sort()
        weights = atom[:, selected].mean(axis=1)
        ipr = float(weights @ weights)
        positive = weights[weights > 0]
        return {
            'orbital_indices': selected, 'num_orbitals': len(selected),
            'energy_min_ev': float(energies[selected].min()),
            'energy_max_ev': float(energies[selected].max()),
            'atom_populations': weights.tolist(),
            'fragment_populations': {name: float(weights[indices].sum()) for name, indices in groups.items()},
            'inverse_participation_ratio': ipr, 'participation_ratio': 1. / ipr,
            'entropy_nats': float(-(positive * np.log(positive)).sum()),
            'atom_centroid_angstrom': (weights @ positions).tolist(),
        }

    homo, lumo = frontier(occupied, True), frontier(virtual, False)
    overlap = separation = None
    if homo is not None and lumo is not None:
        overlap = float(np.sqrt(np.asarray(homo['atom_populations']) * lumo['atom_populations']).sum())
        separation = float(np.linalg.norm(np.asarray(homo['atom_centroid_angstrom']) - lumo['atom_centroid_angstrom']))
    assigned = set(np.concatenate(list(groups.values())).tolist()) if groups else set()
    return {
        'scheme': 'lowdin', 'occupation_threshold': float(occupation_threshold),
        'degeneracy_ev': float(degeneracy_ev), 'occupation_tolerance': float(occupation_tolerance),
        'homo': homo, 'lumo': lumo,
        'homo_lumo_atom_overlap': overlap,
        'homo_lumo_atom_centroid_separation_angstrom': separation,
        'unassigned_atom_indices': sorted(set(range(mol.GetNumAtoms())) - assigned),
        'orthonormality_error': populations['orthonormality_error'],
        'interpretation': 'ground-state frontier subspace populations and atom-centroid proxies; not excitation or transition-density descriptors',
    }
