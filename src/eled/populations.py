"""Same-basis orbital population analysis in a nonorthogonal AO basis.

These are ground-state MO descriptors, not transition densities or excited
states. Native coefficients alone do not supply the required overlap metric.
"""
from __future__ import annotations

import numpy as np

from .electronic_descriptors import (
    _as_result,
    _fragments,
    _indices,
    _nonnegative,
    _orbital_basis,
    _orbitals,
)


def CalcOrbitalPopulations(result, *, orbIdx=None, scheme='lowdin'):
    """Return atom-by-orbital fractions and occupation-weighted populations.

    Lowdin: sum ``|S^(1/2) C|^2`` on each atom (nonnegative).
    Mulliken: sum ``C * (S C)`` on each atom (signed; may be negative).
    Requires overlap and AO mapping from the SAME calculation as ``C``
    (``XTBCalculator(...).Calculate(..., keepOverlap=True)``).
    """
    result = _as_result(result)
    if scheme not in {'lowdin', 'mulliken'}:
        raise ValueError("scheme must be 'lowdin' or 'mulliken'")
    basis = _orbital_basis(result)
    n = len(basis.energies)
    selected = np.arange(n) if orbIdx is None else _indices(orbIdx, n, 'orbIdx')
    c = basis.coefficients[:, selected]
    sc = basis.overlap @ c
    if scheme == 'lowdin':
        ao = (basis.overlap_half @ c) ** 2
    else:
        ao = c * sc
    ao /= np.einsum('ij,ij->j', c, sc)
    atom = np.zeros((basis.natoms, len(selected)))
    np.add.at(atom, basis.ao_atom, ao)
    return {
        'scheme': scheme,
        'axes': ['atom', 'orbital'],
        'orb_idx': selected.tolist(),
        'orbital_occupations': basis.occupations[selected].tolist(),
        'atom_pop': atom.tolist(),
        'mo_n_e': (atom @ basis.occupations[selected]).tolist(),
        'orthonormality_error': basis.orthonormality_error,
        'normalization': 'unit norm per MO; electron populations include occupations',
    }


def CalcFrontierDescriptors(
    mol,
    result,
    *,
    fragments=None,
    degeneracyEV=1e-3,
    occupationTolerance=1e-6,
):
    """Lowdin HOMO/LUMO group descriptors, with atom-centroid spatial proxies."""
    from .mapping import _validate_mapping
    result = _as_result(result)
    _validate_mapping(mol, result)
    _nonnegative(degeneracyEV, 'degeneracyEV')
    _nonnegative(occupationTolerance, 'occupationTolerance')
    energies, occupations, occupied, virtual = _orbitals(result)
    populations = CalcOrbitalPopulations(result, scheme='lowdin')
    atom = np.asarray(populations['atom_pop'])
    groups = _fragments(fragments, mol.GetNumAtoms()) if fragments is not None else {}
    positions = np.asarray(mol.GetConformer(result.metadata['conformer_id']).GetPositions())

    def frontier(indices, last):
        if not len(indices):
            return None
        selected = _frontier_group(
            energies, occupations, occupied, virtual, last, degeneracyEV, occupationTolerance
        )
        weights = atom[:, selected].mean(axis=1)
        ipr = float(weights @ weights)
        positive = weights[weights > 0]
        return {
            'orb_idx': selected,
            'n_orb': len(selected),
            'emin_ev': float(energies[selected].min()),
            'emax_ev': float(energies[selected].max()),
            'atom_pop': weights.tolist(),
            'frag_pop': {name: float(weights[group_idx].sum()) for name, group_idx in groups.items()},
            'ipr': ipr,
            'pr': 1.0 / ipr,
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
        'scheme': 'lowdin',
        'frontier_definition': 'aufbau electron count',
        'degeneracy_ev': float(degeneracyEV),
        'occupation_tolerance': float(occupationTolerance),
        'homo': homo,
        'lumo': lumo,
        'hl_ovlp': overlap,
        'hl_dist_a': separation,
        'unassigned_atom_indices': sorted(set(range(mol.GetNumAtoms())) - assigned),
        'orthonormality_error': populations['orthonormality_error'],
        'interpretation': (
            'ground-state frontier subspace populations and atom-centroid proxies; '
            'not excitation or transition-density descriptors'
        ),
    }


def _frontier_group(energies, occ, occupied, virtual, last, degeneracyEV, occupationTolerance):
    indices = occupied if last else virtual
    if not len(indices):
        return []
    anchor = int(indices[-1] if last else indices[0])
    selected = [anchor]
    candidates = indices[-2::-1] if last else indices[1:]
    for i in candidates:
        if (abs(energies[i] - energies[anchor]) > degeneracyEV
                or abs(occ[i] - occ[anchor]) > occupationTolerance
                or abs(int(i) - selected[-1]) != 1):
            break
        selected.append(int(i))
    return sorted(selected)
