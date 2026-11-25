"""Canonical atom, bond, ring and reaction-site selections."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping

from .electronic_descriptors import _indices


@dataclass(frozen=True)
class SelectionSet:
    """Validated selections in the original RDKit atom order."""

    atoms: dict[str, tuple[int, ...]]
    bonds: dict[str, tuple[tuple[int, int], ...]]
    rings: dict[str, tuple[int, ...]]
    reactionSites: dict[str, tuple[int, ...]]

    def ToDict(self) -> dict[str, dict[str, list]]:
        return {
            'atoms': {name: list(values) for name, values in self.atoms.items()},
            'bonds': {name: [list(pair) for pair in values] for name, values in self.bonds.items()},
            'rings': {name: list(values) for name, values in self.rings.items()},
            'reaction_sites': {name: list(values) for name, values in self.reactionSites.items()},
        }


def _named_mapping(values: Any, field: str) -> Mapping:
    if values is None:
        return {}
    if not isinstance(values, Mapping) or not values:
        raise ValueError(f'{field} must be a nonempty mapping when supplied')
    return values


def _name(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{field} names must be nonempty strings')
    return value


def _rings(mol, values):
    result = {}
    for raw_name, raw_indices in _named_mapping(values, 'rings').items():
        name = _name(raw_name, 'ring')
        indices = tuple(int(value) for value in _indices(raw_indices, mol.GetNumAtoms(), f'ring {name!r}'))
        if len(indices) < 3:
            raise ValueError('rings require at least three atoms')
        for first, second in zip(indices, indices[1:] + indices[:1]):
            if mol.GetBondBetweenAtoms(first, second) is None:
                raise ValueError(f'ring {name!r} must be an ordered bonded closed path')
        result[name] = indices
    return result


def _bonds(mol, values):
    result = {}
    for raw_name, raw_pairs in _named_mapping(values, 'bonds').items():
        name = _name(raw_name, 'bond')
        pairs = []
        for pair in raw_pairs:
            if len(pair) != 2:
                raise ValueError(f'bond selection {name!r} requires atom pairs')
            indices = _indices(pair, mol.GetNumAtoms(), f'bond selection {name!r}')
            first, second = map(int, indices)
            if first == second or mol.GetBondBetweenAtoms(first, second) is None:
                raise ValueError(f'bond selection {name!r} contains a nonbonded or self pair')
            pairs.append((first, second))
        if not pairs:
            raise ValueError(f'bond selection {name!r} cannot be empty')
        result[name] = tuple(pairs)
    return result


def _atoms(mol, values, field):
    result = {}
    label = field[:-1] if field.endswith('s') else field
    for raw_name, raw_indices in _named_mapping(values, field).items():
        name = _name(raw_name, label)
        indices = _indices(raw_indices, mol.GetNumAtoms(), f'{label} {name!r}')
        result[name] = tuple(map(int, indices))
    return result


def NormalizeSelections(mol, selections: Mapping[str, Any] | SelectionSet | None = None) -> SelectionSet:
    """Validate one shared selection schema for fragment and local analyses.

    The mapping may contain ``atoms``, ``bonds``, ``rings`` and
    ``reactionSites`` (or JSON wire key ``reaction_sites``).
    """
    if isinstance(selections, SelectionSet):
        selections = selections.ToDict()
    if selections is None:
        selections = {}
    if not isinstance(selections, Mapping):
        raise TypeError('selections must be a mapping or SelectionSet')
    unknown = set(selections) - {'atoms', 'bonds', 'rings', 'reactionSites', 'reaction_sites'}
    if unknown:
        raise ValueError(f'unknown selection fields: {sorted(unknown)}')
    if 'reactionSites' in selections and 'reaction_sites' in selections:
        raise ValueError("specify only one of 'reactionSites' or 'reaction_sites'")
    raw_sites = selections.get('reactionSites', selections.get('reaction_sites'))
    return SelectionSet(
        atoms=_atoms(mol, selections.get('atoms'), 'atoms'),
        bonds=_bonds(mol, selections.get('bonds')),
        rings=_rings(mol, selections.get('rings')),
        reactionSites=_atoms(mol, raw_sites, 'reactionSites'),
    )
