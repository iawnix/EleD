"""Canonical atom, bond, ring and reaction-site selections."""
from __future__ import annotations

from dataclasses import dataclass
from numbers import Integral
from typing import Any, Mapping

from .electronic_descriptors import _indices


@dataclass(frozen=True)
class SelectionSet:
    """Validated selections in the original RDKit atom order."""

    atoms: dict[str, tuple[int, ...]]
    bonds: dict[str, tuple[tuple[int, int], ...]]
    rings: dict[str, tuple[int, ...]]
    reaction_sites: dict[str, tuple[int, ...]]

    def to_dict(self) -> dict[str, dict[str, list]]:
        return {
            'atoms': {name: list(values) for name, values in self.atoms.items()},
            'bonds': {name: [list(pair) for pair in values] for name, values in self.bonds.items()},
            'rings': {name: list(values) for name, values in self.rings.items()},
            'reaction_sites': {name: list(values) for name, values in self.reaction_sites.items()},
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
    for raw_name, raw_indices in _named_mapping(values, field).items():
        name = _name(raw_name, field[:-1] if field.endswith('s') else field)
        indices = _indices(raw_indices, mol.GetNumAtoms(), f'{field[:-1]} {name!r}')
        result[name] = tuple(map(int, indices))
    return result


def NormalizeSelections(mol, selections: Mapping[str, Any] | SelectionSet | None = None) -> SelectionSet:
    """Validate one shared selection schema for fragment and local analyses.

    The mapping may contain ``atoms``, ``bonds``, ``rings`` and
    ``reaction_sites``. Values are named mappings to original atom indices.
    ``rings`` are ordered closed paths. ``bonds`` are atom pairs.
    """
    if isinstance(selections, SelectionSet):
        # Selection validity depends on the molecule being analyzed. Re-run
        # the inexpensive index/topology checks when a reusable selection set
        # is applied to another molecule.
        selections = selections.to_dict()
    if selections is None:
        selections = {}
    if not isinstance(selections, Mapping):
        raise TypeError('selections must be a mapping or SelectionSet')
    unknown = set(selections) - {'atoms', 'bonds', 'rings', 'reaction_sites'}
    if unknown:
        raise ValueError(f'unknown selection fields: {sorted(unknown)}')
    return SelectionSet(
        atoms=_atoms(mol, selections.get('atoms'), 'atoms'),
        bonds=_bonds(mol, selections.get('bonds')),
        rings=_rings(mol, selections.get('rings')),
        reaction_sites=_atoms(mol, selections.get('reaction_sites'), 'reaction_sites'),
    )
