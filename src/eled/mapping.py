"""Map results only to the indexed structure and geometry used for calculation."""

from __future__ import annotations

import json
from typing import Any

from .input import _geometry_hash, _ordered_structure_hash, _xyz_from_mol
from .electronic_descriptors import CalcXTBDescriptors
from .result import XTBResults


def _validate_mapping(mol: Any, result: XTBResults) -> None:
    metadata = result.metadata
    if not metadata.get('ordered_structure_sha256') or not metadata.get('geometry_sha256'):
        raise ValueError('result lacks input provenance; recalculate with XTBCalculator before mapping')
    if metadata.get('atom_count', mol.GetNumAtoms()) != mol.GetNumAtoms():
        raise ValueError('result atom count disagrees with the indexed molecule')
    if len(result.atomicCharges) not in (0, mol.GetNumAtoms()):
        raise ValueError('partial charges do not match the indexed molecule atom count')
    if len(result.cm5Charges) not in (0, mol.GetNumAtoms()):
        raise ValueError(f'CM5 charges do not match atom count: expected {mol.GetNumAtoms()}, got {len(result.cm5Charges)}')
    if metadata['ordered_structure_sha256'] != _ordered_structure_hash(mol):
        raise ValueError('molecule structure or atom order does not match the xTB input')
    try:
        xyz = _xyz_from_mol(mol, metadata['conformer_id'])
    except (ValueError, KeyError) as exc:
        raise ValueError('the original xTB input conformer is missing or invalid') from exc
    if _geometry_hash(xyz) != metadata['geometry_sha256']:
        raise ValueError('molecule coordinates or atom order do not match the xTB input')


def _set_mol_property(mol: Any, name: str, value: Any) -> None:
    if value is None:
        if mol.HasProp(name):
            mol.ClearProp(name)
    elif isinstance(value, int):
        mol.SetIntProp(name, value)
    elif isinstance(value, float):
        mol.SetDoubleProp(name, value)
    else:
        mol.SetProp(name, str(value))


def AnnotateMol(
    mol: Any,
    result: XTBResults,
    *,
    copy: bool = True,
    privateProps: bool = False,
    atomChargeProp: str = 'xtb_partial_charge',
    bondOrderProp: str = 'xtb_wbo',
    atomWibergProp: str = 'xtb_wbo_sum',
) -> Any:
    """Return an annotated copy; ``copy=False`` writes to the validated input.

    Set ``privateProps=True`` to prefix property names with ``_`` (matching
    RDKit's ``_GasteigerCharge`` convention so properties are omitted from
    default ``GetPropsAsDict()`` / ``GetPropNames()`` inspection unless
    ``includePrivate=True``).
    """
    from rdkit import Chem
    from .electronic_descriptors import _as_result, _charges, _wbo
    if not isinstance(copy, bool) or not isinstance(privateProps, bool):
        raise ValueError('copy and privateProps must be bool')
    result = _as_result(result)

    _validate_mapping(mol, result)
    # _charges/_wbo check finiteness, atom count and symmetry; _validate_mapping
    # already matched the atom count and geometry.
    charges = _charges(result)
    matrix = _wbo(result, mol.GetNumAtoms())
    scalars = CalcXTBDescriptors(result)
    metadata_json = json.dumps(dict(result.metadata), sort_keys=True, allow_nan=False)
    if not atomChargeProp or not bondOrderProp or not atomWibergProp:
        raise ValueError('property names must be nonempty')

    def prop_name(name: str) -> str:
        return f'_{name}' if privateProps and not name.startswith('_') else name

    charge_prop = prop_name(atomChargeProp)
    bond_prop = prop_name(bondOrderProp)
    wiberg_prop = prop_name(atomWibergProp)
    if charge_prop == wiberg_prop:
        raise ValueError('atomic charge and Wiberg property names must differ')
    target = Chem.Mol(mol) if copy else mol
    for i, atom in enumerate(target.GetAtoms()):
        _set_mol_property(atom, charge_prop, float(charges[i]) if len(charges) else None)
        _set_mol_property(atom, wiberg_prop, float(matrix[i].sum()) if matrix is not None else None)
    for bond in target.GetBonds():
        value = float(matrix[bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()]) if matrix is not None else None
        _set_mol_property(bond, bond_prop, value)
    for name, value in scalars.items():
        _set_mol_property(target, prop_name(name), value)
    target.SetProp(prop_name('xtb_metadata_json'), metadata_json)
    return target
