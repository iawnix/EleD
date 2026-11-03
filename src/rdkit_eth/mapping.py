"""Map results only to the indexed structure and geometry used for calculation."""

from __future__ import annotations
import json
import math
from typing import Any

from .result import XTBResult
from .input import _ordered_structure_hash, _geometry_hash, _xyz_from_mol
from .xtb import calc_xtb_descriptors


def _validate_mapping(mol: Any, result: XTBResult) -> None:
    metadata = result.metadata
    if not metadata.get('ordered_structure_sha256') or not metadata.get('geometry_sha256'):
        raise ValueError('result lacks input provenance; recalculate with run_xtb before mapping')
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


def annotate_mol(
    mol: Any,
    result: XTBResult,
    *,
    copy: bool = True,
    atom_charge_prop: str = 'xtb_partial_charge',
    bond_order_prop: str = 'xtb_wiberg_order',
    atom_wiberg_prop: str = 'xtb_wiberg_sum',
) -> Any:
    """Return an annotated copy; ``copy=False`` writes to the validated input.

    The original atom order, indexed topology and input conformer must still
    match. Existing RDKit bond types are preserved; WBOs are numeric properties
    on existing bonds. Nonbonded pair values remain in the result matrix.
    """
    from rdkit import Chem
    from .electronic_descriptors import _as_result
    result = _as_result(result)

    _validate_mapping(mol, result)
    natoms = mol.GetNumAtoms()
    charges = tuple(float(v) for v in result.partial_charges)
    matrix = tuple(tuple(float(v) for v in row) for row in result.wiberg_bond_orders)
    if charges and len(charges) != natoms:
        raise ValueError('xTB atomic charges do not match the RDKit atom count')
    if matrix and (len(matrix) != natoms or any(len(row) != natoms for row in matrix)):
        raise ValueError('xTB Wiberg bond-order matrix does not match the RDKit atom count')
    if not all(math.isfinite(v) for v in charges + tuple(v for row in matrix for v in row)):
        raise ValueError('xTB charges and bond orders must be finite')
    if matrix and any(abs(matrix[i][j] - matrix[j][i]) > 1e-8 for i in range(natoms) for j in range(i)):
        raise ValueError('xTB bond-order matrix must be symmetric')
    scalars = calc_xtb_descriptors(result)
    if any(isinstance(v, float) and not math.isfinite(v) for v in scalars.values()):
        raise ValueError('xTB scalar properties must be finite')
    metadata_json = json.dumps(dict(result.metadata), sort_keys=True, allow_nan=False)
    if not atom_charge_prop or not bond_order_prop or not atom_wiberg_prop:
        raise ValueError('property names must be nonempty')
    if atom_charge_prop == atom_wiberg_prop:
        raise ValueError('atomic charge and Wiberg property names must differ')
    target = Chem.Mol(mol) if copy else mol
    for i, atom in enumerate(target.GetAtoms()):
        _set_mol_property(atom, atom_charge_prop, charges[i] if charges else None)
        _set_mol_property(atom, atom_wiberg_prop,
                          sum(value for j, value in enumerate(matrix[i]) if i != j) if matrix else None)
    for bond in target.GetBonds():
        value = matrix[bond.GetBeginAtomIdx()][bond.GetEndAtomIdx()] if matrix else None
        _set_mol_property(bond, bond_order_prop, value)
    for name, value in scalars.items():
        _set_mol_property(target, name, value)
    target.SetProp('xtb_metadata_json', metadata_json)
    return target
