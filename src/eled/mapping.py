"""Map results only to the indexed structure and geometry used for calculation."""

from __future__ import annotations

import json
import math
from typing import Any

from .input import _geometry_hash, _ordered_structure_hash, _xyz_from_mol
from .result import XTBResults
from .xtb import CalcXTBDescriptors


def _validate_mapping(mol: Any, result: XTBResults) -> None:
    metadata = result.metadata
    if not metadata.get('ordered_structure_sha256') or not metadata.get('geometry_sha256'):
        raise ValueError('result lacks input provenance; recalculate with XTBCalculator before mapping')
    if metadata.get('atom_count', mol.GetNumAtoms()) != mol.GetNumAtoms():
        raise ValueError('result atom count disagrees with the indexed molecule')
    if len(result.atomicCharges) not in (0, mol.GetNumAtoms()):
        raise ValueError('partial charges do not match the indexed molecule atom count')
    if result.cm5Charges is not None and len(result.cm5Charges) not in (0, mol.GetNumAtoms()):
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
    from .electronic_descriptors import _as_result
    if not isinstance(copy, bool) or not isinstance(privateProps, bool):
        raise ValueError('copy and privateProps must be bool')
    result = _as_result(result)

    _validate_mapping(mol, result)
    natoms = mol.GetNumAtoms()
    charges = tuple(float(v) for v in result.atomicCharges)
    matrix = tuple(tuple(float(v) for v in row) for row in result.wibergBondOrders)
    if charges and len(charges) != natoms:
        raise ValueError('xTB atomic charges do not match the RDKit atom count')
    if matrix and (len(matrix) != natoms or any(len(row) != natoms for row in matrix)):
        raise ValueError('xTB Wiberg bond-order matrix does not match the RDKit atom count')
    if not all(math.isfinite(v) for v in charges + tuple(v for row in matrix for v in row)):
        raise ValueError('xTB charges and bond orders must be finite')
    if matrix and any(abs(matrix[i][j] - matrix[j][i]) > 1e-8 for i in range(natoms) for j in range(i)):
        raise ValueError('xTB bond-order matrix must be symmetric')
    scalars = CalcXTBDescriptors(result)
    if any(isinstance(v, float) and not math.isfinite(v) for v in scalars.values()):
        raise ValueError('xTB scalar properties must be finite')
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
        _set_mol_property(atom, charge_prop, charges[i] if charges else None)
        _set_mol_property(
            atom,
            wiberg_prop,
            sum(value for j, value in enumerate(matrix[i]) if i != j) if matrix else None,
        )
    for bond in target.GetBonds():
        value = matrix[bond.GetBeginAtomIdx()][bond.GetEndAtomIdx()] if matrix else None
        _set_mol_property(bond, bond_prop, value)
    for name, value in scalars.items():
        _set_mol_property(target, prop_name(name), value)
    target.SetProp(prop_name('xtb_metadata_json'), metadata_json)
    return target
