"""Electronic descriptors derived from a single result, without engine calls.

Indices always refer to the original calculation's atom/MO order. Wiberg
sums include all off-diagonal atom pairs, not just RDKit topology bonds.
"""
from __future__ import annotations

import json
from collections.abc import Mapping
from numbers import Integral

import numpy as np

from .result import XTBResult

AU_DIPOLE_TO_DEBYE = 2.541746473


def _as_result(result):
    """Copy numeric fields from a compiled result; never fabricate provenance."""
    if isinstance(result, XTBResult):
        return result
    # Avoid importing the optional compiled module just for type discovery.
    if (type(result).__module__ != 'eled.rdXTBTools'
            or type(result).__name__ != 'XTBResults'):
        raise TypeError('expected XTBResult or rdXTBTools.XTBResults')
    if not result.success:
        raise ValueError(f'calculation failed: {result.errorMessage}')
    basis = {}
    if result.overlapRetained:
        basis = dict(overlap_matrix=tuple(map(tuple, result.GetOverlapMatrix())),
                     orbital_coefficients=tuple(map(tuple, result.GetOrbitalCoefficients())),
                     ao_atom_indices=tuple(map(int, result.GetAOAtomMap())),
                     ao_shell_indices=tuple(map(int, result.GetAOShellMap())),
                     shell_atom_indices=tuple(map(int, result.GetShellAtomMap())),
                     shell_angular_momenta=tuple(map(int, result.GetShellAngularMomenta())))
    metadata = json.loads(result.inputProvenanceJSON) if result.inputProvenanceJSON else {'atom_count': result.numAtoms}
    metadata = dict(metadata)
    from .scf import scf_record
    metadata['scf'] = scf_record(result.method, getattr(result, 'scfIterations', -1),
        [getattr(result, name, -1.) for name in ('scfEnergyChangeHartree', 'scfMixerRMS',
                                                'scfEnergyThresholdHartree', 'scfMixerThreshold')])
    metadata.update(calculation_status='success', initialization=result.initialization, CalculationSettings={
        'max_iterations': result.maxIterations, 'accuracy': result.accuracy,
        'electronic_temperature': result.electronicTemperature,
        'effective_electronic_temperature': result.effectiveElectronicTemperature},
        engine={'backend': 'compiled', 'path': result.libraryPath, 'sha256': result.librarySHA256,
                'source_commit': result.sourceCommit, 'api_version': result.apiVersion},
        parameter_file={'path': result.parameterFile, 'sha256': result.parameterSHA256} if result.parameterFile else None)
    if hasattr(result, 'GetGradient'):
        metadata['gradient_hartree_per_bohr'] = tuple(float(value) for value in np.asarray(result.GetGradient()).ravel())
    return XTBResult(
        **basis,
        method=result.method, version=None, total_energy_hartree=result.totalEnergyHartree,
        homo_lumo_gap_ev=None, number_of_molecular_orbitals=result.numOrbitals,
        number_of_electrons=result.numElectrons,
        number_of_unpaired_electrons=result.numUnpairedElectrons,
        orbital_energies_ev=tuple(result.GetOrbitalEnergies()),
        orbital_occupations=tuple(result.GetOrbitalOccupancies()),
        partial_charges=tuple(result.GetAtomicCharges()), dipole_e_bohr=tuple(result.GetDipole()),
        wiberg_bond_orders=tuple(map(tuple, result.GetWibergBondOrders())),
        backend='compiled', metadata=metadata,
    )


def _finite(values, name, ndim):
    raw = np.asarray(values)
    if np.iscomplexobj(raw):
        raise ValueError(f'{name} must be real')
    array = np.asarray(values, dtype=float)
    if array.ndim != ndim or not np.isfinite(array).all():
        raise ValueError(f'{name} must be finite with {ndim} dimensions')
    return array


def _nonnegative(value, name):
    if isinstance(value, bool) or not np.isfinite(value) or value < 0:
        raise ValueError(f'{name} must be finite and nonnegative')


def _orbitals(result, occupation_threshold=1e-7):
    _nonnegative(occupation_threshold, 'occupation_threshold')
    energies = _finite(result.orbital_energies_ev, 'orbital energies', 1)
    occupations = _finite(result.orbital_occupations, 'occupations', 1)
    if energies.shape != occupations.shape or np.any(occupations < 0):
        raise ValueError('orbital energies and nonnegative occupations must have matching lengths')
    if np.any(np.diff(energies) < -1e-10):
        raise ValueError('orbital energies must be sorted')
    if result.number_of_molecular_orbitals is not None and len(energies) != result.number_of_molecular_orbitals:
        raise ValueError('orbital count does not match arrays')
    occupied = np.flatnonzero(occupations > occupation_threshold)
    virtual = np.flatnonzero(occupations <= occupation_threshold)
    return energies, occupations, occupied, virtual


def _charges(result):
    charges = _finite(result.partial_charges, 'partial charges', 1)
    expected = result.metadata.get('atom_count')
    if expected is not None and len(charges) not in (0, expected):
        raise ValueError('partial charges do not match atom count')
    return charges


def _wbo(result, natoms):
    if not len(result.wiberg_bond_orders):
        return None
    wbo = _finite(result.wiberg_bond_orders, 'Wiberg matrix', 2)
    if wbo.shape != (natoms, natoms) or not np.allclose(wbo, wbo.T, atol=1e-8, rtol=0):
        raise ValueError('Wiberg matrix must be symmetric and match atom count')
    wbo = wbo.copy()
    np.fill_diagonal(wbo, 0.)
    return wbo


def _indices(values, limit, name):
    values = list(values)
    if (not values or any(isinstance(v, (bool, np.bool_)) or not isinstance(v, Integral)
                          or v < 0 or v >= limit for v in values)
            or len(set(values)) != len(values)):
        raise ValueError(f'{name} must contain unique integer indices in [0, {limit})')
    return np.asarray(values, dtype=int)


def _fragments(fragments, natoms):
    if not isinstance(fragments, Mapping) or not fragments:
        raise ValueError('fragments must be a nonempty mapping of names to atom indices')
    groups, seen = {}, set()
    for name, atoms in fragments.items():
        if not isinstance(name, str) or not name.strip():
            raise ValueError('fragment names must be nonempty strings')
        indices = _indices(atoms, natoms, f'fragment {name!r}')
        if seen.intersection(indices):
            raise ValueError('fragments must be disjoint; atoms cannot belong to multiple groups')
        seen.update(indices)
        groups[name] = indices
    return groups


def CalcXTBDescriptors(result, *, frontier_depth=2, occupation_threshold=1e-7):
    """Return JSON-compatible molecular scalars; unavailable values are None.

    frontier_depth=2 adds HOMO-1/-2 and LUMO+1/+2. Charge standard deviation
    uses ddof=0. Dipoles refer to the supplied coordinate origin (relevant for
    ions); a.u. denotes e*Bohr. A finite-temperature occupation threshold
    defines the frontier convention, not a predicted excitation energy.
    """
    if isinstance(frontier_depth, bool) or not isinstance(frontier_depth, Integral) or frontier_depth < 0:
        raise ValueError('frontier_depth must be a nonnegative integer')
    result = _as_result(result)
    energies, occupations, occupied, virtual = _orbitals(result, occupation_threshold)
    charges = _charges(result)
    dipole = _finite(result.dipole_e_bohr, 'dipole', 1)
    if len(dipole) not in (0, 3):
        raise ValueError('dipole must have three components')
    energy = result.total_energy_hartree
    if energy is not None and not np.isfinite(energy):
        raise ValueError('total energy must be finite')
    homo = float(energies[occupied[-1]]) if len(occupied) else None
    lumo = float(energies[virtual[0]]) if len(virtual) else None
    output = {
        'xtb_method': result.method, 'xtb_backend': result.backend, 'xtb_version': result.version,
        'xtb_total_energy_hartree': energy, 'xtb_homo_ev': homo, 'xtb_lumo_ev': lumo,
        'xtb_gap_ev': lumo - homo if homo is not None and lumo is not None else None,
        'xtb_num_orbitals': result.number_of_molecular_orbitals,
        'xtb_num_electrons': result.number_of_electrons,
        'xtb_num_unpaired_electrons': result.number_of_unpaired_electrons,
        'xtb_occupation_threshold': float(occupation_threshold),
        'xtb_homo_index': int(occupied[-1]) if len(occupied) else None,
        'xtb_lumo_index': int(virtual[0]) if len(virtual) else None,
    }
    for offset in range(1, frontier_depth + 1):
        output[f'xtb_homo_minus_{offset}_ev'] = float(energies[occupied[-1-offset]]) if len(occupied) > offset else None
        output[f'xtb_lumo_plus_{offset}_ev'] = float(energies[virtual[offset]]) if len(virtual) > offset else None
    for i, axis in enumerate('xyz'):
        output[f'xtb_dipole_{axis}_e_bohr'] = float(dipole[i]) if len(dipole) else None
    norm = float(np.linalg.norm(dipole)) if len(dipole) else None
    output['xtb_dipole_norm_e_bohr'] = norm
    output['xtb_dipole_norm_debye'] = norm * AU_DIPOLE_TO_DEBYE if norm is not None else None
    for label, func in [('sum', np.sum), ('min', np.min), ('max', np.max),
                        ('mean', np.mean), ('std', np.std), ('abs_sum', lambda q: np.abs(q).sum())]:
        output[f'xtb_charge_{label}_e'] = float(func(charges)) if len(charges) else None
    wbo = _wbo(result, len(charges))
    output['xtb_wiberg_pair_sum'] = float(np.triu(wbo, 1).sum()) if wbo is not None else None
    return output


def CalcXTBAtomDescriptors(result):
    """Return rows in calculation atom order: partial charge and total WBO.

    All off-diagonal pairs contribute to each atom's WBO sum, including pairs
    absent from an RDKit bond graph. This sum is not a formal valence.
    """
    result = _as_result(result)
    charges = _charges(result)
    if not len(charges):
        raise ValueError('atomic descriptors require partial charges')
    wbo = _wbo(result, len(charges))
    return [{'atom_index': i, 'partial_charge_e': float(q),
             'wiberg_sum': float(wbo[i].sum()) if wbo is not None else None}
            for i, q in enumerate(charges)]


def CalcXTBBondDescriptors(mol, result):
    """Return existing RDKit bond rows after checking topology/geometry provenance.

    Accepts ctypes and compiled Python results carrying input snapshot
    provenance. Missing provenance or changed topology/geometry is rejected.
    """
    from .mapping import _validate_mapping
    result = _as_result(result)
    _validate_mapping(mol, result)
    wbo = _wbo(result, mol.GetNumAtoms())
    if wbo is None:
        raise ValueError('bond descriptors require a Wiberg matrix')
    return [{'bond_index': b.GetIdx(), 'begin_atom_index': b.GetBeginAtomIdx(),
             'end_atom_index': b.GetEndAtomIdx(),
             'wiberg_order': float(wbo[b.GetBeginAtomIdx(), b.GetEndAtomIdx()])}
            for b in mol.GetBonds()]


def CalcXTBFragmentDescriptors(result, fragments):
    """Aggregate disjoint, explicitly indexed fragments; a full partition is optional.

    Internal WBO counts each pair once. External WBO includes all atoms outside
    the fragment, including unassigned atoms. Pair WBOs are reported once per
    fragment pair. Hydrogens are included only if their indices are supplied.
    """
    result = _as_result(result)
    charges = _charges(result)
    if not len(charges):
        raise ValueError('fragment descriptors require partial charges')
    groups = _fragments(fragments, len(charges))
    wbo = _wbo(result, len(charges))
    records, pairs = {}, []
    for name, indices in groups.items():
        outside = np.setdiff1d(np.arange(len(charges)), indices)
        records[name] = {
            'atom_indices': indices.tolist(), 'num_atoms': len(indices),
            'charge_e': float(charges[indices].sum()),
            'charge_min_e': float(charges[indices].min()),
            'charge_max_e': float(charges[indices].max()),
            'charge_std_e': float(charges[indices].std()),
            'internal_wiberg_sum': float(np.triu(wbo[np.ix_(indices, indices)], 1).sum()) if wbo is not None else None,
            'external_wiberg_sum': float(wbo[np.ix_(indices, outside)].sum()) if wbo is not None else None,
        }
    names = list(groups)
    for i, first in enumerate(names):
        for second in names[i+1:]:
            pairs.append({'fragment_a': first, 'fragment_b': second,
                          'wiberg_sum': float(wbo[np.ix_(groups[first], groups[second])].sum()) if wbo is not None else None})
    assigned = set(np.concatenate(list(groups.values())).tolist())
    return {'fragments': records, 'pairs': pairs,
            'unassigned_atom_indices': sorted(set(range(len(charges))) - assigned)}
