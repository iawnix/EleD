"""Electronic descriptors derived from a single result, without engine calls.

Indices always refer to the original calculation's atom/MO order. Wiberg
sums include all off-diagonal atom pairs, not just RDKit topology bonds.
"""
from __future__ import annotations

from collections.abc import Mapping
import json
from numbers import Integral

import numpy as np

from .charges import CHARGE_MODEL_SCC, _calc_cm5_charges, _charge_model_for_method
from .errors import DescriptorUnavailableError
from .result import XTBResults, _aufbau_occupied_count

AU_DIPOLE_TO_DEBYE = 2.541746473


def _as_result(result) -> XTBResults:
    """Copy numeric fields from a compiled result; never fabricate provenance."""
    if isinstance(result, XTBResults):
        if not result.success:
            detail = result.metadata.get('error_message', 'calculation failed')
            raise ValueError(f'calculation failed: {detail}')
        return result
    if hasattr(result, 'GetNumAtoms') and hasattr(result, 'GetConformers'):
        raise TypeError(
            'expected XTBResults or rdXTBTools.XTBResults, got an RDKit Mol; '
            'pass a calculation result or call CalcMolXTBDescriptors(mol)'
        )
    if (type(result).__module__ != 'eled.rdXTBTools'
            or type(result).__name__ != 'XTBResults'):
        raise TypeError('expected XTBResults or rdXTBTools.XTBResults')
    if not result.success:
        raise ValueError(f'calculation failed: {result.errorMessage}')
    basis = {}
    if result.overlapRetained:
        from .basis import _native_p_directions
        ao_shells = tuple(map(int, result.GetAOShellMap()))
        shell_l = tuple(map(int, result.GetShellAngularMomenta()))
        basis = dict(
            overlapMatrix=tuple(map(tuple, result.GetOverlapMatrix())),
            orbitalCoefficients=tuple(map(tuple, result.GetOrbitalCoefficients())),
            aoAtomMap=tuple(map(int, result.GetAOAtomMap())),
            aoShellMap=ao_shells,
            shellAtomMap=tuple(map(int, result.GetShellAtomMap())),
            shellAngularMomenta=shell_l,
            aoPDirections=_native_p_directions(ao_shells, shell_l),
        )
    metadata = (
        json.loads(result.inputProvenanceJSON)
        if result.inputProvenanceJSON
        else {'atom_count': result.numAtoms}
    )
    metadata = dict(metadata)
    charge_model = getattr(result, 'chargeModel', None) or metadata.get(
        'charge_model', _charge_model_for_method(result.method)
    )
    from .scf import scf_record
    metadata['scf'] = scf_record(
        result.method,
        getattr(result, 'scfIterations', -1),
        [
            getattr(result, name, -1.0)
            for name in (
                'scfEnergyChangeHartree',
                'scfMixerRMS',
                'scfEnergyThresholdHartree',
                'scfMixerThreshold',
            )
        ],
    )
    metadata.update(
        calculation_status='success',
        initialization=result.initialization,
        calculation_settings={
            'max_iterations': result.maxIterations,
            'accuracy': result.accuracy,
            'electronic_temperature': result.electronicTemperature,
            'effective_electronic_temperature': result.effectiveElectronicTemperature,
        },
        engine={
            'backend': 'compiled',
            'path': result.libraryPath,
            'sha256': result.librarySHA256,
            'source_commit': result.sourceCommit,
            'api_version': result.apiVersion,
        },
        parameter_file=(
            {'path': result.parameterFile, 'sha256': result.parameterSHA256}
            if result.parameterFile
            else None
        ),
    )
    metadata['charge_model'] = charge_model
    grad_rows = ()
    if hasattr(result, 'GetGradient'):
        grad = np.asarray(result.GetGradient(), dtype=float)
        grad_rows = tuple(map(tuple, grad))
        metadata['gradient_hartree_per_bohr'] = tuple(float(value) for value in grad.ravel())
    atomic_charges = tuple(map(float, result.GetAtomicCharges()))
    cm5_charges = ()
    if charge_model == 'xtb_scc_mulliken' and hasattr(result, 'GetCM5Charges'):
        cm5_charges = tuple(map(float, result.GetCM5Charges()))
    elif charge_model == 'xtb_scc_mulliken':
        coordinates = metadata.get('coordinates_angstrom')
        atomic_numbers = metadata.get('atomic_numbers', ())
        if coordinates is not None and atomic_numbers is not None and len(atomic_numbers):
            cm5_charges = _calc_cm5_charges(atomic_numbers, coordinates, atomic_charges)
    esp = {}
    if getattr(result, 'espSampled', False):
        esp = dict(
            espPointsAngstrom=tuple(map(tuple, result.GetESPPoints())),
            espValuesAU=tuple(map(float, result.GetESPValues())),
        )
        metadata['esp'] = {
            'api_version': 1,
            'potential_unit': 'hartree/e',
            'coordinate_unit': 'angstrom',
            'source': 'native_density',
        }
    realspace = {}
    if getattr(result, 'realSpaceSampled', False):
        realspace = dict(
            realSpacePointsAngstrom=tuple(map(tuple, result.GetRealSpacePoints())),
            realSpaceDensityAU=tuple(map(float, result.GetRealSpaceDensity())),
            realSpaceALIEEV=tuple(map(float, result.GetRealSpaceALIE())),
            realSpaceVirtualDensityAU=tuple(map(float, result.GetRealSpaceVirtualDensity())),
            realSpaceLEAEV=tuple(map(float, result.GetRealSpaceLEA())),
        )
        metadata['realspace'] = {
            'api_version': 1,
            'density_unit': 'e/bohr^3',
            'alie_unit': 'eV',
            'virtual_density_unit': '1/bohr^3',
            'lea_unit': 'eV',
            'coordinate_unit': 'angstrom',
            'source': 'native_cartesian_gto',
        }
    return XTBResults(
        **esp,
        **realspace,
        **basis,
        gradientHartreePerBohr=grad_rows,
        method=result.method,
        version=None,
        totalEnergyHartree=result.totalEnergyHartree,
        homoLumoGapEV=None,
        numOrbitals=result.numOrbitals,
        numElectrons=result.numElectrons,
        numUnpairedElectrons=result.numUnpairedElectrons,
        orbitalEnergiesEV=tuple(map(float, result.GetOrbitalEnergies())),
        orbitalOccupancies=tuple(map(float, result.GetOrbitalOccupancies())),
        atomicCharges=atomic_charges,
        chargeModel=charge_model,
        cm5Charges=cm5_charges,
        dipoleAtomicUnits=tuple(map(float, result.GetDipole())),
        wibergBondOrders=tuple(map(tuple, result.GetWibergBondOrders())),
        backend='compiled',
        metadata=metadata,
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


def _orbitals(result: XTBResults):
    """Return energies, occupations and the aufbau occupied/virtual partition.

    ``occupied`` holds the lowest ``n_alpha`` spatial orbitals and ``virtual``
    the rest; occupations are not thresholded (see
    :func:`eled.result._aufbau_occupied_count`).
    """
    energies = _finite(result.orbitalEnergiesEV, 'orbital energies', 1)
    occupations = _finite(result.orbitalOccupancies, 'occupations', 1)
    if energies.shape != occupations.shape or np.any(occupations < 0):
        raise ValueError('orbital energies and nonnegative occupations must have matching lengths')
    if np.any(np.diff(energies) < -1e-10):
        raise ValueError('orbital energies must be sorted')
    if result.numOrbitals is not None and len(energies) != result.numOrbitals:
        raise ValueError('orbital count does not match arrays')
    n_occ = _aufbau_occupied_count(
        len(energies), occupations, result.numElectrons, result.numUnpairedElectrons
    )
    order = np.arange(len(energies))
    return energies, occupations, order[:n_occ], order[n_occ:]


def _charges(result: XTBResults):
    charges = _finite(result.atomicCharges, 'partial charges', 1)
    expected = result.metadata.get('atom_count')
    if len(charges):
        if expected is not None and len(charges) != expected:
            raise ValueError(f'partial charges do not match atom count: expected {expected}, got {len(charges)}')
        if result.method in ('GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'):
            if result.chargeModel != _charge_model_for_method(result.method):
                raise ValueError(f'{result.chargeModel} is inconsistent with {result.method}')
    return charges


def _cm5_charges(result: XTBResults):
    if result.cm5Charges is None:
        return np.empty(0)
    values = _finite(result.cm5Charges, 'CM5 charges', 1)
    if len(values):
        if result.method not in ('GFN1-xTB', 'GFN2-xTB') or result.chargeModel != CHARGE_MODEL_SCC:
            raise ValueError(f'CM5 charges are unavailable for {result.method} with {result.chargeModel}')
        for expected in (result.metadata.get('atom_count'), len(result.atomicCharges) or None):
            if expected is not None and len(values) != expected:
                raise ValueError(f'CM5 charges do not match atom count: expected {expected}, got {len(values)}')
    return values


def CalcCM5Charges(result):
    """Return one derived CM5 charge per atom in elementary-charge units."""
    values = _cm5_charges(_as_result(result))
    if not len(values):
        raise DescriptorUnavailableError('CM5 charges are unavailable; GFN1/2 SCC charges and geometry are required')
    return tuple(map(float, values))


def _has_population_basis(result: XTBResults):
    fields = (result.overlapMatrix, result.orbitalCoefficients, result.aoAtomMap)
    return all(value is not None and len(value) > 0 for value in fields)


def _orbital_basis_arrays(result: XTBResults, nOrbitals, nAtoms):
    coefficients = _finite(result.orbitalCoefficients, 'orbital coefficients', 2)
    overlap = _finite(result.overlapMatrix, 'overlap', 2)
    shape = (nOrbitals, nOrbitals)
    if not nOrbitals or coefficients.shape != shape or overlap.shape != shape:
        raise ValueError('coefficients and overlap must match the full square AO/MO basis')
    if not np.allclose(overlap, overlap.T, atol=1e-10, rtol=0):
        raise ValueError('overlap must be symmetric')
    mapping = np.asarray(result.aoAtomMap)
    if (mapping.shape != (nOrbitals,) or mapping.dtype.kind not in 'iu'
            or not np.array_equal(np.unique(mapping), np.arange(nAtoms))):
        raise ValueError('AO-to-atom mapping must cover every atom with valid integer indices')
    return coefficients, overlap, mapping


def _wbo(result: XTBResults, natoms):
    if not len(result.wibergBondOrders):
        return None
    wbo = _finite(result.wibergBondOrders, 'Wiberg matrix', 2)
    if wbo.shape != (natoms, natoms) or not np.allclose(wbo, wbo.T, atol=1e-8, rtol=0):
        raise ValueError('Wiberg matrix must be symmetric and match atom count')
    wbo = wbo.copy()
    np.fill_diagonal(wbo, 0.0)
    return wbo


def _indices(values, limit, name):
    values = list(values)
    if (not values or any(isinstance(v, (bool, np.bool_)) or not isinstance(v, Integral)
                          or v < 0 or v >= limit for v in values)
            or len(set(values)) != len(values)):
        raise ValueError(f'{name} must contain unique integer indices in [0, {limit})')
    return np.asarray(values, dtype=int)


def _fragments(fragments, natoms):
    if hasattr(fragments, 'atoms'):
        fragments = fragments.atoms
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


def CalcXTBDescriptors(result, *, frontierDepth=2):
    """Return JSON-compatible molecular scalars; unavailable values are None.

    HOMO/LUMO follow the aufbau electron count, not an occupation threshold.
    ``frontierDepth=2`` adds HOMO-1/-2 and LUMO+1/+2. Charge standard deviation
    uses ``ddof=0``. Dipoles refer to the supplied coordinate origin (relevant
    for ions); a.u. denotes ``e*Bohr``.
    """
    if isinstance(frontierDepth, bool) or not isinstance(frontierDepth, Integral) or frontierDepth < 0:
        raise ValueError('frontierDepth must be a nonnegative integer')
    result = _as_result(result)
    energies, occupations, occupied, virtual = _orbitals(result)
    charges = _charges(result)
    dipole = _finite(result.dipoleAtomicUnits, 'dipole', 1)
    if len(dipole) not in (0, 3):
        raise ValueError('dipole must have three components')
    energy = result.totalEnergyHartree
    if energy is not None and not np.isfinite(energy):
        raise ValueError('total energy must be finite')
    homo = float(energies[occupied[-1]]) if len(occupied) else None
    lumo = float(energies[virtual[0]]) if len(virtual) else None
    output = {
        'xtb_method': result.method,
        'xtb_backend': result.backend,
        'xtb_version': result.version,
        'xtb_energy_ha': energy,
        'xtb_homo_ev': homo,
        'xtb_lumo_ev': lumo,
        'xtb_gap_ev': lumo - homo if homo is not None and lumo is not None else None,
        'xtb_n_orb': result.numOrbitals,
        'xtb_n_e': result.numElectrons,
        'xtb_n_unpaired': result.numUnpairedElectrons,
        'xtb_homo_idx': int(occupied[-1]) if len(occupied) else None,
        'xtb_lumo_idx': int(virtual[0]) if len(virtual) else None,
    }
    for offset in range(1, frontierDepth + 1):
        output[f'xtb_homo_m{offset}_ev'] = float(energies[occupied[-1 - offset]]) if len(occupied) > offset else None
        output[f'xtb_lumo_p{offset}_ev'] = float(energies[virtual[offset]]) if len(virtual) > offset else None
    for i, axis in enumerate('xyz'):
        output[f'xtb_dip_{axis}_ebohr'] = float(dipole[i]) if len(dipole) else None
    norm = float(np.linalg.norm(dipole)) if len(dipole) else None
    output['xtb_dip_ebohr'] = norm
    output['xtb_dip_debye'] = norm * AU_DIPOLE_TO_DEBYE if norm is not None else None
    for label, func in [
        ('sum', np.sum),
        ('min', np.min),
        ('max', np.max),
        ('mean', np.mean),
        ('sd', np.std),
        ('abs', lambda q: np.abs(q).sum()),
    ]:
        output[f'xtb_q_{label}_e'] = float(func(charges)) if len(charges) else None
    cm5 = _cm5_charges(result)
    if len(cm5):
        for label, func in [
            ('sum', np.sum),
            ('min', np.min),
            ('max', np.max),
            ('mean', np.mean),
            ('sd', np.std),
            ('abs', lambda q: np.abs(q).sum()),
        ]:
            output[f'xtb_cm5_q_{label}_e'] = float(func(cm5))
    wbo = _wbo(result, len(charges))
    output['xtb_wbo_sum'] = float(np.triu(wbo, 1).sum()) if wbo is not None else None
    return output


def CalcXTBAtomDescriptors(result):
    """Return rows in calculation atom order: partial charge and total WBO."""
    result = _as_result(result)
    charges = _charges(result)
    if not len(charges):
        raise ValueError('atomic descriptors require partial charges')
    wbo = _wbo(result, len(charges))
    cm5 = _cm5_charges(result)
    rows = []
    for i, q in enumerate(charges):
        row = {
            'atom_index': i,
            'q_e': float(q),
            'wbo_sum': float(wbo[i].sum()) if wbo is not None else None,
        }
        if len(cm5):
            row['cm5_q_e'] = float(cm5[i])
        rows.append(row)
    return rows


def CalcXTBBondDescriptors(mol, result):
    """Return existing RDKit bond rows after checking topology/geometry provenance."""
    from .mapping import _validate_mapping
    result = _as_result(result)
    _validate_mapping(mol, result)
    wbo = _wbo(result, mol.GetNumAtoms())
    if wbo is None:
        raise ValueError('bond descriptors require a Wiberg matrix')
    return [
        {
            'bond_index': b.GetIdx(),
            'begin_atom_index': b.GetBeginAtomIdx(),
            'end_atom_index': b.GetEndAtomIdx(),
            'wbo': float(wbo[b.GetBeginAtomIdx(), b.GetEndAtomIdx()]),
        }
        for b in mol.GetBonds()
    ]


def CalcXTBFragmentDescriptors(result, fragments):
    """Aggregate disjoint, explicitly indexed fragments; a full partition is optional."""
    result = _as_result(result)
    charges = _charges(result)
    if not len(charges):
        raise ValueError('fragment descriptors require partial charges')
    groups = _fragments(fragments, len(charges))
    wbo = _wbo(result, len(charges))
    records, pairs = {}, []
    cm5 = _cm5_charges(result)
    for name, indices in groups.items():
        outside = np.setdiff1d(np.arange(len(charges)), indices)
        records[name] = {
            'atom_indices': indices.tolist(),
            'num_atoms': len(indices),
            'q_e': float(charges[indices].sum()),
            'q_min_e': float(charges[indices].min()),
            'q_max_e': float(charges[indices].max()),
            'q_sd_e': float(charges[indices].std()),
            'wbo_intra': float(np.triu(wbo[np.ix_(indices, indices)], 1).sum()) if wbo is not None else None,
            'wbo_inter': float(wbo[np.ix_(indices, outside)].sum()) if wbo is not None else None,
        }
        if len(cm5):
            records[name].update(
                cm5_q_e=float(cm5[indices].sum()),
                cm5_q_min_e=float(cm5[indices].min()),
                cm5_q_max_e=float(cm5[indices].max()),
                cm5_q_sd_e=float(cm5[indices].std()),
            )
    names = list(groups)
    for i, first in enumerate(names):
        for second in names[i + 1:]:
            pairs.append({
                'fragment_a': first,
                'fragment_b': second,
                'wbo': float(wbo[np.ix_(groups[first], groups[second])].sum()) if wbo is not None else None,
            })
    assigned = set(np.concatenate(list(groups.values())).tolist())
    return {
        'fragments': records,
        'pairs': pairs,
        'unassigned_atom_indices': sorted(set(range(len(charges))) - assigned),
    }
