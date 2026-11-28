"""Vertical charge-state energy differences and condensed Fukui descriptors.

No geometry relaxation or spin-state search. GFN1/2 only: GFN0 is not an SCF
method. Finite-temperature energies and unbound attachment are labeled, not
interpreted as guaranteed gas-phase experimental IP/EA values.
"""
from __future__ import annotations

from contextlib import closing
from copy import deepcopy
from dataclasses import replace
import json

import numpy as np

from .backends.native import _native_charge_states
from .constants import HARTREE_TO_EV
from .calculator import XTBCalculator
from .charges import CHARGE_MODEL_CM5, CHARGE_MODEL_SCC, _charge_model_for_method
from .electronic_descriptors import CalcCM5Charges, _as_result, _finite, _fragments
from .errors import XTBEnvironmentError, XTBRunError
from .input import _effective_charge_uhf, _xyz_from_mol
from .mapping import _validate_mapping
from .options import XTBOptions
from .populations import CalcOrbitalPopulations
from .reports import DESCRIPTOR_VERSION, SCHEMA_VERSION, AssessDescriptorQuality
from .selections import NormalizeSelections
from .xtb import _resolve_options

REACTIVITY_UNITS = {
    'ip_v_ev': 'eV',
    'ea_v_ev': 'eV',
    'mu_ev': 'eV',
    'chi_ev': 'eV',
    'eta_ev': 'eV',
    's_1ev': '1/eV',
    'omega_ev': 'eV',
    'omega_minus_ev': 'eV',
    'omega_plus_ev': 'eV',
}

REACTIVITY_PROTOCOL_VERSION = '2.2'


class ChargeStateCalculationError(XTBRunError):
    """A required charge state failed; ``states`` contains completed diagnostics."""

    def __init__(self, state, cause, states):
        super().__init__(f'{state} charge-state calculation failed: {cause}')
        self.state = state
        self.states = deepcopy(states)


def _summary(result):
    return {
        'status': 'success',
        'method': result.method,
        'charge': result.metadata['charge'],
        'uhf': result.metadata['uhf'],
        'energy_hartree': result.totalEnergyHartree,
        'num_electrons': result.numElectrons,
        'partial_charges_e': list(result.atomicCharges),
        'charge_model': result.chargeModel,
        'cm5_charges_e': list(result.cm5Charges),
        'provenance': deepcopy(dict(result.metadata)),
    }


def _empty_report(recordId, scheme, states):
    if recordId is not None and not isinstance(recordId, str):
        raise ValueError('recordId must be a string or None')
    scheme = _scheme(scheme)
    charge_model = CHARGE_MODEL_CM5 if scheme == 'cm5' else CHARGE_MODEL_SCC
    population_model = _population_model(scheme)
    return {
        'schema_version': SCHEMA_VERSION,
        'descriptor_version': DESCRIPTOR_VERSION,
        'report_type': 'vertical_reactivity',
        'record_id': recordId,
        'values': dict.fromkeys(REACTIVITY_UNITS),
        'units': dict(REACTIVITY_UNITS),
        'analysis_settings': {
            'geometry': 'fixed_reference_geometry',
            'population_scheme': scheme,
            'charge_model': charge_model,
            'population_model': population_model,
            'energy_definition': 'libxtb_total_energy',
            'hardness': '(IP-EA)/2',
            'softness': '1/(IP-EA)',
            'electrophilicity': 'mu^2/(2*hardness)',
            'electrodonating_power': '(3*IP+EA)^2/(16*(IP-EA))',
            'electroaccepting_power': '(IP+3*EA)^2/(16*(IP-EA))',
            'f_plus': 'population(N+1)-population(N)',
            'f_minus': 'population(N)-population(N-1)',
            'dual': 'f_plus-f_minus',
            'omega_minus_loc_ev': 'omega_minus_ev*f_minus',
            'omega_plus_loc_ev': 'omega_plus_ev*f_plus',
            'delta_omega_loc_ev': 'omega_plus_ev*f_plus-omega_minus_ev*f_minus',
            's2_loc_1ev2': 's_1ev^2*dual',
            'spin_search': False,
        },
        'states': states,
        'provenance': {'states': states},
        'quality': {},
        'missing': {},
        'error': None,
        'atoms': None,
        'fragments': None,
        'atom_units': {
            'f_plus': 'dimensionless',
            'f_minus': 'dimensionless',
            'f_radical': 'dimensionless',
            'dual': 'dimensionless',
            's_loc_1ev': '1/eV',
            'omega_loc_ev': 'eV',
            'n_loc_1ev': '1/eV',
            'omega_minus_loc_ev': 'eV',
            'omega_plus_loc_ev': 'eV',
            'delta_omega_loc_ev': 'eV',
            's2_loc_1ev2': '1/eV^2',
        },
        'atom_invariance': 'atom index and molecular geometry fixed; values invariant to rigid translation/rotation',
        'atom_missing_semantics': (
            'None means the required global softness/electrophilicity is undefined '
            'or the state failed; never zero-filled'
        ),
    }


def _scheme(scheme):
    if scheme not in ('xtb_scc', 'cm5', 'mulliken', 'lowdin'):
        raise ValueError("populationScheme must be 'xtb_scc', 'cm5', 'mulliken' or 'lowdin'")
    return scheme


def _population_model(scheme):
    return {
        'xtb_scc': 'atomic_charge_difference',
        'cm5': 'cm5_atomic_charge_difference',
        'mulliken': 'mulliken_orbital_population',
        'lowdin': 'lowdin_orbital_population',
    }[scheme]


def _result_charge_model(result):
    model = getattr(result, 'chargeModel', None)
    if model and model != 'unknown':
        return model
    return _charge_model_for_method(result.method)


def CalcReactivityDescriptors(
    mol,
    reference,
    ionized,
    attached,
    *,
    fragments=None,
    reactionSites=None,
    populationScheme='xtb_scc',
    recordId=None,
):
    """Analyze three compatible successful results, without further engine calls."""
    populationScheme = _scheme(populationScheme)
    results = {
        k: _as_result(v)
        for k, v in [('reference', reference), ('ionized', ionized), ('attached', attached)]
    }
    ref = results['reference']
    if ref.method not in ('GFN1-xTB', 'GFN2-xTB'):
        raise ValueError('vertical Delta-SCF requires self-consistent GFN1/2, not GFN0')
    base = ref.metadata
    for label, result in results.items():
        if any(result.metadata.get(key) is None for key in ('charge', 'uhf')):
            raise ValueError(f'{label} lacks charge/spin provenance')
    if not base.get('engine', {}).get('sha256') or not base.get('calculation_settings'):
        raise ValueError('charge-state analysis requires engine identity and calculation settings')

    def identity(result):
        meta = result.metadata
        return (
            result.method,
            meta.get('engine', {}).get('sha256'),
            (meta.get('parameter_file') or {}).get('sha256'),
            meta.get('calculation_settings'),
            meta.get('solvent'),
        )

    states = {}
    flags = []
    for label, result in results.items():
        _validate_mapping(mol, result)
        if result.metadata['geometry_sha256'] != base['geometry_sha256']:
            raise ValueError('charge states must use the same fixed geometry')
        if result.metadata.get('calculation_status') != 'success':
            raise ValueError('charge states must record successful calculations')
        if identity(result) != identity(ref):
            raise ValueError('charge states must share method, engine, parameters and numerical settings')
        if result.numOrbitals != ref.numOrbitals:
            raise ValueError('charge states must use the same orbital basis size')
        delta = {'reference': 0, 'ionized': 1, 'attached': -1}[label]
        if result.metadata['charge'] != base['charge'] + delta:
            raise ValueError(f'{label} has the wrong molecular charge')
        if result.numElectrons != ref.numElectrons - delta:
            raise ValueError(f'{label} has the wrong modeled electron count')
        if result.numUnpairedElectrons != result.metadata['uhf']:
            raise ValueError('charge-state spin count disagrees with its provenance')
        energy = result.totalEnergyHartree
        if energy is None or not np.isfinite(energy):
            raise ValueError('all three energies must be finite')
        charges = _finite(result.atomicCharges, 'charges', 1)
        if _result_charge_model(result) != CHARGE_MODEL_SCC:
            raise ValueError('vertical GFN1/2 charge states must expose SCC Mulliken charges')
        if len(charges) != mol.GetNumAtoms():
            raise ValueError('charges must match the indexed atoms')
        quality = AssessDescriptorQuality(result)
        if quality['status'] == 'error':
            raise ValueError(f'{label} failed charge/electron conservation checks')
        flags.extend(dict(flag, state=label) for flag in quality['flags'])
        states[label] = _summary(result)
    groups = _fragments(fragments, mol.GetNumAtoms()) if fragments is not None else {}
    if reactionSites is None:
        sites = {}
    elif hasattr(reactionSites, 'reactionSites'):
        sites = NormalizeSelections(mol, reactionSites).reactionSites
    else:
        sites = NormalizeSelections(mol, {'reactionSites': reactionSites}).reactionSites
    report = _empty_report(recordId, populationScheme, states)
    ip = (results['ionized'].totalEnergyHartree - ref.totalEnergyHartree) * HARTREE_TO_EV
    ea = (ref.totalEnergyHartree - results['attached'].totalEnergyHartree) * HARTREE_TO_EV
    mu = -(ip + ea) / 2
    eta = (ip - ea) / 2
    values = report['values']
    values.update(ip_v_ev=ip, ea_v_ev=ea, mu_ev=mu, chi_ev=-mu, eta_ev=eta)
    if eta > 1e-8:
        values.update(
            s_1ev=1 / (2 * eta),
            omega_ev=mu * mu / (2 * eta),
            omega_minus_ev=((3 * ip + ea) ** 2) / (32 * eta),
            omega_plus_ev=((ip + 3 * ea) ** 2) / (32 * eta),
        )
    else:
        flags.append({
            'code': 'nonpositive_or_tiny_hardness',
            'severity': 'warning',
            'detail': 'Reciprocal hardness descriptors are undefined or numerically unstable.',
        })
        report['missing'].update(
            s_1ev='nonpositive_or_tiny_hardness',
            omega_ev='nonpositive_or_tiny_hardness',
            omega_minus_ev='nonpositive_or_tiny_hardness',
            omega_plus_ev='nonpositive_or_tiny_hardness',
        )
    if ea < 0:
        flags.append({
            'code': 'negative_vertical_ea',
            'severity': 'warning',
            'detail': 'Attachment is energetically unfavorable at this geometry; no bound-anion claim.',
        })
    if ip < 0:
        flags.append({
            'code': 'negative_vertical_ip',
            'severity': 'warning',
            'detail': 'Ionization lowers the calculated energy; inspect charge/spin states.',
        })
    if base['calculation_settings']['effective_electronic_temperature'] > 1e-6:
        flags.append({
            'code': 'finite_electronic_temperature',
            'severity': 'warning',
            'detail': 'Energy differences use the configured finite-temperature libxtb energies.',
        })
    if populationScheme == 'xtb_scc':
        q0, qp, qm = [np.asarray(results[k].atomicCharges) for k in ('reference', 'ionized', 'attached')]
        plus, minus = q0 - qm, qp - q0
    elif populationScheme == 'cm5':
        q0, qp, qm = [np.asarray(CalcCM5Charges(results[k])) for k in ('reference', 'ionized', 'attached')]
        plus, minus = q0 - qm, qp - q0
    else:
        n0, np_, nm = [
            np.asarray(CalcOrbitalPopulations(results[k], scheme=populationScheme)['mo_n_e'])
            for k in ('reference', 'ionized', 'attached')
        ]
        plus, minus = nm - n0, n0 - np_
    dual, radical = plus - minus, (plus + minus) / 2
    errors = {'f_plus': float(plus.sum() - 1), 'f_minus': float(minus.sum() - 1), 'dual': float(dual.sum())}
    if max(abs(v) for v in errors.values()) > 1e-5:
        raise ValueError('condensed Fukui sums violate electron conservation')
    softness = values.get('s_1ev')
    electrophilicity = values.get('omega_ev')
    omega_minus = values.get('omega_minus_ev')
    omega_plus = values.get('omega_plus_ev')
    softness_sq = softness * softness if softness is not None else None
    atoms = []
    for i in range(len(plus)):
        om_loc = float(omega_minus * minus[i]) if omega_minus is not None else None
        op_loc = float(omega_plus * plus[i]) if omega_plus is not None else None
        atoms.append({
            'atom_index': i,
            'f_plus': float(plus[i]),
            'f_minus': float(minus[i]),
            'f_radical': float(radical[i]),
            'dual': float(dual[i]),
            's_loc_1ev': float(softness * plus[i]) if softness is not None else None,
            'omega_loc_ev': float(electrophilicity * plus[i]) if electrophilicity is not None else None,
            'n_loc_1ev': float(softness * minus[i]) if softness is not None else None,
            'omega_minus_loc_ev': om_loc,
            'omega_plus_loc_ev': op_loc,
            'delta_omega_loc_ev': float(op_loc - om_loc) if (op_loc is not None and om_loc is not None) else None,
            's2_loc_1ev2': float(softness_sq * dual[i]) if softness_sq is not None else None,
        })
    report['atoms'] = atoms
    report['fragments'] = {
        name: {
            'atom_indices': indices.tolist(),
            'f_plus': float(plus[indices].sum()),
            'f_minus': float(minus[indices].sum()),
            'f_radical': float(radical[indices].sum()),
            'dual': float(dual[indices].sum()),
            's_loc_1ev': float(softness * plus[indices].sum()) if softness is not None else None,
            'omega_loc_ev': float(electrophilicity * plus[indices].sum()) if electrophilicity is not None else None,
            'n_loc_1ev': float(softness * minus[indices].sum()) if softness is not None else None,
            'omega_minus_loc_ev': float(omega_minus * minus[indices].sum()) if omega_minus is not None else None,
            'omega_plus_loc_ev': float(omega_plus * plus[indices].sum()) if omega_plus is not None else None,
            'delta_omega_loc_ev': (
                float(omega_plus * plus[indices].sum() - omega_minus * minus[indices].sum())
                if (omega_plus is not None and omega_minus is not None)
                else None
            ),
            's2_loc_1ev2': float(softness_sq * dual[indices].sum()) if softness_sq is not None else None,
        }
        for name, indices in groups.items()
    }
    report['reaction_sites'] = {
        name: {
            'atom_indices': list(indices),
            'f_plus': float(plus[list(indices)].sum()),
            'f_minus': float(minus[list(indices)].sum()),
            'f_radical': float(radical[list(indices)].sum()),
            'dual': float(dual[list(indices)].sum()),
            's_loc_1ev': float(softness * plus[list(indices)].sum()) if softness is not None else None,
            'omega_loc_ev': float(electrophilicity * plus[list(indices)].sum()) if electrophilicity is not None else None,
            'n_loc_1ev': float(softness * minus[list(indices)].sum()) if softness is not None else None,
            'omega_minus_loc_ev': float(omega_minus * minus[list(indices)].sum()) if omega_minus is not None else None,
            'omega_plus_loc_ev': float(omega_plus * plus[list(indices)].sum()) if omega_plus is not None else None,
            'delta_omega_loc_ev': (
                float(omega_plus * plus[list(indices)].sum() - omega_minus * minus[list(indices)].sum())
                if (omega_plus is not None and omega_minus is not None)
                else None
            ),
            's2_loc_1ev2': float(softness_sq * dual[list(indices)].sum()) if softness_sq is not None else None,
        }
        for name, indices in sites.items()
    }
    report['unassigned_atom_indices'] = sorted(
        set(range(len(plus))) - set(i for group in groups.values() for i in group)
    )
    report['quality'] = {
        'status': 'warning' if any(f['severity'] == 'warning' for f in flags) else 'ok',
        'flags': flags,
        'fukui_sum_errors': errors,
        'spin_states_verified_as_ground_states': False,
    }
    json.dumps(report, allow_nan=False)
    return report


def RunReactivity(
    mol,
    *,
    options=None,
    ionizedUnpaired=None,
    attachedUnpaired=None,
    fragments=None,
    reactionSites=None,
    populationScheme='xtb_scc',
    interface='native',
    recordId=None,
    returnErrors=False,
    initialization='reference',
):
    """Calculate fixed-geometry Q, Q+1, Q-1 states with explicit spin policy."""
    from rdkit import Chem
    populationScheme = _scheme(populationScheme)
    if interface not in ('native', 'compiled'):
        raise ValueError('interface must be native or compiled')
    if initialization not in ('reference', 'independent'):
        raise ValueError('initialization must be reference or independent')
    if not isinstance(returnErrors, bool):
        raise ValueError('returnErrors must be Boolean')
    if recordId is not None and not isinstance(recordId, str):
        raise ValueError('recordId must be a string or None')
    snapshot = Chem.Mol(mol)
    config = options or XTBOptions()
    _xyz_from_mol(snapshot, config.confId)
    charge, unpaired = _effective_charge_uhf(snapshot, config.charge, config.numUnpairedElectrons)
    if config.method not in ('GFN1-xTB', 'GFN2-xTB'):
        raise ValueError('vertical Delta-SCF requires self-consistent GFN1/2, not GFN0')
    if unpaired and (ionizedUnpaired is None or attachedUnpaired is None):
        raise ValueError(
            'open-shell references require explicit ionizedUnpaired and attachedUnpaired; no automatic spin search'
        )
    spins = {
        'reference': unpaired,
        'ionized': 1 if ionizedUnpaired is None else ionizedUnpaired,
        'attached': 1 if attachedUnpaired is None else attachedUnpaired,
    }
    charges = {'reference': charge, 'ionized': charge + 1, 'attached': charge - 1}
    for label in spins:
        _, spins[label] = _effective_charge_uhf(snapshot, charges[label], spins[label])
    if fragments is not None:
        _fragments(fragments, snapshot.GetNumAtoms())
    if reactionSites is not None:
        if hasattr(reactionSites, 'reactionSites'):
            NormalizeSelections(snapshot, reactionSites)
        else:
            NormalizeSelections(snapshot, {'reactionSites': reactionSites})
    config = _resolve_options(
        replace(
            config,
            charge=charge,
            numUnpairedElectrons=unpaired,
            keepOverlap=config.keepOverlap or populationScheme in ('mulliken', 'lowdin'),
        )
    )
    protocol = {
        'requested_options': config.ToDict(),
        'interface': interface,
        'protocol_version': REACTIVITY_PROTOCOL_VERSION,
        'initialization': initialization,
        'charge_model': CHARGE_MODEL_CM5 if populationScheme == 'cm5' else CHARGE_MODEL_SCC,
        'population_model': _population_model(populationScheme),
        'charge_states': charges,
        'unpaired_electrons': spins,
        'spin_selection': (
            'closed_shell_doublet_defaults'
            if not unpaired and ionizedUnpaired is None and attachedUnpaired is None
            else 'explicit_adjacent_spins'
        ),
    }
    configurations = [
        replace(config, charge=charges[k], numUnpairedElectrons=spins[k])
        for k in ('reference', 'ionized', 'attached')
    ]
    results, states = {}, {}
    calculations = _charge_state_results(snapshot, configurations, interface, initialization)
    with closing(calculations):
        for label in ('reference', 'ionized', 'attached'):
            try:
                result = next(calculations)
                results[label] = result
                states[label] = _summary(result)
            except XTBEnvironmentError:
                raise
            except XTBRunError as exc:
                states[label] = {
                    'status': 'error',
                    'charge': charges[label],
                    'uhf': spins[label],
                    'error': {'type': type(exc).__name__, 'message': str(exc)},
                }
                for remaining in charges:
                    if remaining not in states:
                        states[remaining] = {
                            'status': 'not_run',
                            'charge': charges[remaining],
                            'uhf': spins[remaining],
                        }
                error = ChargeStateCalculationError(label, exc, states)
                if not returnErrors:
                    raise error from exc
                report = _empty_report(recordId, populationScheme, states)
                report['analysis_settings'].update(protocol)
                report['error'] = {'type': type(error).__name__, 'message': str(error), 'state': label}
                report['quality'] = {
                    'status': 'error',
                    'flags': [{'code': 'charge_state_failed', 'severity': 'error', 'detail': str(error)}],
                }
                report['missing'] = dict.fromkeys(REACTIVITY_UNITS, 'required_charge_state_failed')
                return report
    report = CalcReactivityDescriptors(
        snapshot,
        **results,
        fragments=fragments,
        reactionSites=reactionSites,
        populationScheme=populationScheme,
        recordId=recordId,
    )
    report['analysis_settings'].update(protocol)
    return report


def _charge_state_results(mol, configurations, interface, initialization):
    if interface == 'compiled':
        from . import rdXTBTools
        raw_results = rdXTBTools._RunChargeStates(
            mol,
            configurations[0],
            configurations[1].numUnpairedElectrons,
            configurations[2].numUnpairedElectrons,
            initialization == 'reference',
        )
        for raw in raw_results:
            if not raw.success:
                error_type = XTBRunError if raw.failureKind == 'calculation' else XTBEnvironmentError
                raise error_type(raw.errorMessage)
            yield _as_result(raw)
    elif initialization == 'reference':
        with closing(_native_charge_states(mol, configurations)) as calculations:
            yield from calculations
    else:
        for selected in configurations:
            yield XTBCalculator(**selected.ToDict()).Calculate(mol)
