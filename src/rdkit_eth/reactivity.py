"""Vertical charge-state energy differences and condensed Fukui descriptors.

No geometry relaxation or spin-state search. GFN1/2 only: GFN0 is not an SCF
method. Finite-temperature energies and unbound attachment are labeled, not
interpreted as guaranteed gas-phase experimental IP/EA values.
"""
from __future__ import annotations

from copy import deepcopy
from contextlib import closing
from dataclasses import replace
import json

import numpy as np

from .backends.native import HARTREE_TO_EV, _native_charge_states
from .electronic_descriptors import _as_result, _finite, _fragments
from .errors import XTBRunError
from .input import _effective_charge_uhf, _xyz_from_mol
from .mapping import _validate_mapping
from .options import XTBOptions
from .populations import calc_orbital_populations
from .reports import SCHEMA_VERSION, DESCRIPTOR_VERSION, assess_descriptor_quality
from .xtb import run_xtb, _resolve_options

REACTIVITY_UNITS = {'vertical_ip_eV':'eV', 'vertical_ea_eV':'eV',
                    'chemical_potential_eV':'eV', 'electronegativity_eV':'eV',
                    'hardness_eV':'eV', 'softness_per_eV':'1/eV', 'electrophilicity_eV':'eV'}

REACTIVITY_PROTOCOL_VERSION = '2.0'


class ChargeStateCalculationError(XTBRunError):
    """A required charge state failed; states contains completed diagnostics."""
    def __init__(self, state, cause, states):
        super().__init__(f'{state} charge-state calculation failed: {cause}')
        self.state = state
        self.states = deepcopy(states)


def _summary(result):
    return {'status':'success', 'method':result.method, 'charge':result.metadata['charge'],
            'uhf':result.metadata['uhf'], 'energy_hartree':result.total_energy_eh,
            'num_electrons':result.number_of_electrons,
            'partial_charges_e':list(result.partial_charges), 'provenance':deepcopy(dict(result.metadata))}


def _empty_report(record_id, scheme, states):
    if record_id is not None and not isinstance(record_id, str): raise ValueError('record_id must be a string or None')
    return {'schema_version':SCHEMA_VERSION, 'descriptor_version':DESCRIPTOR_VERSION,
            'report_type':'vertical_reactivity', 'record_id':record_id,
            'values':dict.fromkeys(REACTIVITY_UNITS), 'units':dict(REACTIVITY_UNITS),
            'analysis_settings':{'geometry':'fixed_reference_geometry', 'population_scheme':scheme,
                'energy_definition':'libxtb_total_energy', 'hardness':'(IP-EA)/2',
                'softness':'1/(IP-EA)', 'electrophilicity':'mu^2/(2*hardness)',
                'f_plus':'population(N+1)-population(N)',
                'f_minus':'population(N)-population(N-1)', 'dual':'f_plus-f_minus',
                'spin_search':False},
            'states':states, 'provenance':{'states':states}, 'quality':{},
            'missing':{}, 'error':None, 'atoms':None, 'fragments':None}


def _scheme(scheme):
    if scheme not in ('charges','mulliken','lowdin'):
        raise ValueError("population_scheme must be 'charges', 'mulliken' or 'lowdin'")


def calc_reactivity_descriptors(mol, reference, ionized, attached, *, fragments=None,
                                population_scheme='charges', record_id=None):
    """Analyze three compatible successful results, without further engine calls.

    ionized has Q+1 / N-1; attached has Q-1 / N+1. All use the exact same
    indexed input geometry, method, library, parameters and numerical settings.
    'charges' takes differences of the method's engine partial charges. Orbital
    schemes use occupation-weighted populations, not normalized frontier MOs.
    """
    _scheme(population_scheme)
    results = {k:_as_result(v) for k,v in [('reference',reference),('ionized',ionized),('attached',attached)]}
    ref = results['reference']
    if ref.method not in ('GFN1-xTB','GFN2-xTB'): raise ValueError('vertical Delta-SCF requires self-consistent GFN1/2, not GFN0')
    base = ref.metadata
    if not base.get('engine',{}).get('sha256') or not base.get('calculation_settings'):
        raise ValueError('charge-state analysis requires engine identity and calculation settings')
    def identity(result):
        meta = result.metadata
        return (result.method, meta.get('engine',{}).get('sha256'),
                (meta.get('parameter_file') or {}).get('sha256'), meta.get('calculation_settings'))
    states = {}; flags = []
    for label, result in results.items():
        _validate_mapping(mol, result)
        if result.metadata['geometry_sha256'] != base['geometry_sha256']:
            raise ValueError('charge states must use the same fixed geometry')
        if result.metadata.get('calculation_status') != 'success':
            raise ValueError('charge states must record successful calculations')
        if identity(result) != identity(ref):
            raise ValueError('charge states must share method, engine, parameters and numerical settings')
        if result.number_of_molecular_orbitals != ref.number_of_molecular_orbitals:
            raise ValueError('charge states must use the same orbital basis size')
        delta = {'reference':0,'ionized':1,'attached':-1}[label]
        if result.metadata.get('charge') != base['charge'] + delta:
            raise ValueError(f'{label} has the wrong molecular charge')
        if result.number_of_electrons != ref.number_of_electrons - delta:
            raise ValueError(f'{label} has the wrong modeled electron count')
        _effective_charge_uhf(mol, result.metadata['charge'], result.metadata['uhf'])
        if result.number_of_unpaired_electrons != result.metadata['uhf']:
            raise ValueError('charge-state spin count disagrees with its provenance')
        energy = result.total_energy_eh
        if energy is None or not np.isfinite(energy): raise ValueError('all three energies must be finite')
        charges = _finite(result.partial_charges,'charges',1)
        if len(charges) != mol.GetNumAtoms(): raise ValueError('charges must match the indexed atoms')
        quality = assess_descriptor_quality(result)
        if quality['status'] == 'error': raise ValueError(f'{label} failed charge/electron conservation checks')
        flags.extend(dict(flag, state=label) for flag in quality['flags'])
        states[label] = _summary(result)
    groups = _fragments(fragments,mol.GetNumAtoms()) if fragments is not None else {}
    report = _empty_report(record_id,population_scheme,states)
    ip = (results['ionized'].total_energy_eh-ref.total_energy_eh)*HARTREE_TO_EV
    ea = (ref.total_energy_eh-results['attached'].total_energy_eh)*HARTREE_TO_EV
    mu = -(ip+ea)/2; eta = (ip-ea)/2
    values = report['values']; values.update(vertical_ip_eV=ip, vertical_ea_eV=ea,
            chemical_potential_eV=mu, electronegativity_eV=-mu, hardness_eV=eta)
    if eta > 1e-8:
        values.update(softness_per_eV=1/(2*eta), electrophilicity_eV=mu*mu/(2*eta))
    else:
        flags.append({'code':'nonpositive_or_tiny_hardness','severity':'warning',
                      'detail':'Reciprocal hardness descriptors are undefined or numerically unstable.'})
        report['missing'].update(softness_per_eV='nonpositive_or_tiny_hardness',electrophilicity_eV='nonpositive_or_tiny_hardness')
    if ea < 0:
        flags.append({'code':'negative_vertical_ea','severity':'warning',
                      'detail':'Attachment is energetically unfavorable at this geometry; no bound-anion claim.'})
    if ip < 0:
        flags.append({'code':'negative_vertical_ip','severity':'warning',
                      'detail':'Ionization lowers the calculated energy; inspect charge/spin states.'})
    if base['calculation_settings']['effective_electronic_temperature'] > 1e-6:
        flags.append({'code':'finite_electronic_temperature','severity':'warning',
                      'detail':'Energy differences use the configured finite-temperature libxtb energies.'})
    if population_scheme == 'charges':
        q0, qp, qm = [np.asarray(results[k].partial_charges) for k in ('reference','ionized','attached')]
        plus, minus = q0-qm, qp-q0
    else:
        n0, np_, nm = [np.asarray(calc_orbital_populations(results[k],scheme=population_scheme)['electron_populations_e'])
                       for k in ('reference','ionized','attached')]
        plus, minus = nm-n0, n0-np_
    dual, radical = plus-minus, (plus+minus)/2
    errors = {'f_plus':float(plus.sum()-1), 'f_minus':float(minus.sum()-1), 'dual':float(dual.sum())}
    if max(abs(v) for v in errors.values()) > 1e-5: raise ValueError('condensed Fukui sums violate electron conservation')
    report['atoms'] = [{'atom_index':i,'f_plus':float(plus[i]),'f_minus':float(minus[i]),
                        'f_radical':float(radical[i]),'dual':float(dual[i])} for i in range(len(plus))]
    report['fragments'] = {name:{'atom_indices':indices.tolist(),
        'f_plus':float(plus[indices].sum()), 'f_minus':float(minus[indices].sum()),
        'f_radical':float(radical[indices].sum()), 'dual':float(dual[indices].sum())}
        for name, indices in groups.items()}
    report['unassigned_atom_indices'] = sorted(set(range(len(plus)))-set(i for group in groups.values() for i in group))
    report['quality'] = {'status':'warning' if any(f['severity']=='warning' for f in flags) else 'ok',
                          'flags':flags, 'fukui_sum_errors':errors,
                          'spin_states_verified_as_ground_states':False}
    json.dumps(report, allow_nan=False)
    return report


def run_reactivity(mol, *, options=None, ionized_uhf=None, attached_uhf=None,
                   fragments=None, population_scheme='charges', interface='native',
                   record_id=None, return_errors=False, initialization='reference'):
    """Calculate fixed-geometry Q, Q+1, Q-1 states with explicit spin policy.

    Defaults use 300 K and independent copies of the reference wavefunction for
    both adjacent states. initialization='independent' uses fresh guesses.
    A closed-shell reference defaults to doublet ions. For an open-shell reference,
    both neighboring spin counts must be specified. No search for the lowest spin
    state is implied. Calculation failures raise ChargeStateCalculationError, or
    return an error report with all derived values null when return_errors=True.
    """
    from rdkit import Chem
    _scheme(population_scheme)
    if interface not in ('native','compiled'): raise ValueError('interface must be native or compiled')
    if initialization not in ('reference','independent'):
        raise ValueError('initialization must be reference or independent')
    if not isinstance(return_errors,bool): raise ValueError('return_errors must be Boolean')
    _empty_report(record_id,population_scheme,{})
    snapshot = Chem.Mol(mol)
    config = options or XTBOptions()
    _xyz_from_mol(snapshot, config.conf_id)
    charge, uhf = _effective_charge_uhf(snapshot,config.charge,config.uhf)
    if config.method not in ('GFN1-xTB','GFN2-xTB'): raise ValueError('vertical Delta-SCF requires self-consistent GFN1/2, not GFN0')
    if uhf and (ionized_uhf is None or attached_uhf is None):
        raise ValueError('open-shell references require explicit ionized_uhf and attached_uhf; no automatic spin search')
    spins = {'reference':uhf, 'ionized':1 if ionized_uhf is None else ionized_uhf,
             'attached':1 if attached_uhf is None else attached_uhf}
    charges = {'reference':charge,'ionized':charge+1,'attached':charge-1}
    for label in spins:
        _, spins[label] = _effective_charge_uhf(snapshot,charges[label],spins[label])
    if fragments is not None: _fragments(fragments,snapshot.GetNumAtoms())
    config = _resolve_options(replace(config,charge=charge,uhf=uhf,
        keep_overlap=config.keep_overlap or population_scheme != 'charges'))
    protocol = {'requested_options': config.to_dict(), 'interface': interface,
                'protocol_version': REACTIVITY_PROTOCOL_VERSION,
                'initialization': initialization,
                'charge_states': charges, 'unpaired_electrons': spins,
                'spin_selection': 'closed_shell_doublet_defaults'
                    if not uhf and ionized_uhf is None and attached_uhf is None
                    else 'explicit_adjacent_spins'}
    configurations = [replace(config, charge=charges[k], uhf=spins[k])
                      for k in ('reference','ionized','attached')]
    results, states = {}, {}
    calculations = _charge_state_results(snapshot, configurations, interface, initialization)
    with closing(calculations):
        for label in ('reference','ionized','attached'):
            try:
                result = next(calculations)
                results[label] = result; states[label] = _summary(result)
            except Exception as exc:
                states[label] = {'status':'error','charge':charges[label],'uhf':spins[label],
                                 'error':{'type':type(exc).__name__,'message':str(exc)}}
                for remaining in charges:
                    if remaining not in states: states[remaining] = {'status':'not_run','charge':charges[remaining],'uhf':spins[remaining]}
                error = ChargeStateCalculationError(label,exc,states)
                if not return_errors: raise error from exc
                report = _empty_report(record_id,population_scheme,states)
                report['analysis_settings'].update(protocol)
                report['error'] = {'type':type(error).__name__,'message':str(error),'state':label}
                report['quality'] = {'status':'error','flags':[{'code':'charge_state_failed','severity':'error','detail':str(error)}]}
                report['missing'] = dict.fromkeys(REACTIVITY_UNITS,'required_charge_state_failed')
                return report
    report = calc_reactivity_descriptors(snapshot,**results,fragments=fragments,
              population_scheme=population_scheme,record_id=record_id)
    report['analysis_settings'].update(protocol)
    return report


def _charge_state_results(mol, configurations, interface, initialization):
    if interface == 'compiled':
        from . import rdXTBTools
        raw_results = rdXTBTools._RunChargeStates(mol, configurations[0],
            configurations[1].uhf, configurations[2].uhf, initialization == 'reference')
        for raw in raw_results:
            if not raw.success: raise XTBRunError(raw.errorMessage)
            yield _as_result(raw)
    elif initialization == 'reference':
        with closing(_native_charge_states(mol, configurations)) as calculations:
            yield from calculations
    else:
        for selected in configurations:
            yield run_xtb(mol, **selected.to_dict())
