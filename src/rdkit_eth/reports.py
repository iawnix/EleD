"""Versioned descriptor reports, observable quality checks and table export."""
from __future__ import annotations

from copy import deepcopy
import csv
import io
import json
import os
from pathlib import Path
import tempfile

import numpy as np

from .electronic_descriptors import (
    _as_result, _orbitals, _nonnegative, calc_xtb_descriptors,
    calc_xtb_atom_descriptors, calc_xtb_bond_descriptors, calc_xtb_fragment_descriptors,
)
from .result import XTBResult
from .populations import calc_orbital_populations, calc_frontier_descriptors
from .density_descriptors import (calc_density_matrix, calc_shell_populations,
    calc_mayer_bond_descriptors, calc_frontier_bond_descriptors, calc_pi_descriptors)

DESCRIPTOR_VERSION = '1.1'
SCHEMA_VERSION = 2


def descriptor_columns(frontier_depth=2):
    """Ordered scalar-column -> unit mapping, fixed for a chosen frontier depth."""
    empty = XTBResult('', None, None, None, None, None, None)
    values = calc_xtb_descriptors(empty, frontier_depth=frontier_depth)
    units = {}
    for name in values:
        if name.endswith('_hartree'): unit = 'Hartree'
        elif name.endswith('_ev'): unit = 'eV'
        elif name.endswith('_e_bohr'): unit = 'e*Bohr'
        elif name.endswith('_debye'): unit = 'Debye'
        elif '_charge_' in name: unit = 'e'
        elif name == 'xtb_occupation_threshold': unit = 'dimensionless'
        elif name.startswith('xtb_num_'): unit = 'count'
        elif name.endswith('_index'): unit = 'index'
        elif name == 'xtb_wiberg_pair_sum': unit = 'dimensionless'
        else: unit = None
        units[name] = unit
    for label in ('homo', 'lumo'):
        for field, unit in [('num_orbitals','count'), ('energy_min_ev','eV'), ('energy_max_ev','eV'),
                            ('inverse_participation_ratio','dimensionless'), ('participation_ratio','dimensionless'),
                            ('entropy_nats','nat'), *[(f'centroid_{a}_angstrom','angstrom') for a in 'xyz']]:
            units[f'frontier_{label}_{field}'] = unit
    units['frontier_homo_lumo_atom_overlap'] = 'dimensionless'
    units['frontier_homo_lumo_atom_centroid_separation_angstrom'] = 'angstrom'
    return units


def assess_descriptor_quality(result, *, occupation_threshold=1e-7, fractional_tolerance=1e-4,
                              conservation_tolerance=1e-5):
    """Report observable issues and optional verified engine SCC diagnostics."""
    result = _as_result(result)
    for name, value in [('fractional_tolerance', fractional_tolerance),
                        ('conservation_tolerance', conservation_tolerance)]:
        _nonnegative(value, name)
    _, occ, occupied, virtual = _orbitals(result, occupation_threshold)
    if np.any(occ > 2 + conservation_tolerance):
        raise ValueError('spin-summed orbital occupations cannot exceed two')
    flags = []
    def flag(code, severity, detail):
        flags.append({'code': code, 'severity': severity, 'detail': detail})
    fractional = np.flatnonzero(np.abs(occ-np.round(occ)) > fractional_tolerance).tolist()
    if fractional:
        flag('fractional_occupations', 'warning', 'Orbital occupations differ from integer values.')
    counts = [int(np.sum(occ > t)) for t in
              (occupation_threshold * .1, occupation_threshold, max(occupation_threshold*10, 1e-8))]
    if len(set(counts)) > 1:
        flag('frontier_threshold_sensitive', 'warning', 'Occupied/virtual classification changes within the tested threshold window.')
    if not len(occupied) or not len(virtual):
        flag('frontier_unavailable', 'warning', 'No occupied or no virtual orbital under this threshold.')
    if result.number_of_unpaired_electrons:
        flag('open_shell_spin_summed', 'info', 'No separate alpha/beta population analysis or spin-state search.')
    expected_charge = result.metadata.get('charge')
    charge_error = None
    if len(result.partial_charges) and expected_charge is not None:
        charge_error = float(sum(result.partial_charges) - expected_charge)
        if not np.isfinite(charge_error): raise ValueError('charges must be finite')
        if abs(charge_error) > conservation_tolerance:
            flag('charge_conservation_failed', 'error', 'Atomic charges do not sum to the requested molecular charge.')
    if expected_charge:
        flag('charged_dipole_origin_dependent', 'info', 'Dipole components and norm depend on the coordinate origin.')
    electron_error = None
    if result.number_of_electrons is not None:
        electron_error = float(sum(occ)-result.number_of_electrons)
        if abs(electron_error) > conservation_tolerance:
            flag('electron_conservation_failed', 'error', 'Occupations do not sum to the reported modeled electron count.')
    settings = result.metadata.get('calculation_settings')
    if not settings:
        flag('calculation_settings_unknown', 'warning', 'Result does not record the numerical calculation settings.')
    elif settings['electronic_temperature'] != settings['effective_electronic_temperature']:
        flag('temperature_floor_applied', 'info', 'xTB uses a minimum electronic temperature of 1e-6 K.')
    scf = result.metadata.get('scf', {})
    if scf.get('available') and not scf['criteria_satisfied']:
        flag('scf_criteria_not_satisfied', 'error', 'Final SCC diagnostics do not satisfy engine thresholds.')
    severity = {f['severity'] for f in flags}
    return {'status': 'error' if 'error' in severity else 'warning' if 'warning' in severity else 'ok',
            'flags': flags, 'fractional_orbital_indices': fractional,
            'occupied_counts_threshold_window': counts,
            'charge_sum_error_e': charge_error, 'electron_sum_error_e': electron_error,
            'scf_iterations': scf.get('iterations'), 'scf_residual': scf.get('mixer_rms'),
            'scf_details': scf or None,
            'scf_details_unavailable_reason': None if scf.get('available') else scf.get('reason', 'not_exposed_by_libxtb_c_api'),
            'calculation_status': result.metadata.get('calculation_status', 'unknown')}


def descriptor_sensitivity(reference, *, gfn1=None, gfn2=None,
                           temperature_300=None, temperature_1000=None):
    """Compare already-produced reports without invoking an engine.

    The returned per-column diagnostics deliberately keep the three effects
    separate: conformer dispersion, GFN1/GFN2 absolute difference, and the
    300 K/1000 K absolute difference.  ``None`` means the comparison report
    was not supplied or the value is missing; it is never treated as zero.
    """
    if not isinstance(reference, dict):
        raise TypeError('reference must be a descriptor report or conformer report')
    if 'statistics' in reference:
        base = {k: v.get('std') for k, v in reference['statistics'].items()}
    else:
        base = {k: None for k in reference.get('values', {})}
    def values(report):
        if report is None: return {}
        return report.get('values', {}) if 'values' in report else {k: v.get('mean') for k, v in report.get('statistics', {}).items()}
    ref_values = values(reference)
    g1, g2 = values(gfn1), values(gfn2)
    t3, t10 = values(temperature_300), values(temperature_1000)
    keys = list(ref_values) if ref_values else list(base)
    output = {}
    for key in keys:
        def delta(a, b):
            try:
                return abs(float(a) - float(b)) if a is not None and b is not None else None
            except (TypeError, ValueError): return None
        output[key] = {
            'conformer_dispersion': base.get(key),
            'gfn1_vs_gfn2_abs_difference': delta(g1.get(key), g2.get(key)),
            'electronic_temperature_300_vs_1000_abs_difference': delta(t3.get(key), t10.get(key)),
            'missing': {
                'gfn1_vs_gfn2': gfn1 is None or g2.get(key) is None or g1.get(key) is None,
                'electronic_temperature': temperature_300 is None or temperature_1000 is None or t3.get(key) is None or t10.get(key) is None,
            },
        }
    return {'descriptor_version': DESCRIPTOR_VERSION,
            'definition': 'absolute differences; conformer dispersion is population std',
            'values': output}


def _base_report(record_id, columns, analysis):
    if record_id is not None and not isinstance(record_id, str):
        raise ValueError('record_id must be a string or None')
    return {'schema_version': SCHEMA_VERSION, 'descriptor_version': DESCRIPTOR_VERSION,
            'report_type': 'single_point', 'record_id': record_id,
            'units': columns, 'values': dict.fromkeys(columns), 'analysis_settings': analysis,
            'provenance': {}, 'quality': {}, 'missing': {}, 'error': None,
            'atoms': None, 'bonds': None, 'fragments': None, 'populations': None, 'frontier': None,
            'density': None, 'shell_populations': None, 'mayer_bonds': None,
            'frontier_bonds': None, 'pi': None}


def _analysis_settings(*, frontier_depth=2, occupation_threshold=1e-7,
                       degeneracy_ev=1e-3, occupation_tolerance=1e-6,
                       include_orbitals=True, fractional_tolerance=1e-4,
                       conservation_tolerance=1e-5):
    descriptor_columns(frontier_depth)
    if not isinstance(include_orbitals, bool):
        raise ValueError('include_orbitals must be Boolean')
    values = dict(occupation_threshold=occupation_threshold, degeneracy_ev=degeneracy_ev,
                  occupation_tolerance=occupation_tolerance, fractional_tolerance=fractional_tolerance,
                  conservation_tolerance=conservation_tolerance)
    for key, value in values.items():
        _nonnegative(value, key)
    return dict(frontier_depth=int(frontier_depth), include_orbitals=include_orbitals,
                **{key:float(value) for key,value in values.items()})


def calc_descriptor_report(mol, result, *, record_id=None, fragments=None, frontier_depth=2,
                           occupation_threshold=1e-7, degeneracy_ev=1e-3,
                           occupation_tolerance=1e-6, include_orbitals=True,
                           fractional_tolerance=1e-4, conservation_tolerance=1e-5):
    """Extract all descriptor sections without running an engine.

    Missing retained basis gives explicit nulls/reasons; inconsistent arrays or
    mismatched molecule provenance raise errors. CSV scalar columns are fixed
    for frontier_depth; variable-size atom/fragment/orbital data remain nested.
    """
    from .mapping import _validate_mapping
    analysis = _analysis_settings(frontier_depth=frontier_depth, occupation_threshold=occupation_threshold,
        degeneracy_ev=degeneracy_ev, occupation_tolerance=occupation_tolerance,
        include_orbitals=include_orbitals, fractional_tolerance=fractional_tolerance,
        conservation_tolerance=conservation_tolerance)
    result = _as_result(result)
    _validate_mapping(mol, result)
    report = _base_report(record_id, descriptor_columns(frontier_depth), analysis)
    report['provenance'] = deepcopy(dict(result.metadata))
    report['quality'] = assess_descriptor_quality(result, occupation_threshold=occupation_threshold,
            fractional_tolerance=fractional_tolerance, conservation_tolerance=conservation_tolerance)
    scalars = calc_xtb_descriptors(result, frontier_depth=frontier_depth, occupation_threshold=occupation_threshold)
    report['values'].update(scalars)
    for key, value in scalars.items():
        if value is None: report['missing'][key] = 'unavailable_from_result_or_orbital_range'
    if len(result.partial_charges):
        report['atoms'] = calc_xtb_atom_descriptors(result)
    else: report['missing']['atoms'] = 'charges_not_reported'
    if len(result.wiberg_bond_orders):
        report['bonds'] = calc_xtb_bond_descriptors(mol, result)
    else: report['missing']['bonds'] = 'bond_orders_not_reported'
    if fragments is not None:
        report['fragments'] = calc_xtb_fragment_descriptors(result, fragments)
    else: report['missing']['fragments'] = 'not_requested'
    if include_orbitals and len(result.overlap_matrix):
        report['populations'] = {scheme: calc_orbital_populations(result, scheme=scheme)
                                for scheme in ('lowdin', 'mulliken')}
        frontier = report['frontier'] = calc_frontier_descriptors(mol, result, fragments=fragments,
            occupation_threshold=occupation_threshold, degeneracy_ev=degeneracy_ev,
            occupation_tolerance=occupation_tolerance)
        for label in ('homo', 'lumo'):
            group = frontier[label]
            if group is None: continue
            for field in ('num_orbitals','energy_min_ev','energy_max_ev','inverse_participation_ratio',
                          'participation_ratio','entropy_nats'):
                report['values'][f'frontier_{label}_{field}'] = group[field]
            for axis, value in zip('xyz', group['atom_centroid_angstrom']):
                report['values'][f'frontier_{label}_centroid_{axis}_angstrom'] = value
        for key in ('homo_lumo_atom_overlap','homo_lumo_atom_centroid_separation_angstrom'):
            report['values']['frontier_' + key] = frontier[key]
        for name, function in (
            ('density', calc_density_matrix),
            ('shell_populations', calc_shell_populations),
            ('mayer_bonds', lambda r: calc_mayer_bond_descriptors(mol, r)),
            ('frontier_bonds', lambda r: calc_frontier_bond_descriptors(mol, r)),
            ('pi', lambda r: calc_pi_descriptors(mol, r)),
        ):
            try:
                report[name] = function(result)
            except Exception as exc:
                report['missing'][name] = f'{type(exc).__name__}: {exc}'
                report['quality']['flags'].append({'code': f'{name}_unavailable',
                    'severity': 'warning', 'detail': str(exc)})
                if report['quality']['status'] == 'ok': report['quality']['status'] = 'partial'
    else:
        reason = 'basis_not_retained' if include_orbitals else 'not_requested'
        report['missing'].update(populations=reason, frontier=reason, density=reason,
                                 shell_populations=reason, mayer_bonds=reason,
                                 frontier_bonds=reason, pi=reason)
        if include_orbitals and report['quality']['status'] == 'ok': report['quality']['status'] = 'partial'
    for key, value in report['values'].items():
        if key.startswith('frontier_') and value is None:
            report['missing'][key] = report['missing'].get('frontier', 'frontier_group_unavailable')
    json.dumps(report, allow_nan=False)
    return report


def run_descriptor_batch(molecules, *, options=None, cache=None, record_ids=None, fragments=None,
                         max_workers=1, timeout=None, isolate=False, **analysis):
    """Ordered calculation/report rows; failures contain null features and diagnostics."""
    from .batch import run_xtb_batch
    from .options import XTBOptions
    mols = list(molecules)
    ids = list(record_ids) if record_ids is not None else list(map(str, range(len(mols))))
    if len(ids) != len(mols) or any(v is not None and not isinstance(v, str) for v in ids):
        raise ValueError('record_ids must match molecules and contain strings or None')
    analysis = _analysis_settings(**analysis)
    columns = descriptor_columns(analysis['frontier_depth'])
    config = options or XTBOptions(keep_overlap=True)
    results = run_xtb_batch(mols, options=config,
                            cache=cache, return_exceptions=True,
                            max_workers=max_workers, timeout=timeout, isolate=isolate)
    reports = []
    for mol, result, identifier in zip(mols, results, ids):
        try:
            if isinstance(result, Exception): raise result
            report = calc_descriptor_report(mol, result, record_id=identifier, fragments=fragments, **analysis)
        except Exception as exc:
            report = _base_report(identifier, columns.copy(), dict(analysis))
            report['quality'] = {'status':'error', 'flags':[{'code':'calculation_or_analysis_failed',
                                   'severity':'error', 'detail': str(exc)}]}
            report['error'] = {'type':type(exc).__name__, 'message':str(exc)}
            report['missing'] = dict.fromkeys(columns, 'calculation_or_analysis_failed')
            report['provenance'] = {'requested_options': config.to_dict()}
            if not isinstance(result, Exception):
                report['provenance'].update(deepcopy(dict(result.metadata)))
        reports.append(report)
    return reports


def export_descriptor_reports(reports, path, *, format=None):
    """Write consistent reports to CSV (scalar table) or JSONL (all sections).

    CSV includes JSON provenance, settings, flags and missing reasons per row.
    It never silently concatenates incompatible descriptor versions/columns.
    """
    reports = list(reports)
    if not reports: raise ValueError('at least one report is required')
    path = Path(path); kind = format or path.suffix.lstrip('.').lower()
    if kind not in ('csv', 'jsonl'): raise ValueError('format must be csv or jsonl')
    first = reports[0]
    signature = (first['schema_version'], first['descriptor_version'], first['report_type'], first['units'])
    for report in reports:
        if signature != (report['schema_version'],report['descriptor_version'],report['report_type'],report['units']):
            raise ValueError('reports have incompatible versions, types or scalar schemas')
        if set(report['values']) != set(first['units']): raise ValueError('values do not match scalar schema')
        json.dumps(report, allow_nan=False)
    if kind == 'jsonl':
        text = ''.join(json.dumps(r, allow_nan=False) + '\n' for r in reports)
    else:
        output = io.StringIO(newline='')
        headers = ['record_id','schema_version','descriptor_version','report_type','quality_status',
                   *first['units'], 'units_json','analysis_settings_json','provenance_json','quality_json','missing_json','error_json','execution_json']
        writer = csv.DictWriter(output, fieldnames=headers); writer.writeheader()
        for report in reports:
            row = {k:report[k] for k in ('record_id','schema_version','descriptor_version','report_type')}
            row.update(report['values']); row['quality_status'] = report['quality']['status']
            row.update({k+'_json':json.dumps(report[k],allow_nan=False) for k in
                        ('units','analysis_settings','provenance','quality','missing','error')})
            row['execution_json'] = json.dumps(report.get('execution'), allow_nan=False)
            writer.writerow(row)
        text = output.getvalue()
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', newline='', dir=path.parent, delete=False) as handle:
        temporary = Path(handle.name)
        try: handle.write(text)
        except BaseException:
            temporary.unlink(missing_ok=True); raise
    try: os.replace(temporary, path)
    finally: temporary.unlink(missing_ok=True)
    return path
