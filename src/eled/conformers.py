"""Uniform descriptive sensitivity across explicitly supplied conformers.

These statistics describe the supplied geometries, not a thermodynamic ensemble.
No embedding, optimization, deduplication, alignment or Boltzmann weighting.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import json
from numbers import Integral

import numpy as np
from rdkit import Chem

from .batch import _prepare, _key
from .errors import XTBEnvironmentError, XTBRunError
from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .reports import DESCRIPTOR_VERSION, _analysis_settings, _base_report, calc_descriptor_report, descriptor_columns
from .calculator import XTBCalculator

_KB_EH_PER_K = 3.166811563e-6


def conformer_descriptor_columns(*, charge=0, frontier_depth=2):
    """Ordered units for frame-invariant scalar sensitivity statistics.

    Charged molecules omit dipole norms because those depend on coordinate origin.
    Orbital indices, settings, categories and Cartesian components are excluded.
    Counts are retained as diagnostics of frontier classification changes.
    """
    if isinstance(charge, bool) or not isinstance(charge, Integral):
        raise ValueError('charge must be an integer')
    return {k: u for k, u in descriptor_columns(frontier_depth).items()
            if u is not None and u != 'index' and k != 'xtb_occupation_threshold'
            and not any(k == f'xtb_dipole_{a}_e_bohr' or k.endswith(f'centroid_{a}_angstrom') for a in 'xyz')
            and not (charge and k.startswith('xtb_dipole_'))}


def run_conformer_analysis(mol, *, conf_ids=None, options=None, cache=None,
                           record_id=None, columns=None, return_errors=True,
                           boltzmann_weighted=False, boltzmann_temperature=300.0,
                           **analysis):
    """Calculate one molecule's conformers and retain reports and coverage.

    ``conf_ids=None`` uses all actual IDs in RDKit order (not 0..N-1).
    ``options.conf_id`` must remain -1; choose conformers with ``conf_ids``.
    Only calculation failures become failed rows. Invalid inputs, environment,
    storage and analysis errors propagate. All geometries are validated first.
    Default columns are invariant under rigid transformations; no coordinate
    alignment is required. Each supplied conformer has equal weight, including
    duplicate geometries. ``std`` is population standard deviation (ddof=0).
    """
    if not isinstance(return_errors, bool): raise ValueError('return_errors must be Boolean')
    if not isinstance(boltzmann_weighted, bool): raise ValueError('boltzmann_weighted must be Boolean')
    if not np.isfinite(boltzmann_temperature) or boltzmann_temperature <= 0:
        raise ValueError('boltzmann_temperature must be finite and positive')
    if record_id is not None and not isinstance(record_id, str):
        raise ValueError('record_id must be a string or None')
    config = options or XTBOptions(keep_overlap=True)
    if config.conf_id != -1: raise ValueError('use conf_ids, not options.conf_id')
    settings = _analysis_settings(**analysis)
    snapshot = Chem.Mol(mol)
    ids = list(conf_ids) if conf_ids is not None else [c.GetId() for c in snapshot.GetConformers()]
    if not ids or any(isinstance(i, bool) or not isinstance(i, Integral) or i < 0 for i in ids):
        raise ValueError('conf_ids must contain nonnegative actual conformer IDs')
    ids = [int(i) for i in ids]
    if len(set(ids)) != len(ids): raise ValueError('duplicate conformer IDs')
    charge, uhf = _effective_charge_uhf(snapshot, config.charge, config.uhf)
    available = conformer_descriptor_columns(charge=charge, frontier_depth=settings['frontier_depth'])
    chosen = list(available) if columns is None else list(columns)
    if not chosen or len(set(chosen)) != len(chosen) or any(k not in available for k in chosen):
        raise ValueError('columns must be unique supported invariant descriptor names')
    units = {k: available[k] for k in chosen}
    inputs = [_input_metadata(snapshot, i, charge, uhf, _xyz_from_mol(snapshot, i)) for i in ids]
    reports, signature = [], None
    all_units = descriptor_columns(settings['frontier_depth'])
    for conf_id, metadata in zip(ids, inputs):
        resolved, _ = _prepare(snapshot, replace(config, conf_id=conf_id))
        key = _key(resolved, metadata) if cache is not None else None
        # Keep the catch inside the calculation callback: cache/storage failures
        # must never masquerade as a failed single point.
        def calculate():
            try:
                return XTBCalculator(**resolved.to_dict()).calculate(snapshot)
            except XTBEnvironmentError:
                raise
            except XTBRunError as exc:
                if not return_errors: raise
                return exc
        result = cache.get(key) if cache is not None else None
        if result is None:
            result = calculate()
            if cache is not None and not isinstance(result, Exception): cache.set(key, result)
        if isinstance(result, XTBRunError):
            report = _base_report(str(conf_id), all_units, deepcopy(settings))
            report['provenance'] = {**metadata, 'requested_options': resolved.to_dict()}
            report['quality'] = {'status': 'error', 'flags': [
                {'code': 'calculation_failed', 'severity': 'error', 'detail': str(result)}]}
            report['error'] = {'type': type(result).__name__, 'message': str(result)}
            report['missing'] = dict.fromkeys(all_units, 'calculation_failed')
        else:
            if (result.method != resolved.method or any(result.metadata.get(k) != metadata[k]
                    for k in ('conformer_id', 'charge', 'uhf', 'geometry_sha256', 'ordered_structure_sha256'))):
                raise ValueError('result does not match the requested conformer/protocol')
            current = {'method': result.method, 'backend': result.backend, 'version': result.version,
                       **{k: result.metadata.get(k) for k in ('engine', 'parameter_file',
                          'calculation_settings', 'initialization', 'optimize', 'solvent')}}
            if signature is not None and signature != current:
                raise ValueError('conformer results have inconsistent calculation provenance')
            signature = deepcopy(current)
            report = calc_descriptor_report(snapshot, result, record_id=str(conf_id), **settings)
        reports.append(report)
    failed = [r['error'] is not None or r['quality']['status'] == 'error' for r in reports]
    statistics = {}
    for name in chosen:
        values = [r['values'][name] for r, bad in zip(reports, failed) if not bad and r['values'][name] is not None]
        numbers = np.array(values, dtype=float)
        if not np.all(np.isfinite(numbers)): raise ValueError('nonfinite descriptor values')
        statistics[name] = {'n_total': len(ids), 'n_valid': len(numbers), 'n_failed': sum(failed),
            'n_missing': len(ids)-sum(failed)-len(numbers),
            **{key: float(func(numbers)) if len(numbers) else None for key, func in (
                ('mean', np.mean), ('std', np.std), ('min', np.min), ('max', np.max), ('range', np.ptp))}}
    # Variable-size sections are summarized by fixed atom/bond indices.  This
    # is well-defined because RDKit atom and bond indices are preserved across
    # conformers. Frontier group composition is tracked separately because a
    # near-degenerate group can change between geometries.
    def nested_statistics(sections, key_fields):
        records = {}
        for report, bad in zip(reports, failed):
            if bad: continue
            rows = []
            for section in sections:
                value = report.get(section)
                if isinstance(value, list): rows.extend((section, row) for row in value)
                elif isinstance(value, dict) and isinstance(value.get('atoms'), list):
                    rows.extend((section, row) for row in value['atoms'])
                elif isinstance(value, dict) and isinstance(value.get('bonds'), list):
                    rows.extend((section, row) for row in value['bonds'])
            for section, row in rows:
                identity = row.get('atom_index', row.get('bond_index'))
                if identity is None: continue
                for key, value in row.items():
                    if key in key_fields or not isinstance(value, (int, float)) or isinstance(value, bool): continue
                    if not np.isfinite(value): continue
                    records.setdefault((int(identity), f'{section}.{key}'), []).append(float(value))
        output = {}
        for (identity, key), values in sorted(records.items()):
            values = np.asarray(values, dtype=float)
            output.setdefault(str(identity), {})[key] = {
                'n_total': len(ids), 'n_valid': int(len(values)), 'n_failed': sum(failed),
                'n_missing': len(ids) - sum(failed) - len(values),
                'mean': float(values.mean()), 'std': float(values.std()),
                'min': float(values.min()), 'max': float(values.max())}
        return output
    atom_statistics = nested_statistics(('atoms', 'shell_populations', 'pi'), {'atom_index'})
    bond_statistics = nested_statistics(('bonds', 'mayer_bonds', 'frontier_bonds', 'pi'),
                                        {'bond_index', 'begin_atom_index', 'end_atom_index'})
    frontier_compositions = []
    for report, bad in zip(reports, failed):
        if bad or not report.get('frontier_bonds'): continue
        frontier_compositions.append({row['frontier_group']: tuple(row['orbital_indices'])
                                      for row in report['frontier_bonds'].get('bonds', [])})
    composition_changed = bool(frontier_compositions and
                               any(c != frontier_compositions[0] for c in frontier_compositions[1:]))
    weighted_statistics = None
    if boltzmann_weighted:
        energies = np.asarray([r['values'].get('xtb_total_energy_hartree') if not bad else np.nan
                               for r, bad in zip(reports, failed)], dtype=float)
        valid = np.isfinite(energies)
        weights = np.zeros(len(ids), dtype=float)
        if valid.any():
            shifted = energies[valid] - np.nanmin(energies[valid])
            weights[valid] = np.exp(-shifted / (_KB_EH_PER_K * boltzmann_temperature))
            weights /= weights.sum()
        weighted_statistics = {}
        for name in chosen:
            vals = np.asarray([r['values'][name] if not bad else np.nan for r, bad in zip(reports, failed)], dtype=float)
            good = np.isfinite(vals) & (weights > 0)
            if not good.any():
                weighted_statistics[name] = {'mean': None, 'std': None, 'n_valid': 0,
                                             'weighting': 'electronic_energy_not_free_energy'}
                continue
            w = weights[good] / weights[good].sum(); x = vals[good]; mean = float(np.sum(w*x))
            weighted_statistics[name] = {'mean': mean, 'std': float(np.sqrt(np.sum(w*(x-mean)**2))),
                                         'n_valid': int(good.sum()), 'weighting': 'electronic_energy_not_free_energy'}
    complete = all(v['n_valid'] == len(ids) for v in statistics.values())
    report = {'schema_version': 2, 'descriptor_version': DESCRIPTOR_VERSION,
        'report_type': 'conformer_sensitivity', 'record_id': record_id,
        'status': 'failed' if all(failed) else 'complete' if complete else 'partial',
        'n_conformers': len(ids), 'n_successful': len(ids)-sum(failed), 'n_failed': sum(failed),
        'weighting': 'uniform_supplied_geometries', 'std_ddof': 0,
        'units': units, 'statistics': statistics, 'conformer_ids': ids, 'inputs': inputs,
        'analysis_settings': settings, 'requested_options': config.to_dict(),
        'calculation_provenance': signature, 'reports': reports,
        'atom_statistics': atom_statistics, 'bond_statistics': bond_statistics,
        'frontier_group_composition_changed': composition_changed,
        'sensitivity': {
            'conformer_dispersion': {k: v['std'] for k, v in statistics.items()},
            'gfn1_vs_gfn2': {'status': 'not_requested', 'reason': 'run two explicit analyses to compare methods'},
            'electronic_temperature': {'status': 'not_requested', 'reason': 'no second temperature calculation supplied'},
        },
        'boltzmann': ({'enabled': True, 'temperature_K': float(boltzmann_temperature),
                       'energy_definition': 'xTB electronic energy; not a free energy',
                       'statistics': weighted_statistics} if boltzmann_weighted else
                      {'enabled': False, 'statistics': None}),
        'excluded_columns': [k for k in all_units if k not in chosen]}
    json.dumps(report, allow_nan=False)
    return report
