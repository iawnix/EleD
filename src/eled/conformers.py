"""Descriptor statistics across explicitly supplied conformers.

Equal-weight statistics are the default. Optional Boltzmann-style statistics use
xTB electronic energies as weights and are reported separately from the default
geometry statistics.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
import json
from numbers import Integral

import numpy as np
from rdkit import Chem

from .batch import _key, _prepare
from .calculator import XTBCalculator
from .density_descriptors import _pi_options
from .errors import XTBEnvironmentError, XTBRunError
from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .reports import (
    DESCRIPTOR_VERSION,
    CalcDescriptorReport,
    DescriptorColumns,
    _analysis_settings,
    _base_report,
)

_KB_EH_PER_K = 3.166811563e-6


def ConformerDescriptorColumns(*, charge=0, frontierDepth=2):
    """Ordered units for frame-invariant scalar sensitivity statistics."""
    if isinstance(charge, bool) or not isinstance(charge, Integral):
        raise ValueError('charge must be an integer')
    return {
        k: u
        for k, u in DescriptorColumns(frontierDepth).items()
        if u is not None
        and u != 'index'
        and k != 'xtb_occ_cut'
        and not any(k == f'xtb_dip_{a}_ebohr' or k.endswith(f'c{a}_a') for a in 'xyz')
        and not (charge and k.startswith('xtb_dip_'))
    }


def RunConformerAnalysis(
    mol,
    *,
    confIds=None,
    options=None,
    cache=None,
    recordId=None,
    columns=None,
    returnErrors=True,
    boltzmannWeighted=False,
    boltzmannTemperature=300.0,
    **analysis,
):
    """Calculate one molecule's conformers and retain reports and coverage."""
    if not isinstance(returnErrors, bool):
        raise ValueError('returnErrors must be Boolean')
    if not isinstance(boltzmannWeighted, bool):
        raise ValueError('boltzmannWeighted must be Boolean')
    if isinstance(boltzmannTemperature, bool) or not np.isfinite(boltzmannTemperature) or boltzmannTemperature <= 0:
        raise ValueError('boltzmannTemperature must be finite and positive')
    if recordId is not None and not isinstance(recordId, str):
        raise ValueError('recordId must be a string or None')
    config = options or XTBOptions(keepOverlap=True)
    if config.confId != -1:
        raise ValueError('use confIds, not options.confId')
    settings = _analysis_settings(**analysis)
    snapshot = Chem.Mol(mol)
    _pi_options(
        settings['piScheme'],
        settings['piPlaneTolAngstrom'],
        settings['piRef'],
        natoms=snapshot.GetNumAtoms(),
    )
    ids = list(confIds) if confIds is not None else [c.GetId() for c in snapshot.GetConformers()]
    if not ids or any(isinstance(i, bool) or not isinstance(i, Integral) or i < 0 for i in ids):
        raise ValueError('confIds must contain nonnegative actual conformer IDs')
    ids = [int(i) for i in ids]
    if len(set(ids)) != len(ids):
        raise ValueError('duplicate conformer IDs')
    charge, unpaired = _effective_charge_uhf(snapshot, config.charge, config.numUnpairedElectrons)
    available = ConformerDescriptorColumns(charge=charge, frontierDepth=settings['frontierDepth'])
    chosen = list(available) if columns is None else list(columns)
    if not chosen or len(set(chosen)) != len(chosen) or any(k not in available for k in chosen):
        raise ValueError('columns must be unique supported invariant descriptor names')
    units = {k: available[k] for k in chosen}
    inputs = [_input_metadata(snapshot, i, charge, unpaired, _xyz_from_mol(snapshot, i)) for i in ids]
    reports, signature = [], None
    all_units = DescriptorColumns(settings['frontierDepth'])
    for conf_id, metadata in zip(ids, inputs):
        resolved, _ = _prepare(snapshot, replace(config, confId=conf_id))
        key = _key(resolved, metadata) if cache is not None else None

        def calculate():
            try:
                return XTBCalculator(**resolved.ToDict()).Calculate(snapshot)
            except XTBEnvironmentError:
                raise
            except XTBRunError as exc:
                if not returnErrors:
                    raise
                return exc

        result = cache.Get(key) if cache is not None else None
        if result is None:
            result = calculate()
            if cache is not None and not isinstance(result, Exception):
                cache.Set(key, result)
        if isinstance(result, XTBRunError):
            report = _base_report(str(conf_id), all_units, deepcopy(settings))
            report['provenance'] = {**metadata, 'requested_options': resolved.ToDict()}
            report['quality'] = {
                'status': 'error',
                'flags': [{'code': 'calculation_failed', 'severity': 'error', 'detail': str(result)}],
            }
            report['error'] = {'type': type(result).__name__, 'message': str(result)}
            report['missing'] = dict.fromkeys(all_units, 'calculation_failed')
        else:
            if (result.method != resolved.method or any(
                result.metadata.get(k) != metadata[k]
                for k in ('conformer_id', 'charge', 'uhf', 'geometry_sha256', 'ordered_structure_sha256')
            )):
                raise ValueError('result does not match the requested conformer/protocol')
            current = {
                'method': result.method,
                'backend': result.backend,
                'version': result.version,
                **{
                    k: result.metadata.get(k)
                    for k in ('engine', 'parameter_file', 'calculation_settings', 'initialization', 'optimize', 'solvent')
                },
            }
            if signature is not None and signature != current:
                raise ValueError('conformer results have inconsistent calculation provenance')
            signature = deepcopy(current)
            report = CalcDescriptorReport(snapshot, result, recordId=str(conf_id), **settings)
        reports.append(report)
    failed = [r['error'] is not None or r['quality']['status'] == 'error' for r in reports]
    statistics = {}
    for name in chosen:
        values = [r['values'][name] for r, bad in zip(reports, failed) if not bad and r['values'][name] is not None]
        numbers = np.array(values, dtype=float)
        if not np.all(np.isfinite(numbers)):
            raise ValueError('nonfinite descriptor values')
        statistics[name] = {
            'n_total': len(ids),
            'n_valid': len(numbers),
            'n_failed': sum(failed),
            'n_missing': len(ids) - sum(failed) - len(numbers),
            **{
                key: float(func(numbers)) if len(numbers) else None
                for key, func in (
                    ('mean', np.mean),
                    ('std', np.std),
                    ('min', np.min),
                    ('max', np.max),
                    ('range', np.ptp),
                )
            },
        }
    atom_series = _local_series(reports, failed, 'atoms')
    bond_series = _local_series(reports, failed, 'bonds')
    atom_statistics = _local_statistics(atom_series, failed)
    bond_statistics = _local_statistics(bond_series, failed)
    frontier_compositions = [
        {label: (tuple(report['frontier'][label]['orb_idx']) if report['frontier'][label] else None)
         for label in ('homo', 'lumo')}
        for report, bad in zip(reports, failed)
        if not bad and report.get('frontier')
    ]
    composition_changed = bool(
        frontier_compositions and any(c != frontier_compositions[0] for c in frontier_compositions[1:])
    )
    weighted_statistics = None
    boltzmann_weights = [0.0] * len(ids)
    if boltzmannWeighted:
        energies = np.asarray(
            [r['values'].get('xtb_energy_ha') if not bad else np.nan for r, bad in zip(reports, failed)],
            dtype=float,
        )
        valid = np.isfinite(energies)
        weights = np.zeros(len(ids), dtype=float)
        if valid.any():
            shifted = energies[valid] - np.nanmin(energies[valid])
            weights[valid] = np.exp(-shifted / (_KB_EH_PER_K * boltzmannTemperature))
            weights /= weights.sum()
            boltzmann_weights = weights.tolist()
        weighted_statistics = {}
        for name in chosen:
            vals = [r['values'][name] if not bad else None for r, bad in zip(reports, failed)]
            weighted_statistics[name] = _summary(vals, failed, weights=weights, energy_valid=valid)
        weighted_atoms = _local_statistics(atom_series, failed, weights=weights, energy_valid=valid)
        weighted_bonds = _local_statistics(bond_series, failed, weights=weights, energy_valid=valid)
    complete = all(v['n_valid'] == len(ids) for v in statistics.values())
    report = {
        'schema_version': 2,
        'descriptor_version': DESCRIPTOR_VERSION,
        'report_type': 'conformer_sensitivity',
        'record_id': recordId,
        'status': 'failed' if all(failed) else 'complete' if complete else 'partial',
        'n_conformers': len(ids),
        'n_successful': len(ids) - sum(failed),
        'n_failed': sum(failed),
        'weighting': 'uniform_supplied_geometries',
        'std_ddof': 0,
        'units': units,
        'statistics': statistics,
        'conformer_ids': ids,
        'inputs': inputs,
        'analysis_settings': settings,
        'requested_options': config.ToDict(),
        'calculation_provenance': signature,
        'reports': reports,
        'atom_statistics': atom_statistics,
        'bond_statistics': bond_statistics,
        'atom_units': _local_units(reports, atom_series),
        'bond_units': _local_units(reports, bond_series),
        'frontier_group_composition_changed': composition_changed,
        'boltzmann': (
            {
                'enabled': True,
                'temperature_K': float(boltzmannTemperature),
                'energy_definition': 'xTB electronic energy; not a free energy',
                'weights': boltzmann_weights,
                'statistics': weighted_statistics,
                'atom_statistics': weighted_atoms,
                'bond_statistics': weighted_bonds,
            }
            if boltzmannWeighted
            else {'enabled': False, 'statistics': None}
        ),
        'excluded_columns': [k for k in all_units if k not in chosen],
    }
    json.dumps(report, allow_nan=False)
    return report


def _local_series(reports, failed, entity):
    sections = (
        ('atoms', 'shell_populations', 'pi')
        if entity == 'atoms'
        else ('bonds', 'mayer_bonds', 'frontier_bonds', 'pi')
    )
    index_key = 'atom_index' if entity == 'atoms' else 'bond_index'
    excluded = {'atom_index', 'bond_index', 'begin_atom_index', 'end_atom_index'}
    series = {}
    for i, (report, bad) in enumerate(zip(reports, failed)):
        if bad:
            continue
        for section in sections:
            value = report.get(section)
            rows = value if isinstance(value, list) else value.get(entity, []) if isinstance(value, dict) else []
            for row in rows:
                identity = row[index_key]
                prefix = section + ('.' + row['frontier_group'] if section == 'frontier_bonds' else '')
                for field, value in row.items():
                    if field in excluded or (
                        value is not None and (isinstance(value, bool) or not isinstance(value, (int, float)))
                    ):
                        continue
                    if field in {'normal', 'missing'}:
                        continue
                    key = (str(identity), prefix + '.' + field)
                    series.setdefault(key, [None] * len(reports))[i] = value
    return series


def _local_units(reports, series):
    units = {
        'atoms.q_e': 'e',
        'atoms.cm5_q_e': 'e',
        'atoms.wbo_sum': 'dimensionless',
        'bonds.wbo': 'dimensionless',
    }
    for report in reports:
        for section in ('shell_populations', 'mayer_bonds', 'frontier_bonds', 'pi'):
            data = report.get(section) or {}
            for field, unit in data.get('units', {}).items():
                labels = ('frontier_bonds.homo', 'frontier_bonds.lumo') if section == 'frontier_bonds' else (section,)
                for label in labels:
                    units[label + '.' + field] = unit
    return {field: units[field] for field in sorted({key[1] for key in series})}


def _summary(values, failed, *, weights=None, energy_valid=None):
    x = np.asarray(values, dtype=float)
    present = ~np.asarray(failed) & ~np.isnan(x)
    if np.any(present & ~np.isfinite(x)):
        raise ValueError('nonfinite descriptor values')
    if energy_valid is not None:
        present &= energy_valid
    count = int(present.sum())
    out = {
        'n_total': len(x),
        'n_valid': count,
        'n_failed': sum(failed),
        'n_missing': len(x) - sum(failed) - count,
    }
    if weights is None:
        out.update({
            name: float(func(x[present])) if count else None
            for name, func in [('mean', np.mean), ('std', np.std), ('min', np.min), ('max', np.max)]
        })
    else:
        mass = float(weights[present].sum())
        out.update(weight_sum=mass, weighting='electronic_energy_not_free_energy', mean=None, std=None)
        if mass > 0:
            w = weights[present] / mass
            mean = float(w @ x[present])
            out.update(mean=mean, std=float(np.sqrt(w @ (x[present] - mean) ** 2)))
        else:
            out['missing_reason'] = 'no_positive_weight_for_available_values'
    return out


def _local_statistics(series, failed, **kwargs):
    output = {}
    for (identity, field), values in sorted(series.items()):
        output.setdefault(identity, {})[field] = _summary(values, failed, **kwargs)
    return output
