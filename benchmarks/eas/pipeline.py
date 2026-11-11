"""Feature-archive contract for the EAS model benchmark.

This module deliberately consumes JSON produced by ``tools/run_eas_features``.
It never calls RDKit or xTB and never changes the public descriptor API. Missing
values remain ``None`` with a parallel mask so model code can fit imputation
inside training folds.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping


_ROOT = Path(__file__).resolve().parents[2]
_FEATURE_PROTOCOL = _ROOT / 'protocol.json'
_MODEL_PROTOCOL = Path(__file__).with_name('model_protocol.json')
_HYBRIDIZATION = {
    'UNSPECIFIED': 0, 'S': 1, 'SP': 2, 'SP2': 3, 'SP3': 4,
    'SP3D': 5, 'SP3D2': 6, 'OTHER': 7,
}

_BASE_FIELDS = ('atomic_number', 'formal_charge', 'degree', 'total_hydrogens',
                'is_aromatic', 'is_in_ring', 'hybridization')
_GFN_FIELDS = ('q_e', 'wbo_sum')
_SHELL_FIELDS = (
    'mul_s_e', 'mul_p_e', 'mul_d_e', 'mul_ps',
    'low_s_e', 'low_p_e', 'low_d_e', 'low_ps',
)
_REACTIVITY_FIELDS = ('f_plus', 'f_minus', 'f_radical', 'dual',
                      's_loc_1ev', 'omega_loc_ev', 'n_loc_1ev')
_BOND_FIELDS = {
    'wbo': ('bonds', 'wbo', None),
    'mbo': ('mayer_bonds', 'mbo', None),
    'homo_mul_ovlp_e': ('frontier_bonds', 'mul_ovlp_e', 'homo'),
    'lumo_mul_ovlp_e': ('frontier_bonds', 'mul_ovlp_e', 'lumo'),
    'pi_bo': ('pi', 'pi_bo', None),
}


def _canonical_sha256(payload: Mapping) -> str:
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'),
                         ensure_ascii=False).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


def _read_json(path):
    path = Path(path)
    return json.loads(path.read_text(encoding='utf-8'))


def load_feature_archive(path):
    """Load and validate an archive against the frozen feature protocol."""
    archive = _read_json(path)
    validate_feature_archive(archive)
    return archive


def validate_feature_archive(archive, *, feature_protocol=None, model_protocol=None):
    """Validate archive identity and required per-record provenance.

    The validator intentionally does not require successful calculations for
    every record. Failed records are valid observations for coverage analysis,
    but they are excluded by the model matrix builder.
    """
    if not isinstance(archive, Mapping):
        raise TypeError('feature archive must be an object')
    protocol = _read_json(feature_protocol or _FEATURE_PROTOCOL)
    model = _read_json(model_protocol or _MODEL_PROTOCOL)
    expected_hash = _canonical_sha256(protocol)
    if archive.get('protocol_id') != protocol.get('protocol_id'):
        raise ValueError('feature archive protocol_id does not match protocol.json')
    if archive.get('protocol_sha256') != expected_hash:
        raise ValueError('feature archive protocol_sha256 does not match protocol.json')
    if model.get('feature_protocol_id') != protocol.get('protocol_id'):
        raise ValueError('model protocol targets a different feature protocol')
    if model.get('feature_protocol_sha256') != expected_hash:
        raise ValueError('model protocol feature hash does not match protocol.json')
    records = archive.get('records')
    if not isinstance(records, list):
        raise ValueError('feature archive records must be a list')
    identifiers = set()
    for row in records:
        if not isinstance(row, Mapping):
            raise ValueError('feature archive record must be an object')
        identifier = row.get('id')
        if not isinstance(identifier, str) or not identifier:
            raise ValueError('feature archive record id must be a nonempty string')
        if identifier in identifiers:
            raise ValueError(f'duplicate feature archive record: {identifier}')
        identifiers.add(identifier)
        if not isinstance(row.get('reactive_atoms'), list):
            raise ValueError(f'{identifier}: reactive_atoms must be a list')
        if row.get('split') not in ('train', 'valid', 'test'):
            raise ValueError(f'{identifier}: invalid primary split')
        fold = row.get('cv_fold')
        if isinstance(fold, bool) or not isinstance(fold, int) or fold < 0:
            raise ValueError(f'{identifier}: invalid cv_fold')
    return {'protocol_sha256': expected_hash, 'n_records': len(records),
            'model_protocol_id': model.get('model_protocol_id')}


def _number(value):
    return value if isinstance(value, (int, float)) and not isinstance(value, bool) else None


def _atom_stat(analysis, atom_index, section, field, statistic='mean'):
    values = analysis.get('atom_statistics') or {}
    row = values.get(str(atom_index), {})
    stats = row.get(f'{section}.{field}') or {}
    if statistic == 'range':
        low, high = _number(stats.get('min')), _number(stats.get('max'))
        return high - low if low is not None and high is not None else None
    return _number(stats.get(statistic))


def _reactivity_stat(record, atom_index, field):
    reports = (record.get('analyses') or {}).get('gfn2_reactivity_300K') or []
    values = []
    for report in reports:
        for row in report.get('atoms') or []:
            if row.get('atom_index') == atom_index and _number(row.get(field)) is not None:
                values.append(float(row[field]))
    return sum(values) / len(values) if values else None


def _bond_stat(record, atom_index, field):
    """Mean over available incident bond/conformer pairs, preserving frontier groups."""
    values = []
    section, source_field, frontier = _BOND_FIELDS[field]
    analysis = (record.get('analyses') or {}).get('gfn2_300K') or {}
    for report in analysis.get('reports') or []:
        if report.get('error') is not None or (report.get('quality') or {}).get('status') == 'error':
            continue
        data = report.get(section)
        rows = data if isinstance(data, list) else (data or {}).get('bonds', [])
        for row in rows:
            if frontier is not None and row.get('frontier_group') != frontier:
                continue
            if atom_index in (row.get('begin_atom_index'), row.get('end_atom_index')):
                value = _number(row.get(source_field))
                if value is not None:
                    values.append(float(value))
    return sum(values) / len(values) if values else None


def _structure_features(record, atom_index):
    rows = record.get('atom_features') or []
    if atom_index >= len(rows):
        return {name: None for name in _BASE_FIELDS}
    source = rows[atom_index]
    hybridization = str(source.get('hybridization', 'UNSPECIFIED')).upper()
    return {name: (_HYBRIDIZATION.get(hybridization, _HYBRIDIZATION['OTHER'])
                   if name == 'hybridization' else _number(source.get(name)))
            for name in _BASE_FIELDS}


def atom_rows(archive, *, include_failed=False) -> list[dict]:
    """Return one labeled row per heavy atom and molecule.

    ``include_failed=False`` removes records without a successful primary
    GFN2 analysis. It does not impute missing descriptor values within a
    successful record.
    """
    validate_feature_archive(archive)
    rows = []
    for record in archive['records']:
        analysis = (record.get('analyses') or {}).get('gfn2_300K')
        if not include_failed and (not isinstance(analysis, Mapping) or
                                   analysis.get('status') == 'failed'):
            continue
        n_atoms = len(record.get('atom_features') or [])
        labels = {int(i) for i in record.get('reactive_atoms', [])}
        for atom_index in range(n_atoms):
            rows.append({'record_id': record['id'], 'atom_index': atom_index,
                         'target': int(atom_index in labels),
                         'split': record['split'], 'cv_fold': record['cv_fold'],
                         'scaffold': record.get('scaffold'),
                         'record': record})
    return rows


def feature_matrix(archive, group='A', *, include_failed=False):
    """Build a JSON-serializable atom matrix for one declared feature group."""
    if group not in set('ABCDEFG'):
        raise ValueError("group must be one of 'A'...'G'")
    output = []
    for row in atom_rows(archive, include_failed=include_failed):
        record, atom_index = row['record'], row['atom_index']
        features = {}
        if group in 'ABCDEFG':
            features.update({f'rdkit_{k}': v for k, v in
                             _structure_features(record, atom_index).items()})
        if group in 'BCDEFG':
            analysis = (record.get('analyses') or {}).get('gfn1_300K') or {}
            features.update({f'gfn1_{field}': _atom_stat(analysis, atom_index, 'atoms', field)
                             for field in _GFN_FIELDS})
        if group in 'CDEFG':
            analysis = (record.get('analyses') or {}).get('gfn2_300K') or {}
            features.update({f'gfn2_{field}': _atom_stat(analysis, atom_index, 'atoms', field)
                             for field in _GFN_FIELDS})
        if group in 'DEFG':
            analysis = (record.get('analyses') or {}).get('gfn2_300K') or {}
            features.update({f'gfn2_{field}': _atom_stat(analysis, atom_index, 'shell_populations', field)
                             for field in _SHELL_FIELDS})
        if group in 'EF':
            features.update({f'gfn2_{field}': _bond_stat(record, atom_index, field)
                             for field in _BOND_FIELDS})
        if group == 'F':
            features.update({f'reactivity_{field}': _reactivity_stat(record, atom_index, field)
                             for field in _REACTIVITY_FIELDS})
        if group == 'G':
            analysis = (record.get('analyses') or {}).get('gfn2_300K') or {}
            for field in _GFN_FIELDS + _SHELL_FIELDS:
                for statistic in ('std', 'min', 'max', 'range'):
                    features[f'gfn2_{field}_{statistic}'] = _atom_stat(
                        analysis, atom_index,
                        'atoms' if field in _GFN_FIELDS else 'shell_populations', field, statistic)
        output.append({k: v for k, v in row.items() if k != 'record'} |
                      {'features': features,
                       'missing': {key: value is None for key, value in features.items()}})
    return output
