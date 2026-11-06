"""Input validation and leakage-safe splits for atom-labelled EAS records.

This module intentionally does not download a dataset or run xTB.  It turns a
locally supplied public-data export into a frozen, auditable record table.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable, Mapping

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

_ALLOWED = {1, 5, 6, 7, 8, 9, 14, 15, 16, 17, 35, 53}


def _load_rows(path):
    path = Path(path)
    if path.suffix.lower() == '.csv':
        with path.open(newline='', encoding='utf-8') as handle:
            return list(csv.DictReader(handle))
    text = path.read_text(encoding='utf-8')
    if path.suffix.lower() in ('.jsonl', '.ndjson'):
        return [json.loads(line) for line in text.splitlines() if line.strip()]
    payload = json.loads(text)
    if isinstance(payload, list): return payload
    if isinstance(payload, dict):
        for key in ('records', 'data', 'reactions', 'molecules'):
            if isinstance(payload.get(key), list): return payload[key]
    raise ValueError('input must be a JSON list/object, JSONL, or CSV')


def _atoms(value):
    if value is None or value == '': return []
    if isinstance(value, str):
        value = value.strip()
        if value.startswith('['): value = json.loads(value)
        else: value = [x for x in value.replace(';', ',').split(',') if x.strip()]
    if not isinstance(value, (list, tuple)):
        raise ValueError('reactive atom labels must be a list of integer atom indices')
    out = []
    for item in value:
        if isinstance(item, bool): raise ValueError('Boolean atom labels are invalid')
        out.append(int(item))
    return out


def _canonical_smiles(mol):
    return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)


def _scaffold(mol):
    core = MurckoScaffold.GetScaffoldForMol(Chem.RemoveHs(mol))
    value = Chem.MolToSmiles(core, canonical=True, isomericSmiles=False)
    return value or 'acyclic:' + _canonical_smiles(Chem.RemoveHs(mol))


def normalize_eas_records(rows: Iterable[Mapping], *, source='local'):
    """Validate and normalize public EAS records.

    Accepted label keys are ``reactive_atoms``, ``target_atoms``,
    ``label_atoms`` and ``labels``. Labels are zero-based atom indices in the
    supplied SMILES. Invalid rows are returned in ``excluded`` with a reason;
    no row is silently repaired.
    """
    records, excluded, seen = [], [], set()
    for ordinal, raw in enumerate(rows):
        try:
            if not isinstance(raw, Mapping): raise ValueError('record is not an object')
            identifier = str(raw.get('id', raw.get('record_id', ordinal))).strip()
            smiles = str(raw.get('smiles', '')).strip()
            if not identifier or not smiles: raise ValueError('id and smiles are required')
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: raise ValueError('SMILES sanitization failed')
            if len(Chem.GetMolFrags(mol)) != 1: raise ValueError('disconnected mixture is excluded')
            if any(atom.GetAtomicNum() not in _ALLOWED for atom in mol.GetAtoms()):
                raise ValueError('molecule contains an element outside the frozen allowed set')
            map_key = next((k for k in ('reactive_atom_maps', 'target_atom_maps') if k in raw), None)
            if map_key is not None:
                requested_maps = _atoms(raw[map_key])
                by_map = {atom.GetAtomMapNum(): atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomMapNum()}
                if any(value not in by_map for value in requested_maps):
                    raise ValueError('target atom map number is absent from the molecule')
                labels = sorted(set(by_map[value] for value in requested_maps))
                key = map_key
            else:
                key = next((k for k in ('reactive_atoms', 'target_atoms', 'label_atoms', 'labels') if k in raw), None)
                labels = None
            if key is None: raise ValueError('missing atom-level target labels')
            if labels is None: labels = sorted(set(_atoms(raw[key])))
            if not labels or any(i < 0 or i >= mol.GetNumAtoms() for i in labels):
                raise ValueError('target atom index is outside the molecule')
            if identifier in seen: raise ValueError('duplicate record id')
            seen.add(identifier)
            canonical = _canonical_smiles(mol)
            records.append({'id': identifier, 'smiles': smiles,
                            'canonical_smiles': canonical,
                            'reactive_atoms': labels,
                            'scaffold': _scaffold(mol),
                            'source': source,
                            'source_row': ordinal})
        except Exception as exc:
            excluded.append({'source_row': ordinal, 'reason': str(exc), 'record': dict(raw) if isinstance(raw, Mapping) else raw})
    return records, excluded


def _union_find(n):
    parent = list(range(n))
    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]; i = parent[i]
        return i
    def union(a, b):
        a, b = find(a), find(b)
        if a != b: parent[b] = a
    return parent, find, union


def cluster_keys(records, threshold=0.7):
    """Return single-linkage Morgan clusters at the frozen Tanimoto threshold."""
    if not 0 < threshold <= 1: raise ValueError('threshold must be in (0, 1]')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=False)
    fps = []
    for row in records:
        mol = Chem.MolFromSmiles(row['smiles'])
        fps.append(generator.GetFingerprint(mol))
    parent, find, union = _union_find(len(records))
    for i in range(len(fps)):
        for j in range(i):
            if DataStructs.TanimotoSimilarity(fps[i], fps[j]) >= threshold: union(i, j)
    roots = {}
    keys = []
    for i in range(len(records)):
        root = find(i); keys.append(roots.setdefault(root, f'cluster-{len(roots):05d}'))
    return keys


def assign_eas_splits(records, *, seed=20261106, train_fraction=0.8, valid_fraction=0.1):
    """Assign scaffold-plus-cluster groups to train/valid/test deterministically."""
    if not 0 < train_fraction < 1 or not 0 <= valid_fraction < 1 or train_fraction + valid_fraction >= 1:
        raise ValueError('split fractions must leave a positive test fraction')
    clusters = cluster_keys(records)
    parent, find, union = _union_find(len(records))
    owners = {}
    for index, (row, cluster) in enumerate(zip(records, clusters)):
        # Shared scaffolds and shared fingerprint clusters form connected
        # components. This enforces both leakage constraints even when a
        # scaffold contains multiple fingerprint clusters (or vice versa).
        for key in (('scaffold', row['scaffold']), ('cluster', cluster)):
            if key in owners: union(index, owners[key])
            else: owners[key] = index
    groups = {}
    for index, row in enumerate(records): groups.setdefault(find(index), []).append(row['id'])
    ordered = sorted(groups.values(), key=lambda ids: (-len(ids), hashlib.sha256(f'{seed}:{ids}'.encode()).hexdigest()))
    totals = {'train': 0, 'valid': 0, 'test': 0}; n = len(records)
    result = {}
    for ids in ordered:
        candidates = ('train', 'valid', 'test')
        split = min(candidates, key=lambda name: (totals[name] / max(n, 1) /
            ({'train': train_fraction, 'valid': valid_fraction, 'test': 1-train_fraction-valid_fraction}[name]),
            candidates.index(name)))
        for identifier in ids: result[identifier] = split
        totals[split] += len(ids)
    return result


def prepare_eas_dataset(path, *, source='local', seed=20261106):
    rows, excluded = normalize_eas_records(_load_rows(path), source=source)
    splits = assign_eas_splits(rows, seed=seed) if rows else {}
    for row in rows: row['split'] = splits[row['id']]
    return {'protocol_id': 'eas_regioselectivity_v1', 'source': source,
            'input_path': str(Path(path)), 'n_input': len(rows) + len(excluded),
            'n_selected': len(rows), 'n_excluded': len(excluded),
            'records': rows, 'excluded': excluded,
            'split_counts': {name: sum(r['split'] == name for r in rows) for name in ('train','valid','test')}}

