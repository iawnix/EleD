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

_PROTOCOL_PATHS = (Path(__file__).with_name('protocol.json'), Path(__file__).resolve().parents[2] / 'protocol.json')


def load_frozen_protocol():
    """Load the repository's frozen EAS protocol.

    The protocol is an input artifact, rather than a second set of Python
    defaults.  Keeping this loader here lets manifests carry a digest of the
    exact JSON used for selection and splitting.
    """
    for protocol_path in _PROTOCOL_PATHS:
        if protocol_path.exists():
            payload = json.loads(protocol_path.read_text(encoding='utf-8'))
            break
    else:
        raise RuntimeError(f'frozen protocol is missing: {_PROTOCOL_PATHS[-1]}')
    if payload.get('protocol_id') != 'eas_regioselectivity_v1':
        raise ValueError('unsupported frozen EAS protocol')
    return payload


def frozen_protocol_sha256(protocol=None):
    """Return the SHA256 of canonical protocol JSON (without filesystem paths)."""
    payload = load_frozen_protocol() if protocol is None else protocol
    encoded = json.dumps(payload, sort_keys=True, separators=(',', ':'), ensure_ascii=False).encode('utf-8')
    return hashlib.sha256(encoded).hexdigest()


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
    protocol = load_frozen_protocol()
    allowed = {Chem.GetPeriodicTable().GetAtomicNumber(symbol)
               for symbol in protocol['selection']['allowed_elements']}
    records, excluded, seen, seen_ids = [], [], set(), set()
    for ordinal, raw in enumerate(rows):
        try:
            if not isinstance(raw, Mapping): raise ValueError('record is not an object')
            identifier = str(raw.get('id', raw.get('record_id', ordinal))).strip()
            smiles = str(raw.get('smiles', '')).strip()
            if not identifier or not smiles: raise ValueError('id and smiles are required')
            mol = Chem.MolFromSmiles(smiles)
            if mol is None: raise ValueError('SMILES sanitization failed')
            if len(Chem.GetMolFrags(mol)) != 1: raise ValueError('disconnected mixture is excluded')
            if any(atom.GetAtomicNum() not in allowed for atom in mol.GetAtoms()):
                raise ValueError('molecule contains an element outside the frozen allowed set')
            event = raw.get('event_type', raw.get('reaction_class'))
            if event is not None and str(event).strip().lower() not in {
                    'eas', 'electrophilic aromatic substitution', 'electrophilic_aromatic_substitution'}:
                raise ValueError('record is not an electrophilic aromatic substitution event')
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
            canonical = _canonical_smiles(mol)
            scaffold = _scaffold(mol)
            reaction_id = str(raw.get('reaction_id', identifier)).strip()
            formal_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
            requested_charge = raw.get('charge', raw.get('formal_charge', formal_charge))
            requested_uhf = raw.get('uhf', raw.get('spin_uhf'))
            if isinstance(requested_charge, bool):
                raise ValueError('molecular charge must be an integer')
            requested_charge = int(requested_charge)
            if requested_charge != formal_charge:
                raise ValueError('explicit charge does not match the SMILES formal charge')
            radical_uhf = sum(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms())
            if requested_uhf is None:
                if radical_uhf:
                    raise ValueError('open-shell records must provide explicit uhf')
                requested_uhf = radical_uhf
            if isinstance(requested_uhf, bool) or int(requested_uhf) < 0:
                raise ValueError('uhf must be a nonnegative integer')
            requested_uhf = int(requested_uhf)
            ionized_uhf = raw.get('ionized_uhf')
            attached_uhf = raw.get('attached_uhf')
            if requested_uhf and (ionized_uhf is None or attached_uhf is None):
                raise ValueError('open-shell records must provide ionized_uhf and attached_uhf')
            for name, value in (('ionized_uhf', ionized_uhf), ('attached_uhf', attached_uhf)):
                if value is not None:
                    if isinstance(value, bool) or int(value) < 0:
                        raise ValueError(f'{name} must be a nonnegative integer')
                    value = int(value)
                if name == 'ionized_uhf': ionized_uhf = value
                else: attached_uhf = value
            dedup_key = (reaction_id, scaffold)
            if identifier in seen_ids:
                raise ValueError('duplicate record id')
            if dedup_key in seen: raise ValueError('duplicate reaction/scaffold record')
            seen.add(dedup_key); seen_ids.add(identifier)
            records.append({'id': identifier, 'smiles': smiles,
                            'canonical_smiles': canonical,
                            'reactive_atoms': labels,
                            'scaffold': scaffold,
                            'reaction_id': reaction_id,
                            'charge': requested_charge,
                            'uhf': requested_uhf,
                            'ionized_uhf': ionized_uhf,
                            'attached_uhf': attached_uhf,
                            'dedup_key': f'{reaction_id}:{scaffold}',
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


def cluster_keys(records, threshold=0.7, *, radius=2, fp_size=1024, include_chirality=False):
    """Return single-linkage Morgan clusters at the frozen Tanimoto threshold."""
    if not 0 < threshold <= 1: raise ValueError('threshold must be in (0, 1]')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size,
                                                          includeChirality=include_chirality)
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


def _leakage_groups(records, *, threshold=0.7, radius=2, fp_size=1024,
                    include_chirality=False):
    clusters = cluster_keys(records, threshold=threshold, radius=radius,
                            fp_size=fp_size, include_chirality=include_chirality)
    parent, find, union = _union_find(len(records))
    owners = {}
    for index, (row, cluster) in enumerate(zip(records, clusters)):
        for key in (('scaffold', row['scaffold']), ('cluster', cluster)):
            if key in owners:
                union(index, owners[key])
            else:
                owners[key] = index
    groups = {}
    for index, row in enumerate(records):
        groups.setdefault(find(index), []).append(row['id'])
    return list(groups.values())


def _frozen_cluster_settings():
    cluster = load_frozen_protocol()['split']['cluster']
    return dict(threshold=cluster['threshold'], radius=cluster['radius'],
                fp_size=cluster['fp_size'], include_chirality=cluster['include_chirality'])


def assign_eas_splits(records, *, seed=20261106, train_fraction=0.8, valid_fraction=0.1):
    """Assign scaffold-plus-cluster groups to train/valid/test deterministically."""
    if not 0 < train_fraction < 1 or not 0 <= valid_fraction < 1 or train_fraction + valid_fraction >= 1:
        raise ValueError('split fractions must leave a positive test fraction')
    groups = _leakage_groups(records, **_frozen_cluster_settings())
    ordered = sorted(groups, key=lambda ids: (-len(ids), hashlib.sha256(f'{seed}:{ids}'.encode()).hexdigest()))
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


def assign_eas_folds(records, *, seed=20261106, n_splits=5):
    """Assign deterministic leakage-safe grouped folds for model selection."""
    if isinstance(n_splits, bool) or not isinstance(n_splits, int) or n_splits < 2:
        raise ValueError('n_splits must be an integer >= 2')
    groups = _leakage_groups(records, **_frozen_cluster_settings())
    ordered = sorted(groups, key=lambda ids: (-len(ids), hashlib.sha256(
        f'{seed}:{ids}'.encode()).hexdigest()))
    totals = [0] * n_splits
    result = {}
    for ids in ordered:
        fold = min(range(n_splits), key=lambda i: (totals[i], i))
        for identifier in ids:
            result[identifier] = fold
        totals[fold] += len(ids)
    return result


def prepare_eas_dataset(path, *, source='local', seed=None):
    protocol = load_frozen_protocol()
    if seed is None:
        seed = protocol['split']['seed']
    elif seed != protocol['split']['seed']:
        raise ValueError(f'seed must equal frozen protocol value {protocol["split"]["seed"]}')
    rows, excluded = normalize_eas_records(_load_rows(path), source=source)
    fractions = protocol['split']['fractions']
    splits = assign_eas_splits(rows, seed=seed, train_fraction=fractions['train'],
                               valid_fraction=fractions['valid']) if rows else {}
    folds = assign_eas_folds(rows, seed=seed, n_splits=protocol['split']['secondary_folds']) if rows else {}
    for row in rows:
        row['split'] = splits[row['id']]
        row['cv_fold'] = folds[row['id']]
    return {'protocol_id': 'eas_regioselectivity_v1', 'source': source,
            'protocol_sha256': frozen_protocol_sha256(protocol),
            'split_seed': seed,
            'input_path': str(Path(path)), 'n_input': len(rows) + len(excluded),
            'n_selected': len(rows), 'n_excluded': len(excluded),
            'records': rows, 'excluded': excluded,
            'split_counts': {name: sum(r['split'] == name for r in rows) for name in ('train','valid','test')}}

