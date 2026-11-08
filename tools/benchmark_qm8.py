#!/usr/bin/env python3
"""Frozen QM8 E1-CC2 pilot; optional sklearn tooling, never a runtime backend.

No tuning on test data. Source geometries, not xTB-optimized structures.
Downloads stay in --data; outputs include cohort, splits, failures and predictions.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
import platform
import re
import tarfile
import time
import urllib.request

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from eled import XTBEnvironmentError, XTBOptions, XTBRunError, CalcDescriptorReport, XTBCalculator
from eled.input import _xyz_from_mol, _geometry_hash
from eled.persistence import _implementation_identity
from eled.backends.native import HARTREE_TO_EV

URL = 'https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/gdb8.tar.gz'
SHA256 = 'e668b58f8dfba329bc20cd7e0eb3944e8a73a5d9db37160f686a8b31b20548a0'
SEEDS = (17, 29, 43)
TWO_D = ('MolWt', 'MolLogP', 'TPSA', 'NumHDonors', 'NumHAcceptors',
         'NumRotatableBonds', 'RingCount', 'NumAromaticRings', 'NumAliphaticRings',
         'FractionCSP3', 'HeavyAtomCount', 'NHOHCount', 'NOCount', 'LabuteASA',
         'BalabanJ', 'BertzCT', 'Chi0v', 'Chi1v', 'Kappa1', 'Kappa2', 'Kappa3',
         'NumValenceElectrons', 'NumHeteroatoms', 'MolMR')
BASIC = ('xtb_total_energy_hartree', 'xtb_homo_ev', 'xtb_lumo_ev', 'xtb_gap_ev',
         'xtb_dipole_norm_debye', 'xtb_charge_min_e', 'xtb_charge_max_e',
         'xtb_charge_std_e', 'xtb_wiberg_pair_sum')
# Immutable 0.9 evidence stores the former column spellings. This tuple is
# only for validating that archival manifest; new calculations use BASIC.
FROZEN_BASIC = ('xtb_total_energy_Eh', 'xtb_homo_eV', 'xtb_lumo_eV', 'xtb_gap_eV',
                'xtb_dipole_norm_debye', 'xtb_charge_min_e', 'xtb_charge_max_e',
                'xtb_charge_std_e', 'xtb_wiberg_pair_sum')
ORBITAL = tuple(f'frontier_{label}_{field}' for label in ('homo', 'lumo')
                for field in ('inverse_participation_ratio', 'participation_ratio', 'entropy_nats')) + (
    'frontier_homo_lumo_atom_overlap', 'frontier_homo_lumo_atom_centroid_separation_angstrom')


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def read_source(root):
    root.mkdir(parents=True, exist_ok=True)
    path = root / 'gdb8.tar.gz'
    if not path.exists():
        temporary = path.with_suffix('.download')
        urllib.request.urlretrieve(URL, temporary)
        if hashlib.sha256(temporary.read_bytes()).hexdigest() != SHA256:
            raise ValueError('download checksum differs from the frozen source')
        temporary.replace(path)
    if hashlib.sha256(path.read_bytes()).hexdigest() != SHA256:
        raise ValueError('source checksum mismatch')
    with tarfile.open(path) as archive:
        sdf = archive.extractfile('qm8.sdf').read().decode()
        properties = archive.extractfile('qm8.sdf.csv').read().decode()
    reader = csv.reader(io.StringIO(properties)); header = next(reader)
    # Other repeated PBE0 names refer to different basis sets; never DictReader.
    assert header[:2] == ['gdb9_index', 'E1-CC2'] and header.count('E1-CC2') == 1
    targets = {}
    for row in reader:
        identifier = str(int(float(row[0])))
        if identifier in targets: raise ValueError('duplicate target ID')
        targets[identifier] = float(row[1]) * HARTREE_TO_EV
    candidates, excluded, seen_ids, seen_chemistry = [], [], set(), set()
    for block in sdf.split('$$$$'):
        block = block.lstrip('\r\n')
        if not block.strip(): continue
        match = re.match(r'gdb\s+(\d+)\b', block)
        if not match: raise ValueError('missing source ID')
        identifier = match[1]
        if identifier in seen_ids: raise ValueError('duplicate geometry ID')
        seen_ids.add(identifier)
        if identifier not in targets: raise ValueError('geometry has no target')
        mol = Chem.MolFromMolBlock(block, removeHs=False)
        if mol is None:
            excluded.append({'id': identifier, 'reason': 'RDKit_sanitization_failed'}); continue
        heavy = Chem.RemoveHs(mol)
        # Conservative deduplication also groups stereoisomers; choose first source ID.
        identity = Chem.MolToSmiles(heavy, isomericSmiles=False)
        reason = None
        if identity in seen_chemistry: reason = 'duplicate_nonstereo_identity'
        elif Chem.GetFormalCharge(mol) or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
            reason = 'outside_neutral_closed_shell_protocol'
        elif not np.isfinite(targets[identifier]): reason = 'nonfinite_target'
        else:
            try: _xyz_from_mol(mol, -1)
            except ValueError as exc: reason = 'invalid_geometry: ' + str(exc)
        if reason:
            excluded.append({'id': identifier, 'reason': reason}); continue
        seen_chemistry.add(identity)
        scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol=heavy, includeChirality=False)
        candidates.append(({'id': identifier, 'smiles': Chem.MolToSmiles(heavy),
                            'identity': identity, 'scaffold': scaffold,
                            'target_eV': targets[identifier],
                            'geometry_sha256': _geometry_hash(_xyz_from_mol(mol, -1))}, mol))
    if seen_ids != set(targets): raise ValueError('unmatched source IDs')
    candidates.sort(key=lambda pair: hashlib.sha256(('qm8-pilot-v1:' + pair[0]['identity']).encode()).hexdigest())
    return candidates, excluded, len(targets)


def split_indices(rows, kind, seed):
    from sklearn.model_selection import GroupShuffleSplit
    groups = [r['identity'] if kind == 'random' else r['scaffold'] for r in rows]
    # All acyclic compounds share the empty scaffold; never split this group.
    train, test = next(GroupShuffleSplit(n_splits=1, test_size=.2, random_state=seed).split(rows, groups=groups))
    assert not set(groups[i] for i in train) & set(groups[i] for i in test)
    return train, test


def validate_splits(rows, splits):
    """Validate frozen, exhaustive partitions without resampling missing records."""
    index = {r['id']: i for i, r in enumerate(rows)}
    if len(index) != len(rows): raise ValueError('duplicate row IDs')
    seen = set()
    for split in splits:
        key = split['kind'], split['seed']
        if key in seen or split['kind'] not in ('random', 'scaffold'):
            raise ValueError('invalid or repeated split')
        seen.add(key)
        train, test = split['train'], split['test']
        if (not train or not test or len(set(train)) != len(train) or len(set(test)) != len(test)
                or set(train) & set(test) or set(train) | set(test) != set(index)):
            raise ValueError('split must partition the exact cohort without duplicates')
        group = 'identity' if split['kind'] == 'random' else 'scaffold'
        if {rows[index[i]][group] for i in train} & {rows[index[i]][group] for i in test}:
            raise ValueError('chemical group crosses train/test')
    if not seen: raise ValueError('at least one split is required')


def fit_ablation(rows, output, *, splits=None):
    import sklearn
    import scipy
    from sklearn.ensemble import ExtraTreesRegressor
    from sklearn.impute import SimpleImputer
    from sklearn.pipeline import make_pipeline
    features = ['rdkit', 'rdkit_basic', 'rdkit_basic_orbital']
    arrays = [np.asarray([r['rdkit'] for r in rows], dtype=float)]
    arrays.append(np.column_stack([arrays[0], [r['basic'] for r in rows]]))
    arrays.append(np.column_stack([arrays[1], [r['orbital'] for r in rows]]))
    y = np.array([r['target_eV'] for r in rows])
    if splits is None:
        splits = []
        for kind in ('random', 'scaffold'):
            for seed in SEEDS:
                train, test = split_indices(rows, kind, seed)
                splits.append({'kind': kind, 'seed': seed, 'train': [rows[i]['id'] for i in train],
                               'test': [rows[i]['id'] for i in test]})
    validate_splits(rows, splits)
    index = {r['id']: i for i, r in enumerate(rows)}
    metrics, predictions = [], []
    for split in splits:
        kind, seed = split['kind'], split['seed']
        train = np.array([index[i] for i in split['train']])
        test = np.array([index[i] for i in split['test']])
        for name, x in zip(features, arrays):
            model = make_pipeline(SimpleImputer(strategy='median', keep_empty_features=True),
                ExtraTreesRegressor(n_estimators=256, min_samples_leaf=2, max_features=1.,
                                    random_state=seed, n_jobs=2))
            started = time.perf_counter(); model.fit(x[train], y[train]); pred = model.predict(x[test])
            metrics.append({'split': kind, 'seed': seed, 'features': name, 'n_features': x.shape[1],
                'n_train': len(train), 'n_test': len(test), 'mae_eV': float(np.mean(abs(pred-y[test]))),
                'rmse_eV': float(np.sqrt(np.mean((pred-y[test])**2))),
                'fit_predict_seconds': time.perf_counter()-started})
            predictions.extend({'split': kind, 'seed': seed, 'features': name, 'id': rows[i]['id'],
                                'target_eV': y[i], 'prediction_eV': float(p)} for i, p in zip(test, pred))
            print(json.dumps(metrics[-1]), flush=True)
    write_json(output/'splits.json', splits)
    write_json(output/'predictions.json', predictions)
    return {'model': 'ExtraTreesRegressor; 256 trees, min_samples_leaf=2, max_features=1.0; no tuning',
            'imputation': 'train-only median, keep_empty_features=True', 'metrics': metrics,
            'versions': {'python': platform.python_version(), 'rdkit': rdBase.rdkitVersion,
                         'numpy': np.__version__, 'sklearn': sklearn.__version__, 'scipy': scipy.__version__}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=Path('build/qm8'))
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--size', type=int, default=1024)
    args = parser.parse_args()
    if args.size < 32: parser.error('size must be >=32')
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()): raise ValueError('output directory must be empty; preserve prior evidence')
    candidates, excluded, source_count = read_source(args.data)
    if args.size > len(candidates): raise ValueError('size exceeds eligible source')
    selected = candidates[:args.size]
    options = XTBOptions(keep_overlap=True)
    manifest = {'protocol': 'qm8-pilot-v1', 'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'source': URL, 'sha256': SHA256,
        'reference_doi': '10.1063/1.4928757', 'target': 'E1-CC2, RI-CC2/def2-TZVP, Hartree converted to eV',
        'hartree_to_eV': HARTREE_TO_EV, 'source_count': source_count, 'eligible_count': len(candidates),
        'exclusions': excluded, 'selection': 'first N SHA256(qm8-pilot-v1:nonstereo canonical identity); target-blind',
        'geometry': 'provided ground-state geometry, SDF coordinates rounded to 4 decimals; no reoptimization',
        'options': options.to_dict(), 'selected': [row for row, mol in selected],
        'rdkit_features': list(TWO_D) + ['Morgan radius2, 1024bits, no chirality'],
        'basic_features': BASIC, 'orbital_features': ORBITAL,
        'seeds': SEEDS, 'test_fraction_of_groups': .2,
        'implementation': _implementation_identity()}
    write_json(args.output/'manifest.json', manifest)  # freeze before any calculation/model
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=False)
    rows, failures = [], []
    started = time.perf_counter()
    # Large population matrices are deliberately not duplicated in benchmark evidence.
    with gzip.open(args.output/'features.jsonl.gz', 'wt') as handle:
        for n, (metadata, mol) in enumerate(selected):
            heavy = Chem.RemoveHs(mol); row = dict(metadata)
            t = time.perf_counter()
            row['rdkit'] = [float(getattr(Descriptors, name)(heavy)) for name in TWO_D] + fpgen.GetFingerprintAsNumPy(heavy).astype(float).tolist()
            row['rdkit_seconds'] = time.perf_counter()-t
            t = time.perf_counter()
            try:
                result = XTBCalculator(**options.to_dict()).calculate(mol)
            except XTBEnvironmentError: raise
            except XTBRunError as exc:
                failures.append({**metadata, 'reason': str(exc)}); continue
            row['xtb_seconds'] = time.perf_counter()-t
            t = time.perf_counter(); report = CalcDescriptorReport(mol, result)
            row['analysis_seconds'] = time.perf_counter()-t
            if report['quality']['status'] == 'error':
                failures.append({**metadata, 'reason': report['quality']}); continue
            row['basic'] = [report['values'][k] for k in BASIC]
            row['orbital'] = [report['values'][k] for k in ORBITAL]
            if any(v is None or not np.isfinite(v) for v in row['basic'] + row['orbital']):
                failures.append({**metadata, 'reason': 'missing_or_nonfinite_xtb_features'}); continue
            row['quality'] = report['quality']; row['provenance'] = report['provenance']
            handle.write(json.dumps(row, allow_nan=False)+'\n'); handle.flush(); rows.append(row)
            if (n+1) % 100 == 0: print(f'calculated {n+1}/{len(selected)}', flush=True)
    write_json(args.output/'failures.json', failures)
    summary = {'scope': 'pilot, source geometries; no claim of excited-state xTB or full QM8 accuracy',
        'n_selected': len(selected), 'n_common_complete': len(rows), 'n_failed': len(failures),
        'failure_policy': 'same complete-case cohort for all groups, no replacement; exclusions retained',
        'featurization_wall_seconds': time.perf_counter()-started,
        'cost_seconds': {k: sum(r[k] for r in rows) for k in ('rdkit_seconds', 'xtb_seconds', 'analysis_seconds')},
        'source_exclusion_counts': dict(Counter(r['reason'] for r in excluded)),
        **fit_ablation(rows, args.output)}
    write_json(args.output/'summary.json', summary)


if __name__ == '__main__': main()
