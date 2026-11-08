#!/usr/bin/env python3
"""Refit the frozen QM8 pilot after historical stereo audit or checked preparation.

These are separate, predeclared protocols. Original evidence is never modified.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import time

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from eled import PrepareConformers, XTBCalculator, CalcDescriptorReport, XTBRunError, XTBEnvironmentError
from eled.geometry import _check_geometry, _stereo_constraints
from eled.input import _geometry_hash, _xyz_from_mol
from eled.persistence import _implementation_identity

_spec = importlib.util.spec_from_file_location('qm8_geometry', Path(__file__).with_name('benchmark_qm8_geometry.py'))
study = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(study)
pilot = study.pilot
ROOT = Path(__file__).resolve().parents[1]
PREPARATION = dict(num_conformers=8, method='ETKDGv3', fallback_methods=['ETKDGv4', 'srETKDGv3'],
    attempts_per_strategy=2, retry_random_coords=True, max_embedding_iterations=300,
    optimize='MMFF94', optimization_iterations=500, require_assigned_stereo=False)


def read_checked(root, name, *, lines=False):
    expected = json.loads((root/'SHA256SUMS.json').read_text())[name]
    if hashlib.sha256((root/name).read_bytes()).hexdigest() != expected:
        raise ValueError('evidence checksum mismatch: ' + name)
    opener = gzip.open if name.endswith('.gz') else open
    with opener(root/name, 'rt') as handle:
        value = [json.loads(line) for line in handle] if lines else json.load(handle)
    return value, expected


def audit_record(record):
    mol = Chem.AddHs(Chem.MolFromSmiles(record['smiles']))
    if [a.GetAtomicNum() for a in mol.GetAtoms()] != record['atomic_numbers']:
        raise ValueError('stored atom ordering mismatch')
    if len(record['positions_angstrom']) != mol.GetNumAtoms():
        raise ValueError('coordinate count mismatch')
    conf = Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
    for i, xyz in enumerate(record['positions_angstrom']): conf.SetAtomPosition(i, xyz)
    mol.AddConformer(conf)
    return mol, _check_geometry(mol, 0, _stereo_constraints(mol))


def restrict_splits(rows, splits):
    common = {r['id'] for r in rows}
    fixed = [{**s, 'train': [i for i in s['train'] if i in common],
              'test': [i for i in s['test'] if i in common]} for s in splits]
    pilot.validate_splits(rows, fixed)
    return fixed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', required=True, choices=['historical', 'checked'])
    parser.add_argument('--reference', type=Path, default=ROOT/'docs/validation/0.9.0')
    parser.add_argument('--historical', type=Path, default=ROOT/'docs/validation/0.9.1')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()): raise ValueError('output must be empty')
    frozen, source, splits, hashes = study.read_evidence(args.reference)
    protocol = 'qm8-stereo-refit-v1' if args.mode == 'historical' else 'qm8-checked-preparation-v1'
    history, history_hash = read_checked(args.historical, 'geometry-preparation.jsonl.gz', lines=True)
    old_rows, feature_hash = read_checked(args.historical, 'generated-features.jsonl.gz', lines=True)
    history = {r['id']: r for r in history}; old_rows = {r['id']: r for r in old_rows}
    if set(history) != {r['id'] for r in source}: raise ValueError('historical cohort mismatch')
    manifest = {'protocol': protocol, 'mode': args.mode, 'reference_hashes': hashes,
        'historical_hashes': {'geometry-preparation.jsonl.gz': history_hash,
                              'generated-features.jsonl.gz': feature_hash},
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'helper_sha256': {Path(p).name: hashlib.sha256(Path(p).read_bytes()).hexdigest()
                          for p in (study.__file__, pilot.__file__)},
        'implementation': _implementation_identity(), 'selected_ids': [r['id'] for r in source],
        'target': frozen['target'], 'options': frozen['options'],
        'preparation': PREPARATION if args.mode == 'checked' else 'reuse frozen historical coordinates; independently audit',
        'seed': 'SHA256(protocol:id) first 8 hex digits modulo 2147483646 plus 1',
        'selection': 'lowest MMFF energy among accepted conformers, conformer ID breaks ties',
        'split_policy': 'retain original memberships, refit all models on identical complete cases',
        'model': 'original fixed ExtraTrees and feature sets; no tuning',
        'bootstrap_replicates': study.BOOTSTRAP_REPLICATES,
        'scope': 'same QM8 pilot; not independent validation; source geometries never seed generation'}
    pilot.write_json(args.output/'manifest.json', manifest)
    rows, failures = [], []
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=False)
    started = time.perf_counter()
    with gzip.open(args.output/'preparation.jsonl.gz', 'wt') as prep, gzip.open(args.output/'features.jsonl.gz', 'wt') as out:
        for n, original in enumerate(source):
            identifier = original['id']; t = time.perf_counter()
            if args.mode == 'historical':
                previous = history[identifier]
                if previous['failure']:
                    failures.append({'id': identifier, 'stage': 'historical_geometry', 'reason': previous['failure']})
                    continue
                mol, check = audit_record(previous)
                info = {'id': identifier, 'check': check, 'reused_geometry': True}
                if check['ok']:
                    row = old_rows[identifier]
                    if row['geometry_sha256'] != _geometry_hash(_xyz_from_mol(mol, -1)):
                        raise ValueError('historical feature/geometry mismatch')
                else:
                    failures.append({'id': identifier, 'stage': 'coordinate_audit', 'reason': check['reason']})
                    mol = None
            else:
                seed = int(hashlib.sha256((protocol+':'+identifier).encode()).hexdigest()[:8], 16) % 2147483646 + 1
                prepared = PrepareConformers(original['smiles'], random_seed=seed, **PREPARATION)
                info = {'id': identifier, **prepared.report}
                mol = prepared.molecule if prepared.success else None
                if mol is None:
                    failures.append({'id': identifier, 'stage': 'checked_preparation', 'reason': info['error']['reason']})
                else:
                    chosen = info['recommended_conformer_id']
                    conf = Chem.Conformer(mol.GetConformer(chosen))
                    mol.RemoveAllConformers(); mol.AddConformer(conf, assignId=True)
                    info.update(atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                                positions_angstrom=mol.GetConformer().GetPositions().tolist())
            info['elapsed_seconds'] = time.perf_counter()-t
            prep.write(json.dumps(info, allow_nan=False)+'\n'); prep.flush()
            if mol is None: continue
            if args.mode == 'checked':
                heavy = Chem.RemoveHs(mol)
                two_d = [float(getattr(Descriptors, name)(heavy)) for name in pilot.TWO_D]
                two_d += fpgen.GetFingerprintAsNumPy(heavy).astype(float).tolist()
                if not np.allclose(two_d, original['rdkit'], atol=1e-10, rtol=1e-12):
                    raise ValueError('2D reconstruction mismatch')
                try: result = XTBCalculator(**frozen['options']).calculate(mol)
                except XTBEnvironmentError: raise
                except XTBRunError as exc:
                    failures.append({'id': identifier, 'stage': 'xtb', 'reason': str(exc)}); continue
                if (result.metadata['engine']['sha256'] != original['provenance']['engine']['sha256'] or
                    (result.metadata['parameter_file'] or {}).get('sha256') !=
                    (original['provenance']['parameter_file'] or {}).get('sha256')):
                    raise ValueError('engine identity differs from frozen source features')
                report = CalcDescriptorReport(mol, result)
                basic = [report['values'][k] for k in pilot.BASIC]
                orbital = [report['values'][k] for k in pilot.ORBITAL]
                if report['quality']['status'] == 'error' or any(v is None or not np.isfinite(v) for v in basic+orbital):
                    failures.append({'id': identifier, 'stage': 'features', 'reason': 'quality_or_missing_features'}); continue
                row = {k: original[k] for k in ('id', 'smiles', 'identity', 'scaffold', 'target_eV', 'rdkit')}
                row.update(basic=basic, orbital=orbital, provenance=report['provenance'], quality=report['quality'],
                           geometry_sha256=_geometry_hash(_xyz_from_mol(mol, -1)))
            for key in ('id', 'smiles', 'identity', 'scaffold', 'target_eV', 'rdkit'):
                if row[key] != original[key]: raise ValueError('feature identity mismatch: '+key)
            out.write(json.dumps(row, allow_nan=False)+'\n'); out.flush(); rows.append(row)
            if (n+1) % 100 == 0: print(f'{args.mode}: {n+1}/{len(source)}', flush=True)
    pilot.write_json(args.output/'failures.json', failures)
    common = {r['id'] for r in rows}
    fixed = restrict_splits(rows, splits)
    summary = {'protocol': protocol, 'n_selected': len(source), 'n_common': len(rows),
               'n_failed': len(failures), 'failure_counts': dict(Counter(r['reason'] for r in failures)),
               'preparation_and_features_seconds': time.perf_counter()-started}
    predictions = {}
    for label, features in [('source', [r for r in source if r['id'] in common]), ('generated', rows)]:
        folder = args.output/label; folder.mkdir()
        summary[label] = pilot.fit_ablation(features, folder, splits=fixed)
        predictions[label] = json.loads((folder/'predictions.json').read_text())
    pilot.write_json(args.output/'paired-comparisons.json', study.compare_predictions(rows, predictions))
    pilot.write_json(args.output/'summary.json', summary)
    files = [p for p in args.output.rglob('*') if p.is_file()]
    pilot.write_json(args.output/'SHA256SUMS.json', {
        str(p.relative_to(args.output)): hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(files)})


if __name__ == '__main__': main()
