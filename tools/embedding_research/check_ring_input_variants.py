#!/usr/bin/env python3
"""Supplement retention checks for MMFF bounds and unspecified-stereo SMILES.

Unspecified-stereo inputs are derived stress cases, not new molecule identities.
enforceChirality stays enabled in every executable. No source coordinates enter.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path

from rdkit import Chem

from check_ring_nonregression import measure
from evaluate_ring_angle import SEEDS, digest, invoke
from nonregression_gate import compare_pair, verdict


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cohort', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--planarity-only', type=Path, required=True)
    parser.add_argument('--combined', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()): raise ValueError('output must be empty')
    affected = [r for r in json.loads(args.cohort.read_text())['rows'] if r['group'] == 'affected']
    rows = []
    for row in affected:
        mol = Chem.MolFromSmiles(row['smiles'])
        Chem.RemoveStereochemistry(mol)
        unspecified = Chem.MolToSmiles(mol)
        for label, smiles, ff in (('specified_mmff', row['smiles'], 'MMFF'),
                                  ('unspecified_uff', unspecified, 'UFF'),
                                  ('unspecified_mmff', unspecified, 'MMFF')):
            rows.append({**row, 'smiles': smiles, 'source_smiles': row['smiles'], 'stratum': label, 'force_field': ff})
    (args.output/'inputs.json').write_text(json.dumps(rows, indent=2)+'\n')
    binaries = {'baseline': args.baseline, 'planarity_only': args.planarity_only, 'combined': args.combined}
    builds = {k: json.loads(v.with_name('build-manifest.json').read_text()) for k, v in binaries.items()}
    for k, v in binaries.items():
        if digest(v) != builds[k]['binary_sha256']: raise ValueError('binary identity changed')
    manifest = {'protocol': 'ring-input-retention-v1', 'builds': builds,
        'cohort_sha256': digest(args.cohort), 'inputs_sha256': digest(args.output/'inputs.json'),
        'preset': 'ETKDGv4', 'mode': 'dg', 'seeds': SEEDS, 'max_embedding_attempts': 100,
        'post_mmff_max_iterations': 500, 'energy_tolerance_kcal_mol': .001,
        'script_sha256': digest(Path(__file__)),
        'dependencies': {name: digest(Path(__file__).with_name(name)) for name in (
            'check_ring_nonregression.py', 'evaluate_ring_angle.py', 'nonregression_gate.py',
            'run_diagnosis.py', 'audit_constraints.py')},
        'scope': 'specified/MMFF and unspecified-stereo UFF/MMFF stress inputs on the same 40 molecular identities',
        'source_coordinates_used': False, 'execution_order': 'rotated per paired trial',
        'timing_limit': 'diagnostic process/logging wall times; concurrent experiment load; not a production speed test'}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    counts = defaultdict(Counter)
    comparisons = defaultdict(list)
    successes = Counter()
    failures = []
    trial = 0
    with gzip.open(args.output/'records.jsonl.gz', 'wt') as handle:
        for row in rows:
            for seed in SEEDS:
                labels = list(binaries)
                order = labels[trial % 3:] + labels[:trial % 3]
                reports, seconds, stderr, events_by_label = {}, {}, {}, {}
                for label in order:
                    events, seconds[label], stderr[label] = invoke(binaries[label], row['smiles'], 'ETKDGv4', 'dg', seed, row['force_field'])
                    reports[label] = measure(row['smiles'], events)
                    events_by_label[label] = events
                    counter = counts[row['stratum']+'/'+label]
                    counter['runs'] += 1
                    counter['initial_valid'] += reports[label]['initial_valid']
                    counter['post_valid'] += reports[label]['post_valid']
                original = reports['baseline']
                successes[row['stratum']] += original['initial_valid']
                record = {'id': row['id'], 'stratum': row['stratum'], 'smiles': row['smiles'],
                    'force_field': row['force_field'], 'seed': seed, 'execution_order': order,
                    'reports': reports, 'seconds': seconds, 'stderr': stderr, 'comparisons': {}}
                for label in ('planarity_only', 'combined'):
                    comp = compare_pair(original, reports[label])
                    record['comparisons'][label] = comp
                    comparisons[row['stratum']+'/'+label].append(comp)
                    if comp['failures']:
                        failures.append({'id': row['id'], 'stratum': row['stratum'], 'seed': seed, 'variant': label, **comp})
                        record.setdefault('regression_events', {})[label] = events_by_label[label]
                        record.setdefault('regression_events', {})['baseline'] = events_by_label['baseline']
                handle.write(json.dumps(record, allow_nan=False)+'\n')
                handle.flush()
                trial += 1
            if trial % 75 == 0: print('input-variant paired trials', trial, 'flagged comparisons', len(failures), flush=True)
    summary = {'counts': counts, 'paired_trials': trial, 'failures': failures,
        'verdicts': {k: verdict(v, successes[k.split('/')[0]]) for k, v in comparisons.items()}}
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary['verdicts'], indent=2), flush=True)


if __name__ == '__main__': main()
