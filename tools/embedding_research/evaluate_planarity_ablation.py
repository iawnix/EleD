#!/usr/bin/env python3
"""Complete a two-factor local-angle / local-planarity diagnostic ablation.

Reuses the frozen angle-only experiment; generates no source-derived seeds.
Removing a targeted improper is an explanatory ablation, not a validated model.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import json
from pathlib import Path

from evaluate_ring_angle import SEEDS, PRESETS, MODES, assess, digest, invoke, tally


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--initial', type=Path, required=True)
    parser.add_argument('--planarity-only', type=Path, required=True)
    parser.add_argument('--combined', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise ValueError('output directory must be empty')
    initial_manifest = json.loads((args.initial/'manifest.json').read_text())
    if initial_manifest['script_sha256'] != digest(Path(__file__).with_name('evaluate_ring_angle.py')):
        raise ValueError('initial experiment code changed')
    for name, sha in initial_manifest['dependencies'].items():
        if digest(Path(__file__).with_name(name)) != sha:
            raise ValueError('initial dependency changed: '+name)
    binaries = {'planarity_only': args.planarity_only, 'combined': args.combined}
    builds = {name: json.loads(binary.with_name('build-manifest.json').read_text()) for name, binary in binaries.items()}
    for name, binary in binaries.items():
        if builds[name]['binary_sha256'] != digest(binary):
            raise ValueError('ablation binary changed')
    if builds['planarity_only']['bounds_build'] != initial_manifest['builds']['control']:
        raise ValueError('planarity-only needs the original bounds control')
    if builds['combined']['bounds_build'] != initial_manifest['builds']['candidate']:
        raise ValueError('combined needs the angle-only candidate bounds')
    manifest = {'protocol': 'local-ring-planarity-ablation-v1', 'builds': builds,
        'initial_hashes': {name: digest(args.initial/name) for name in ('manifest.json', 'screen.jsonl.gz', 'sweep.jsonl.gz', 'summary.json')},
        'script_sha256': digest(Path(__file__)),
        'initial_script_sha256': digest(Path(__file__).with_name('evaluate_ring_angle.py')),
        'seeds': SEEDS, 'presets': PRESETS, 'modes': MODES, 'max_iterations': 100,
        'scope': 'two-factor diagnostic ablation on previously identified targets, not a new holdout',
        'source_coordinates_used': False}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    summary = {'screen': defaultdict(Counter), 'screen_gains': defaultdict(list),
               'screen_losses': defaultdict(list), 'unchanged_non_targets': Counter(),
               'sweep': defaultdict(Counter), 'paired_target_gains': Counter(), 'paired_target_losses': Counter()}
    prior_screen = [json.loads(line) for line in gzip.open(args.initial/'screen.jsonl.gz', 'rt')]
    targets = {r['id']: r for r in prior_screen if r['id'] in initial_manifest['target_ids']}
    with gzip.open(args.output/'screen.jsonl.gz', 'wt') as handle:
        for n, row in enumerate(prior_screen):
            record = {'id': row['id'], 'smiles': row['smiles'], 'reports': {}, 'seconds': {}, 'stderr': {}}
            original = row['reports']['baseline']
            old_ok = bool(original['coordinate_check'] and original['coordinate_check']['ok'])
            for name, binary in binaries.items():
                events, seconds, stderr = invoke(binary, row['smiles'], 'ETKDGv4', 'dg', 42, 'UFF')
                report = assess(row['smiles'], events)
                if row['id'] not in targets:
                    if report != original:
                        raise AssertionError('unexpected non-target change: '+row['id'])
                    summary['unchanged_non_targets'][name] += 1
                tally(summary['screen'][name], report)
                ok = bool(report['coordinate_check'] and report['coordinate_check']['ok'])
                if ok and not old_ok: summary['screen_gains'][name].append(row['id'])
                if old_ok and not ok: summary['screen_losses'][name].append(row['id'])
                record['reports'][name] = report
                record['seconds'][name] = seconds
                record['stderr'][name] = stderr
                if row['id'] in targets:
                    record.setdefault('events', {})[name] = events
            handle.write(json.dumps(record, allow_nan=False)+'\n')
            handle.flush()
            if (n+1) % 200 == 0: print('ablation screen', n+1, flush=True)
    prior_sweep = {(r['id'], r['preset'], r['mode'], r['seed'], r['variant']): r
        for r in map(json.loads, gzip.open(args.initial/'sweep.jsonl.gz', 'rt'))}
    with gzip.open(args.output/'sweep.jsonl.gz', 'wt') as handle:
        for identifier, row in targets.items():
            for preset in PRESETS:
                for mode in MODES:
                    for seed in SEEDS:
                        old_report = prior_sweep[identifier, preset, mode, seed, 'baseline']['report']
                        old_ok = bool(old_report['coordinate_check'] and old_report['coordinate_check']['ok'])
                        for name, binary in binaries.items():
                            events, seconds, stderr = invoke(binary, row['smiles'], preset, mode, seed, 'UFF')
                            report = assess(row['smiles'], events, optimize=True)
                            key = '/'.join((identifier, preset, mode, name))
                            tally(summary['sweep'][key], report)
                            ok = bool(report['coordinate_check'] and report['coordinate_check']['ok'])
                            if ok and not old_ok: summary['paired_target_gains'][name] += 1
                            if old_ok and not ok: summary['paired_target_losses'][name] += 1
                            handle.write(json.dumps({'id': identifier, 'smiles': row['smiles'], 'preset': preset,
                                'mode': mode, 'seed': seed, 'variant': name, 'report': report,
                                'events': events, 'seconds': seconds, 'stderr': stderr}, allow_nan=False)+'\n')
                            handle.flush()
            print('ablation sweep', identifier, 'complete', flush=True)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k != 'sweep'}, indent=2), flush=True)


if __name__ == '__main__':
    main()
