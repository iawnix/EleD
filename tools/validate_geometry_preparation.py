#!/usr/bin/env python3
"""Audit the geometry API on frozen QM8 failures plus deterministic controls."""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import time

from rdkit_eth import prepare_conformers
from rdkit_eth.persistence import _implementation_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()): raise ValueError('output must be empty')
    root = Path(__file__).resolve().parents[1]
    inputs = [root/'docs/validation/0.9.0/qm8-manifest.json.gz', root/'docs/validation/0.9.1/failures.json']
    selected = json.load(gzip.open(inputs[0], 'rt'))['selected']
    failed = {r['id'] for r in json.loads(inputs[1].read_text())}
    rows = [dict(r, group='prior_embedding_failure') for r in selected if r['id'] in failed]
    rows += [dict(r, group='prior_success_control') for r in selected if r['id'] not in failed][:40]
    manifest = {'protocol': 'geometry-preparation-validation-v1',
        'source_hashes': {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'implementation': _implementation_identity(), 'ids': [r['id'] for r in rows],
        'settings': {'num_conformers': 1, 'random_seed': 42},
        'scope': 'new preparation protocol; original prediction cohort/results are not changed'}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    summary = {'groups': {}, 'recovered_ids': [], 'remaining_ids': [], 'recovery_methods': {},
               'coordinate_stereo_rejections': 0, 'embedding_failure_counts': {}, 'elapsed_seconds': 0.}
    methods, causes = Counter(), Counter()
    with gzip.open(args.output/'reports.jsonl.gz', 'wt') as handle:
        for row in rows:
            start = time.perf_counter(); result = prepare_conformers(row['smiles'], num_conformers=1, random_seed=42)
            elapsed = time.perf_counter()-start; report = result.report
            group = summary['groups'].setdefault(row['group'], {'total': 0, 'accepted': 0, 'failed': 0})
            group['total'] += 1; group['accepted' if result.success else 'failed'] += 1
            summary['elapsed_seconds'] += elapsed
            if row['group'] == 'prior_embedding_failure':
                summary['recovered_ids' if result.success else 'remaining_ids'].append(row['id'])
                if result.success:
                    a = report['attempts'][report['conformers'][0]['attempt_index']]
                    methods[f"{a['method']}; random_coords={a['use_random_coords']}"] += 1
            for attempt in report['attempts']:
                causes.update(attempt['failure_counts']['named'])
                for c in attempt['candidates']:
                    summary['coordinate_stereo_rejections'] += any(c.get(k, {}).get('reason') == 'stereo_mismatch'
                        for k in ('before_optimization', 'after_optimization'))
            record = {'id': row['id'], 'group': row['group'], 'elapsed_seconds': elapsed, 'report': report,
                'atomic_numbers': [a.GetAtomicNum() for a in result.molecule.GetAtoms()],
                'positions_angstrom': {str(c.GetId()): c.GetPositions().tolist() for c in result.molecule.GetConformers()}}
            handle.write(json.dumps(record, allow_nan=False)+'\n'); handle.flush()
            print(row['id'], row['group'], report['status'], flush=True)
    summary['recovery_methods'] = dict(methods); summary['embedding_failure_counts'] = dict(causes)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
