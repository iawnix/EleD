#!/usr/bin/env python
"""Generate GFN2 density features for a prepared EAS manifest.

The manifest must come from ``prepare_eas_dataset.py``. Failures are retained
per molecule and are never converted into zero-valued atom features.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from rdkit_eth import XTBOptions, prepare_conformers, run_conformer_analysis

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--xtb-library', type=Path)
    parser.add_argument('--num-conformers', type=int, default=10)
    parser.add_argument('--seed', type=int, default=20261106)
    args = parser.parse_args(argv)
    manifest = json.loads(args.manifest.read_text(encoding='utf-8'))
    if manifest.get('protocol_id') != 'eas_regioselectivity_v1':
        raise SystemExit('manifest was not produced by the frozen EAS protocol')
    records = []
    for row in manifest.get('records', []):
        digest = hashlib.sha256(f"{args.seed}:{row['id']}".encode()).hexdigest()
        seed = int(digest[:8], 16) % 2147483646 + 1
        item = {'id': row['id'], 'reactive_atoms': row['reactive_atoms'], 'split': row['split']}
        try:
            prepared = prepare_conformers(row['smiles'], num_conformers=args.num_conformers,
                method='ETKDGv3', optimize='xTB', xtb_method='GFN2-xTB',
                xtb_library=args.xtb_library, random_seed=seed)
            item['preparation'] = prepared.report
            if not prepared.success:
                item['status'] = 'preparation_failed'; records.append(item); continue
            options = XTBOptions(method='GFN2-xTB', xtb_library=args.xtb_library, keep_overlap=True,
                                 electronic_temperature=300.0)
            analysis = run_conformer_analysis(prepared.molecule, options=options,
                record_id=row['id'], return_errors=True)
            item['analysis'] = analysis
            item['status'] = analysis['status']
        except Exception as exc:
            item['status'] = 'failed'
            item['error'] = {'type': type(exc).__name__, 'message': str(exc)}
        records.append(item)
    output = {'protocol_id': manifest['protocol_id'], 'feature_method': 'GFN2-xTB',
              'electronic_temperature_K': 300.0, 'n_records': len(records), 'records': records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'n_records': len(records),
                      'n_complete': sum(r.get('status') == 'complete' for r in records),
                      'n_failed': sum(r.get('status') not in ('complete','partial') for r in records)}, indent=2))
    return 0

if __name__ == '__main__': raise SystemExit(main())
