#!/usr/bin/env python
"""Validate and split a local atom-labelled EAS export.

No network access, conformers or xTB calculations are used here.  The output
is the frozen input manifest consumed by a later feature-generation stage.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from rdkit_eth.eas import prepare_eas_dataset

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--source', default='local')
    parser.add_argument('--seed', type=int, default=None)
    args = parser.parse_args(argv)
    report = prepare_eas_dataset(args.input, source=args.source, seed=args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ('protocol_id','protocol_sha256','n_input','n_selected','n_excluded','split_counts')}, indent=2))
    return 0

if __name__ == '__main__': raise SystemExit(main())
