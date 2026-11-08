#!/usr/bin/env python
"""Run the separate EAS model benchmark on a feature archive."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmarks.eas.pipeline import load_feature_archive
from benchmarks.eas.train import run_benchmark


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('archive', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--groups', default='BCDEFG')
    parser.add_argument('--models', default='logistic_regression,gradient_boosting')
    parser.add_argument('--bootstrap-replicates', type=int, default=2000)
    args = parser.parse_args(argv)
    archive = load_feature_archive(args.archive)
    result = run_benchmark(
        archive, groups=tuple(args.groups),
        model_names=tuple(name for name in args.models.split(',') if name),
        bootstrap_replicates=args.bootstrap_replicates)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'n_results': len(result['results']), 'output': str(args.output)}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
