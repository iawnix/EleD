#!/usr/bin/env python
"""Materialize one EAS atom feature group for the separate benchmark layer."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmarks.eas.pipeline import feature_matrix, load_feature_archive


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('archive', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--group', choices=tuple('ABCDEFG'), default='A')
    args = parser.parse_args(argv)
    archive = load_feature_archive(args.archive)
    rows = feature_matrix(archive, args.group)
    payload = {'protocol_id': archive['protocol_id'],
               'protocol_sha256': archive['protocol_sha256'],
               'model_feature_group': args.group, 'n_rows': len(rows),
               'rows': rows}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(payload, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print(json.dumps({'group': args.group, 'n_rows': len(rows)}, indent=2))
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
