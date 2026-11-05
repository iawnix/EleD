#!/usr/bin/env python3
"""Build a diagnostic local-planarity ablation on an existing bounds variant."""
import argparse
import json
from pathlib import Path
import subprocess

from build_bounds_variant import digest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bounds-build', type=Path, required=True)
    parser.add_argument('--patch', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    base = json.loads((args.bounds_build/'build-manifest.json').read_text())
    command = list(base['command'])
    original = next(Path(arg) for arg in command if arg.endswith('/Embedder.instrumented.cpp'))
    bounds = next(Path(arg) for arg in command if arg.endswith('/BoundsMatrixBuilder.cpp'))
    if digest(original) != base['base_build']['generated_source_sha256'] or digest(bounds) != base['bounds_source_sha256']:
        raise ValueError('source differs from recorded bounds build')
    if digest(Path(__file__).with_name('probe_driver.cpp')) != base['base_build']['driver_sha256']:
        raise ValueError('driver differs from original build')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('output directory must be empty')
    changed = output/'Embedder.instrumented.cpp'
    subprocess.run(['patch', '--batch', '--fuzz=0', '-o', str(changed),
                    str(original), str(args.patch.resolve())], check=True)
    command[command.index(str(original))] = str(changed)
    command[command.index('-o')+1] = str(output/'embedding_probe')
    subprocess.run(command, check=True)
    manifest = {'protocol': 'ring-n-planarity-ablation-build-v1', 'bounds_build': base,
                'patch_sha256': digest(args.patch), 'builder_sha256': digest(Path(__file__)),
                'generated_source_sha256': digest(changed), 'command': command,
                'binary_sha256': digest(output/'embedding_probe')}
    (output/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
