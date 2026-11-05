#!/usr/bin/env python3
"""Link the pinned bounds builder into a diagnostic executable, with an optional patch.

Always build an unmodified bounds control as well as the candidate. The control
must match the installed library before interpreting candidate differences.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--base-build', type=Path, required=True)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--patch', type=Path)
    args = parser.parse_args()
    base = json.loads((args.base_build/'build-manifest.json').read_text())
    revision = subprocess.check_output(['git', '-C', str(args.source), 'rev-parse', 'HEAD'], text=True).strip()
    if revision != base['rdkit_source_commit']:
        raise ValueError('source revision differs from probe build')
    relative = 'Code/GraphMol/DistGeomHelpers/BoundsMatrixBuilder.cpp'
    upstream = (args.source/relative).resolve()
    if subprocess.check_output(['git', '-C', str(args.source), 'diff', 'HEAD', '--', relative], text=True):
        raise ValueError('bounds source must be unmodified')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError('output directory must be empty')
    bounds = output/'BoundsMatrixBuilder.cpp'
    if args.patch:
        subprocess.run(['patch', '--batch', '--fuzz=0', '-o', str(bounds),
                        str(upstream), str(args.patch.resolve())], check=True)
    else:
        bounds.write_bytes(upstream.read_bytes())
    command = list(base['command'])
    # Find by exact generated-source name rather than relying on SDK flags.
    candidates = [Path(arg) for arg in command if arg.endswith('/Embedder.instrumented.cpp')]
    if len(candidates) != 1 or digest(candidates[0]) != base['generated_source_sha256']:
        raise ValueError('base instrumentation differs from its manifest')
    if digest(Path(__file__).with_name('probe_driver.cpp')) != base['driver_sha256']:
        raise ValueError('probe driver differs from its manifest')
    command.insert(command.index(str(candidates[0])) + 1, str(bounds))
    command[command.index('-o')+1] = str(output/'embedding_probe')
    subprocess.run(command, check=True)
    manifest = {'protocol': 'local-ring-angle-build-v1', 'base_build': base,
        'rdkit_source_commit': revision, 'upstream_bounds_sha256': digest(upstream),
        'bounds_source_sha256': digest(bounds), 'patch_sha256': digest(args.patch) if args.patch else None,
        'builder_sha256': digest(Path(__file__)), 'command': command,
        'binary_sha256': digest(output/'embedding_probe')}
    (output/'build-manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')


if __name__ == '__main__':
    main()
