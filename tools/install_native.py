#!/usr/bin/env python3
"""Install the standard Linux profile: matched RDKit, extended libxtb and rdXTBTools.

Requires conda and network access. Creates dedicated locked runtime/toolchain
prefixes, builds and verifies the engine, builds wheels from a clean sdist,
installs the compiled wheel and checks all methods outside the checkout.
No package upload is performed. Keep both prefixes at their installed paths.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run(args, **kwargs):
    print('+ ' + shlex.join(map(str, args)), flush=True)
    return subprocess.run(list(map(str, args)), check=True, **kwargs)


def verify_environment(prefix, lock):
    expected = dict(line.strip().split('#', 1) for line in lock.read_text().splitlines()
                    if line.startswith('https://'))
    records = [json.loads(p.read_text()) for p in (prefix / 'conda-meta').glob('*.json')]
    actual = {r.get('url'): r.get('sha256') for r in records}
    if not expected or actual != expected:
        raise ValueError(f'{prefix} does not match {lock.name}; choose a fresh --directory. '
                         'Existing environments are not solved, reset or overwritten.')


def environment(conda, prefix, lock, env):
    if not (prefix / 'conda-meta/history').exists():
        if prefix.exists() and any(prefix.iterdir()):
            raise ValueError(f'Refusing nonempty environment directory: {prefix}')
        run([conda, 'create', '--yes', '--prefix', prefix, '--file', lock], env=env)
    verify_environment(prefix, lock)


def clean_environment():
    env = os.environ.copy()
    for name in ('XTB_LIBRARY', 'XTBPATH', 'XTBHOME', 'XTB_PREFIX', 'LD_LIBRARY_PATH',
                 'LD_PRELOAD', 'PYTHONPATH', 'PYTHONHOME', 'RDKIT_ETH_BUILD_RDXTB',
                 'RDXTB_XTB_PREFIX', 'RDXTB_RDKIT_PREFIX'):
        env.pop(name, None)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    return env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--directory', type=Path, default=ROOT / '.rdkit-eth',
                        help='dedicated installation root (default: .rdkit-eth in checkout)')
    parser.add_argument('--jobs', type=int, default=min(os.cpu_count() or 1, 6))
    args = parser.parse_args()
    if platform.system() != 'Linux' or platform.machine() not in ('x86_64', 'amd64'):
        parser.error('The locked profile currently supports Linux x86-64 only.')
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    conda = shutil.which('conda')
    if not conda:
        parser.error('Install conda first and make its executable available on PATH.')
    directory = args.directory.resolve()
    runtime, toolchain = directory / 'runtime', directory / 'toolchain'
    env = clean_environment()
    for prefix, name in ((runtime, 'rdxtb-sdk.txt'), (toolchain, 'xtb-toolchain.txt')):
        environment(conda, prefix, ROOT / 'locks/linux-64' / name, env)
    run([conda, 'run', '--no-capture-output', '-p', toolchain, 'python',
         ROOT / 'tools/build_patched_xtb.py', runtime, '--dependency-prefix', toolchain,
         '--jobs', args.jobs], env=env)
    # The engine shares the Python prefix, so neither XTB_LIBRARY nor XTBPATH
    # is necessary. Explicit selection is only a build-time requirement.
    run([runtime / 'bin/python', '-m', 'pip', 'install', '--no-deps', '--require-hashes',
         '-r', ROOT / 'locks/linux-64/release-tools-py311.txt'], env=env)
    artifacts = directory / 'artifacts'
    artifacts.mkdir(exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix='build-', dir=artifacts))
    run([conda, 'run', '--no-capture-output', '-p', runtime, 'python',
         ROOT / 'tools/build_release.py', '--compiled', '--xtb-prefix', runtime,
         '--output', output], env=env)
    wheel, = (output / 'compiled').glob('*.whl')
    run([runtime / 'bin/python', '-m', 'pip', 'install', '--no-deps', '--force-reinstall', wheel], env=env)
    run([runtime / 'bin/python', '-m', 'pip', 'check'], env=env)
    reports = []
    with tempfile.TemporaryDirectory(prefix='rdkit-eth-installed-check-') as neutral:
        for method in ('GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'):
            for backend in ('native', 'compiled'):
                checked = run([runtime / 'bin/python', '-I', '-m', 'rdkit_eth', 'doctor',
                               '--check', backend, '--method', method, '--keep-overlap', '--json',
                               *(['--reactivity'] if method != 'GFN0-xTB' else [])],
                              env=env, cwd=neutral, capture_output=True, text=True)
                report = json.loads(checked.stdout)
                if not Path(report['package_path']).is_relative_to(runtime):
                    raise RuntimeError('Check imported a package outside the installed runtime')
                engine = report['check']['engine']
                selected = engine.get('library', engine.get('path'))
                if Path(selected).resolve() != (runtime / 'lib/libxtb.so').resolve():
                    raise RuntimeError('Check selected an engine outside the standard installation')
                reports.append(report)
                print(f'Installed {backend} {method} basis/populations: OK', flush=True)
    optimization_output = Path(tempfile.mkdtemp(prefix='optimization-', dir=artifacts))
    optimization = run([runtime / 'bin/python', '-I', ROOT / 'tools/benchmark_xtb_optimization.py',
                         '--output', optimization_output, '--smiles', 'O', '--smiles', 'CCO', '--check'],
                        env=env, cwd=directory, capture_output=True, text=True)
    optimization_report = json.loads(optimization.stdout)
    print('Installed in-process xTB optimization comparison: OK', flush=True)
    record = {'schema_version': 1, 'installed_at': datetime.now(timezone.utc).isoformat(),
              'runtime': str(runtime), 'toolchain': str(toolchain), 'wheel': str(wheel),
              'wheel_sha256': digest(wheel), 'source_lock_sha256': digest(ROOT / 'tools/xtb-source.json'),
              'checks': reports, 'optimization_check': optimization_report,
              'scope': 'local Linux x86-64; preserve both environment prefixes'}
    (directory / 'installation.json').write_text(json.dumps(record, indent=2) + '\n')
    print(f'\nVerified installation: {runtime}\nActivate: conda activate {shlex.quote(str(runtime))}\n'
          f'Report: {directory / "installation.json"}')


if __name__ == '__main__':
    main()
