#!/usr/bin/env python3
"""Package a verified runtime and dependency prefix for offline installation.

Build-only dependencies: conda-pack==0.8.1, setuptools==80.9.0 and patchelf==0.17.2.4.
This produces a Linux x86-64 environment bundle, not a manylinux wheel.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import tempfile

import conda_pack

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with open(path, 'rb') as stream:
        h = hashlib.sha256()
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--runtime', type=Path, required=True)
    parser.add_argument('--toolchain', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if platform.system() != 'Linux' or platform.machine() != 'x86_64': raise ValueError('Linux x86-64 required')
    if args.output.exists(): raise ValueError('choose a new output directory')
    patcher = shutil.which('patchelf')
    if not patcher: raise ValueError('patchelf required')
    runtime = args.runtime.resolve(); toolchain = args.toolchain.resolve()
    record = json.loads((runtime/'share/xtb/rdkit-eth/build.json').read_text())
    engine_files = [runtime/'lib/libxtb.so', runtime/'bin/xtb',
                    *(runtime/'share/xtb'/f'param_gfn{i}-xtb.txt' for i in range(3))]
    if any(not path.is_file() for path in engine_files):
        raise ValueError('runtime is missing an installed extended xTB artifact')
    engine_hashes = {str(path.relative_to(runtime)): digest(path) for path in engine_files}
    if any(record['files'].get(relative) != checksum for relative, checksum in engine_hashes.items()):
        raise ValueError('stale build manifest')
    library_hash = engine_hashes['lib/libxtb.so']
    # Record proof of a working installed profile before packing.
    env = __import__('os').environ.copy()
    for key in ('PYTHONPATH','PYTHONHOME','XTB_LIBRARY','XTBPATH','XTB_PREFIX','LD_LIBRARY_PATH','LD_PRELOAD'):
        env.pop(key, None)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1')
    subprocess.run([runtime/'bin/python', '-I', '-m', 'rdkit_eth', 'doctor',
                    '--check', 'compiled', '--keep-overlap'], env=env, check=True)
    args.output.mkdir(parents=True)
    # Archives can be large. Keep conda-pack's staging on the output filesystem,
    # rather than consuming a separate /tmp memory-filesystem quota.
    tempfile.tempdir = str(args.output.resolve())
    rpaths = {}
    for label, prefix in [('runtime', runtime), ('toolchain', toolchain)]:
        environment = conda_pack.CondaEnv.from_prefix(str(prefix))
        environment = environment.exclude('source/*').exclude('build/*').exclude('*.pyc')
        for item in environment.files:
            path = Path(item.source)
            if path.is_symlink() or not path.is_file(): continue
            with path.open('rb') as stream:
                if stream.read(4) != b'\x7fELF': continue
            result = subprocess.run([patcher, '--print-rpath', path], capture_output=True, text=True)
            if result.returncode == 0 and result.stdout.strip():
                rpaths[label+'/'+item.target] = result.stdout.strip()
        environment.pack(output=str(args.output/(label+'.tar.gz')), compress_level=2, n_threads=4)
    shutil.copy2(patcher, args.output/'patchelf')
    shutil.copy2(ROOT/'tools/install_bundle.py', args.output/'install.py')
    shutil.copy2(ROOT/'LICENSE', args.output/'LICENSE')
    # Supply the exact patched xTB source and licenses alongside the binary.
    source = Path(record['source_directory'])
    if not source.is_dir(): raise ValueError('exact engine source is required for distribution')
    with tarfile.open(args.output/'xtb-source.tar.gz', 'w:gz') as archive:
        for path in sorted(source.rglob('*')):
            if '.git' not in path.relative_to(source).parts and path.is_file():
                archive.add(path, arcname='xtb/'+str(path.relative_to(source)))
    manifest = {'schema_version': 1, 'platform': 'linux-x86_64',
        'prefixes': {'runtime': str(runtime), 'toolchain': str(toolchain)},
        'original_library_sha256': library_hash,
        'extended_xtb': {'files': engine_hashes, 'build_manifest_sha256': digest(runtime/'share/xtb/rdkit-eth/build.json')},
        'source_lock': record['source'],
        'rpaths': rpaths, 'build_host': platform.platform(),
        'scope': 'offline conda environment bundle; target glibc/CPU compatibility requires validation',
        'sha256': {p.name: digest(p) for p in sorted(args.output.iterdir()) if p.is_file()}}
    (args.output/'bundle.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print('Bundle:', args.output)


if __name__ == '__main__': main()
