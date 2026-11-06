#!/usr/bin/env python3
"""Install a checksum-verified Linux runtime bundle into a new directory, offline."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import tarfile
import tempfile


def digest(path):
    with open(path, 'rb') as stream:
        h = hashlib.sha256()
        for block in iter(lambda: stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def run(command, **kwargs):
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def install(bundle, destination):
    if platform.system() != 'Linux' or platform.machine() != 'x86_64':
        raise ValueError('this bundle supports Linux x86-64 only')
    destination = destination.resolve()
    if destination.exists(): raise ValueError('destination must not exist; do not overwrite an installation')
    manifest = json.loads((bundle/'bundle.json').read_text())
    if manifest['schema_version'] != 1: raise ValueError('unsupported bundle schema')
    for name, expected in manifest['sha256'].items():
        path = (bundle/name).resolve()
        if not path.is_relative_to(bundle.resolve()) or digest(path) != expected:
            raise ValueError('bundle checksum mismatch: '+name)
    env = os.environ.copy()
    for name in ('PYTHONPATH', 'PYTHONHOME', 'XTB_LIBRARY', 'XTBPATH', 'XTB_PREFIX',
                 'CONDA_PREFIX', 'LD_LIBRARY_PATH', 'LD_PRELOAD'):
        env.pop(name, None)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    destination.mkdir(parents=True)
    # A failed install is left in place for inspection; never recursively remove
    # a caller-selected directory. Restart with a fresh destination.
    for label in ('toolchain', 'runtime'):
        target = destination/label; target.mkdir()
        with tarfile.open(bundle/(label+'.tar.gz')) as archive:
            archive.extractall(target, filter='data')
    runtime = destination/'runtime'
    library = (runtime/'lib/libxtb.so').resolve()
    if digest(library) != manifest['original_library_sha256']:
        raise ValueError('unpacked libxtb differs from the verified build')
    engine = manifest.get('extended_xtb', {})
    for relative, expected in engine.get('files', {}).items():
        artifact = runtime / relative
        if not artifact.is_file() or digest(artifact) != expected:
            raise ValueError('unpacked extended xTB artifact differs from the verified build: ' + relative)
    old_record = runtime/'share/xtb/rdkit-eth/build.json'
    old_digest = digest(old_record)
    for label in ('toolchain', 'runtime'):
        prefix = destination/label
        run([runtime/'bin/python', prefix/'bin/conda-unpack'], env=env, cwd=destination)
    patcher = bundle/'patchelf'
    for relative, original in manifest['rpaths'].items():
        binary = destination/relative
        if not binary.is_file(): raise ValueError('missing ELF file: '+relative)
        # conda-unpack relocates each prefix; also repair cross-prefix references
        # and untracked native binaries compiled by this project.
        current = subprocess.check_output([patcher, '--print-rpath', binary], text=True).strip()
        repaired = current
        for label, old in manifest['prefixes'].items():
            repaired = repaired.replace(old, str(destination/label))
        if repaired != current: run([patcher, '--set-rpath', repaired, binary], env=env)
    record = json.loads(old_record.read_text())
    record['relocation'] = {'original_build_record_sha256': old_digest,
        'original_library_sha256': manifest['original_library_sha256'],
        'bundle_manifest_sha256': digest(bundle/'bundle.json'),
        'operation': 'conda-unpack and RPATH prefix relocation; numerical source unchanged'}
    record['prefix'] = str(runtime)
    record['dependency_prefix'] = str(destination/'toolchain')
    for name in record['files']:
        if (runtime/name).is_file(): record['files'][name] = digest(runtime/name)
    old_record.write_text(json.dumps(record, indent=2)+'\n')
    reports = []
    with tempfile.TemporaryDirectory(prefix='rdkit-eth-bundle-check-') as neutral:
        for method in ('GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'):
            for interface in ('native', 'compiled'):
                check = run([runtime/'bin/python', '-I', '-m', 'rdkit_eth', 'doctor', '--check', interface,
                    '--method', method, '--keep-overlap', '--json',
                    *(['--reactivity'] if method != 'GFN0-xTB' else [])],
                    env=env, cwd=neutral, capture_output=True, text=True)
                report = json.loads(check.stdout)
                engine = report['check']['engine']
                selected = engine.get('path', engine.get('library'))
                if Path(selected).resolve() != library: raise ValueError('check selected an external libxtb')
                if not Path(report['package_path']).is_relative_to(runtime):
                    raise ValueError('check imported an external package')
                reports.append(report)
    dependencies = {}
    candidates = [library, *runtime.glob('lib/python*/site-packages/rdkit_eth/*.so')]
    for binary in candidates:
        listing = subprocess.check_output(['ldd', binary], text=True, env=env)
        if 'not found' in listing or any(old in listing for old in manifest['prefixes'].values()):
            raise ValueError('unresolved/original-prefix ELF dependency: '+str(binary))
        dependencies[str(binary.relative_to(destination))] = listing
    (destination/'installation.json').write_text(json.dumps({
        'schema_version': 1, 'bundle_manifest_sha256': digest(bundle/'bundle.json'),
        'runtime': str(runtime), 'checks': reports, 'ldd': dependencies}, indent=2)+'\n')
    shutil.copy2(bundle/'bundle.json', destination/'bundle.json')
    print('Installed and verified:', runtime)
    print('Use:', runtime/'bin/python')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bundle', type=Path, default=Path(__file__).resolve().parent)
    parser.add_argument('--directory', type=Path, required=True)
    args = parser.parse_args()
    install(args.bundle.resolve(), args.directory)
