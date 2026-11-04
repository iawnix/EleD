"""Build and verify the pinned xTB dependency; Python standard library only."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys

TOOLS = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def run(command, **kwargs):
    print('+ ' + shlex.join(map(str, command)), flush=True)
    return subprocess.run(list(map(str, command)), check=True, **kwargs)


def git_output(source, *args):
    return subprocess.check_output(['git', '-C', str(source), *args], text=True).strip()


def source_tree(prefix, supplied, lock, patch=None):
    source = supplied.resolve() if supplied else prefix / 'source'
    if not source.exists():
        if supplied:
            raise ValueError(f'Source directory does not exist: {source}')
        source.mkdir(parents=True)
        run(['git', 'init', source])
        run(['git', '-C', source, 'remote', 'add', 'origin', lock['repository']])
        run(['git', '-C', source, 'fetch', '--depth=1', 'origin', lock['commit']])
        run(['git', '-C', source, 'checkout', '-b', 'rdkit-eth/native-basis', 'FETCH_HEAD'])
    actual = git_output(source, 'rev-parse', 'HEAD')
    if actual != lock['commit']:
        raise ValueError(f'Source HEAD is {actual}, expected {lock["commit"]}. '
                         'Use a separate prefix or a checkout of the pinned commit.')
    changes = git_output(source, 'diff', 'HEAD', '--name-only').splitlines()
    untracked = git_output(source, 'ls-files', '--others', '--exclude-standard')
    extension = lock.get('basis_extension')
    allowed = ({lock['patched_file']} if patch else set()) | (set(extension['patched_files']) if extension else set())
    if set(changes) - allowed or untracked:
        raise ValueError('Source contains unrelated changes; use a separate checkout. '
                         'The builder never resets or discards source changes.')
    if patch:
        # Historical source locks remain reproducible; current builds use the upstream fix.
        check = subprocess.run(['git', '-C', str(source), 'apply', '--check', str(patch)],
                               capture_output=True, text=True)
        if check.returncode == 0:
            run(['git', '-C', source, 'apply', patch])
        if digest(source / lock['patched_file']) != lock['patched_file_sha256']:
            raise ValueError('Source does not match the reviewed D4 patch.')
    else:
        fix = lock['gfn0_fix']
        if (fix['kind'] != 'upstream' or fix['commit'] != lock['commit']
                or digest(source / fix['file']) != fix['file_sha256']):
            raise ValueError('Source does not match the reviewed upstream GFN0 repair.')
    if extension:
        extra_patch = (TOOLS / extension['patch']).resolve()
        if digest(extra_patch) != extension['patch_sha256']:
            raise ValueError('Basis extension patch hash differs from source lock')
        check = subprocess.run(['git', '-C', str(source), 'apply', '--check', str(extra_patch)],
                               capture_output=True, text=True)
        if check.returncode == 0:
            run(['git', '-C', source, 'apply', extra_patch])
        for filename, expected in extension['patched_files'].items():
            if digest(source / filename) != expected:
                raise ValueError('Source does not match reviewed basis extension: ' + filename)
    return source


def read_cache(build):
    values = {}
    for line in (build / 'CMakeCache.txt').read_text().splitlines():
        if not line or line.startswith(('#', '//')) or '=' not in line:
            continue
        key, value = line.split('=', 1)
        values[key.split(':', 1)[0]] = value
    return values


def verify_library(build, library_dir, parameter_dir, compiler, dependency_prefix, env):
    executable = build / 'verify-gfn0-state'
    command = [compiler, '-std=c11', '-Wall', '-Wextra', '-Werror',
               TOOLS / 'diagnostics/xtb_gfn0_state.c', '-I' + str(parameter_dir.parent.parent / 'include'),
               '-L' + str(library_dir), '-Wl,-rpath,' + str(library_dir)]
    if dependency_prefix:
        command.append('-Wl,-rpath-link,' + str(dependency_prefix / 'lib'))
    run([*command, '-lxtb', '-ldl', '-lm', '-o', executable])
    records = []
    for molecule in ('ethanol', 'water'):
        for method in (0, 1, 2):
            for action in ('load', 'singlepoint'):
                args = [executable, parameter_dir, str(method), '--expect-clean']
                if molecule == 'water':
                    args.append('--water')
                if action == 'singlepoint':
                    args.append('--run-intermediate')
                result = run(args, env=env, cwd=build, capture_output=True, text=True, timeout=60)
                values = [json.loads(line) for line in result.stdout.splitlines() if line.startswith('{')]
                if len(values) != (3 if action == 'singlepoint' else 2):
                    raise ValueError(f'Unexpected diagnostic output: {result.stdout}')
                for value in values:
                    if not math.isfinite(value['energy_Eh']):
                        raise ValueError(f'Nonfinite diagnostic energy: {value}')
                after = values[-1]
                for field, tolerance in [('delta_energy_Eh', 1e-10), ('max_orbital_delta_eV', 1e-8)]:
                    if not math.isfinite(after[field]) or abs(after[field]) > tolerance:
                        raise ValueError(f'Method-switch check failed: {after}')
                records.append(dict(molecule=molecule, method=method, action=action, values=values))
    return records


def verify_basis(build, prefix, compiler, dependency_prefix, env):
    executable = build / 'verify-basis-export'
    command = [compiler, '-std=c11', '-Wall', '-Wextra', '-Werror',
               TOOLS / 'diagnostics/xtb_basis_export.c', '-I' + str(prefix / 'include'),
               '-L' + str(prefix / 'lib'), '-Wl,-rpath,' + str(prefix / 'lib')]
    if dependency_prefix:
        command.append('-Wl,-rpath-link,' + str(dependency_prefix / 'lib'))
    run([*command, '-lxtb', '-lm', '-o', executable])
    result = run([executable, prefix / 'share/xtb'], env=env, cwd=build,
                 capture_output=True, text=True, timeout=60)
    records = [json.loads(line) for line in result.stdout.splitlines() if line.startswith("{")]
    if len(records) != 1:
        raise RuntimeError("basis diagnostic must return exactly one JSON record")
    return records[0]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prefix', nargs='?', type=Path, default=Path('.xtb-native'))
    parser.add_argument('source', nargs='?', type=Path, help='optional clean Git checkout of the pinned commit')
    parser.add_argument('--dependency-prefix', type=Path,
                        help='CMake dependency and runtime prefix, e.g. the compiler conda environment')
    parser.add_argument('--jobs', type=int, default=min(os.cpu_count() or 1, 8))
    parser.add_argument('--cmake-arg', action='append', default=[],
                        help='extra CMake option, e.g. --cmake-arg=-DCMAKE_C_COMPILER=/path/to/gcc')
    args = parser.parse_args()
    if sys.platform != 'linux':
        parser.error('This build workflow currently supports Linux; other platforms are not validated.')
    if args.jobs < 1:
        parser.error('--jobs must be positive')
    for command in ('git', 'cmake', 'ctest', 'ninja'):
        if not shutil.which(command):
            parser.error(f'{command} is required; see environment-xtb-build.yml')

    lock = json.loads((TOOLS / 'xtb-source.json').read_text())
    patch = (TOOLS / lock['patch']).resolve() if 'patch' in lock else None
    if patch and digest(patch) != lock['patch_sha256']:
        raise ValueError('Patch differs from tools/xtb-source.json; review and update the source lock.')
    prefix = args.prefix.resolve()
    dependency_prefix = args.dependency_prefix.resolve() if args.dependency_prefix else None
    source = source_tree(prefix, args.source, lock, patch)
    build = prefix / 'build'
    rpath = ['$ORIGIN', '$ORIGIN/../lib']
    configure = ['cmake', '-S', source, '-B', build, '-G', 'Ninja',
                 '-DCMAKE_BUILD_TYPE=Release', '-DCMAKE_INSTALL_PREFIX=' + str(prefix),
                 '-DCMAKE_INSTALL_LIBDIR=lib', '-DBUILD_SHARED_LIBS=ON', '-DWITH_OBJECT=ON',
                 '-DCMAKE_POSITION_INDEPENDENT_CODE=ON', '-DWITH_TESTS=ON',
                 '-DWITH_TBLITE=OFF', '-DWITH_CPCMX=OFF', '-DCMAKE_C_STANDARD=17',
                 '-DCMAKE_POLICY_VERSION_MINIMUM=3.5', '-DCMAKE_INSTALL_RPATH_USE_LINK_PATH=ON']
    if dependency_prefix:
        if not (dependency_prefix / 'lib').is_dir():
            raise ValueError(f'Dependency prefix has no lib directory: {dependency_prefix}')
        rpath.append(str(dependency_prefix / 'lib'))
        configure += ['-DCMAKE_PREFIX_PATH=' + str(dependency_prefix),
                      '-DCMAKE_EXE_LINKER_FLAGS=-Wl,-rpath-link,' + str(dependency_prefix / 'lib')]
    configure += ['-DCMAKE_INSTALL_RPATH=' + ';'.join(rpath), *args.cmake_arg]
    run(configure, cwd=source)
    cache = read_cache(build)
    # Keep the installation layout used by the verifier/activation file explicit.
    if (Path(cache['CMAKE_INSTALL_PREFIX']).resolve() != prefix
            or cache['CMAKE_INSTALL_LIBDIR'] != 'lib' or cache['WITH_TESTS'] != 'ON'
            or cache['WITH_TBLITE'] != 'OFF' or cache['WITH_CPCMX'] != 'OFF'):
        raise ValueError('Do not override prefix/lib layout, required tests, or the disabled tblite/CPCM-X dependencies.')
    run(['cmake', '--build', build, '--parallel', str(args.jobs)])
    env = os.environ.copy()
    for name in ('LD_LIBRARY_PATH', 'LD_PRELOAD', 'XTBPATH', 'XTBHOME'):
        env.pop(name, None)
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    # Tests must pass before installing or certifying the build. No skip-tests path.
    ctest = run(['ctest', '--test-dir', build, '--output-on-failure', '--no-tests=error',
                 '--parallel', str(args.jobs)], env=env, text=True, capture_output=True)
    print(ctest.stdout, end='', flush=True)
    run(['cmake', '--install', build])
    library = prefix / 'lib/libxtb.so'
    executable = prefix / 'bin/xtb'
    parameters = prefix / 'share/xtb'
    for path in (library, executable, *(parameters / f'param_gfn{i}-xtb.txt' for i in range(3))):
        if not path.is_file():
            raise ValueError(f'Missing installed artifact: {path}')
    checks = verify_library(build, library.parent, parameters, cache['CMAKE_C_COMPILER'],
                            dependency_prefix, env)
    basis_checks = verify_basis(build, prefix, cache['CMAKE_C_COMPILER'], dependency_prefix, env)
    version = run([executable, '--version'], env=env, capture_output=True, text=True).stdout
    record_dir = parameters / 'rdkit-eth'
    record_dir.mkdir(exist_ok=True)
    for filename in ('COPYING', 'COPYING.LESSER'):
        shutil.copy2(source / filename, record_dir / filename)
    if patch:
        shutil.copy2(patch, record_dir / patch.name)
    basis_patch = TOOLS / lock["basis_extension"]["patch"]
    shutil.copy2(basis_patch, record_dir / basis_patch.name)
    source_diff = git_output(source, 'diff', 'HEAD', '--binary') + '\n'
    (record_dir / 'source.diff').write_text(source_diff)
    (record_dir / 'ctest.log').write_text(ctest.stdout + ctest.stderr)
    manifest = dict(schema_version=1, capabilities={'native_gfn0': True, 'orbital_basis_api': 1},
                    source=lock, source_directory=str(source), prefix=str(prefix),
                    built_at=datetime.now(timezone.utc).isoformat(),
                    source_diff_sha256=hashlib.sha256(source_diff.encode()).hexdigest(),
                    configure_command=list(map(str, configure)),
                    dependency_prefix=str(dependency_prefix) if dependency_prefix else None,
                    compilers={k: cache[k] for k in ('CMAKE_C_COMPILER', 'CMAKE_Fortran_COMPILER')},
                    cli_version_output=version.strip(),
                    files={str(p.relative_to(prefix)): digest(p) for p in
                           (library, executable, *(parameters / f'param_gfn{i}-xtb.txt' for i in range(3)))},
                    state_checks=checks, basis_checks=basis_checks)
    (record_dir / 'build.json').write_text(json.dumps(manifest, indent=2) + '\n')
    activation = '# Source this file to select this local xTB installation.\n'
    for key, value in [('XTB_LIBRARY', library), ('XTBPATH', parameters)]:
        activation += f'export {key}={shlex.quote(str(value))}\n'
    (prefix / 'activate.sh').write_text(activation)
    print(f'\nVerified patched xTB: {prefix}\n'
          f'12 method-switch sequences passed; provenance: {record_dir / "build.json"}\n'
          f'Activate with: source {shlex.quote(str(prefix / "activate.sh"))}')


if __name__ == '__main__':
    try:
        main()
    except subprocess.CalledProcessError as exc:
        if exc.stdout:
            print(exc.stdout, file=sys.stderr)
        if exc.stderr:
            print(exc.stderr, file=sys.stderr)
        sys.exit(f'Build/verification failed (exit {exc.returncode}): {shlex.join(map(str, exc.cmd))}')
    except (ValueError, OSError, subprocess.TimeoutExpired) as exc:
        sys.exit(f'Build/verification failed: {exc}')
