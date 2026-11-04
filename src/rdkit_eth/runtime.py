"""Locate xTB and identify the actual installation used for provenance/cache keys."""

from __future__ import annotations
import ctypes
import ctypes.util
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path
import sys

from .errors import XTBNativeNotFoundError, XTBRunError


# The reviewed D4 repair, independent of a particular compiled binary. A
# matching build record is emitted only after installed-library regression.
_UPSTREAM_GFN0_COMMIT = '33fe381adde775c4992610c0f90d8914b7572bb5'
_UPSTREAM_D4_SHA256 = '5c4022078296e46dfe6a21a36505ec07402bb1eee88448d735c588676dd6e61c'

# Accept previously verified 0.4 installations for reproducibility.
_D4_LOCAL_REFERENCES_PATCH = '2f503f401c410b1456824efe906b89d5f932e8a8f9c63f817669a259db04c859'


def _native_gfn0_build(path: str) -> dict | None:
    """Recognize a verified local build without initializing GFN0 to probe it.

    An API/program version alone does not identify the D4 fix. Check the build
    capability and the actual selected library's digest. Missing, stale, or
    malformed records simply leave GFN0 unavailable through the native API.
    """
    try:
        library = Path(path).resolve(strict=True)
        prefix = library.parent.parent
        record = prefix / 'share/xtb/rdkit-eth/build.json'
        data = json.loads(record.read_text())
        if not isinstance(data, dict) or data.get('schema_version') != 1:
            return None
        capabilities = data.get('capabilities')
        source = data.get('source')
        files = data.get('files')
        if not all(isinstance(value, dict) for value in (capabilities, source, files)):
            return None
        fix = source.get('gfn0_fix', {})
        upstream = (isinstance(fix, dict) and fix.get('kind') == 'upstream'
                    and fix.get('commit') == source.get('commit') == _UPSTREAM_GFN0_COMMIT
                    and fix.get('file_sha256') == _UPSTREAM_D4_SHA256)
        legacy = source.get('patch_sha256') == _D4_LOCAL_REFERENCES_PATCH
        if capabilities.get('native_gfn0') is not True or not (upstream or legacy):
            return None
        if (prefix / 'lib/libxtb.so').resolve(strict=True) != library:
            return None
        library_hash = _fingerprint(library)
        if files.get('lib/libxtb.so') != library_hash:
            return None
        return {'manifest': str(record), 'manifest_sha256': _fingerprint(record),
                'library_sha256': library_hash, 'patch_sha256': source.get('patch_sha256'),
                'fix_kind': 'upstream' if upstream else 'local-patch',
                'source_commit': source.get('commit')}
    except (OSError, ValueError):
        return None


_BASIS_EXPORT_PATCH = '30f93f19b52a688727cb8af30cd0b856fe47977a96bef3c52436be707bcb75a2'


def _native_basis_build(path: str) -> dict | None:
    """Recognize the reviewed optional basis ABI for the selected binary."""
    try:
        library = Path(path).resolve(strict=True)
        prefix = library.parent.parent
        record = prefix / 'share/xtb/rdkit-eth/build.json'
        data = json.loads(record.read_text())
        extension = data['source']['basis_extension']
        if (data['schema_version'] != 1 or data['capabilities']['orbital_basis_api'] != 1
                or extension['api_version'] != 1 or extension['patch_sha256'] != _BASIS_EXPORT_PATCH
                or (prefix / 'lib/libxtb.so').resolve(strict=True) != library):
            return None
        fingerprint = _fingerprint(library)
        if data['files']['lib/libxtb.so'] != fingerprint:
            return None
        return {'api_version': 1, 'library_sha256': fingerprint,
                'patch_sha256': extension['patch_sha256'], 'manifest': str(record),
                'manifest_sha256': _fingerprint(record), 'source_commit': data['source']['commit']}
    except (OSError, ValueError, KeyError, TypeError):
        return None


def _find_xtb_library(library=None) -> str:
    candidate = os.fspath(library) if library is not None else os.environ.get('XTB_LIBRARY')
    patterns = ('libxtb.so*', 'libxtb*.dylib', 'libxtb*.dll', 'xtb.dll')
    candidates = []
    if candidate:
        root = Path(candidate)
        if root.is_dir():
            candidates = [str(p) for pattern in patterns for p in sorted(root.glob(pattern))]
        else:
            candidates = [candidate]
    else:
        for prefix in (sys.prefix, os.environ.get('CONDA_PREFIX'), os.environ.get('XTB_PREFIX')):
            if prefix:
                for subdir in ('lib', 'Library/bin', 'bin'):
                    for pattern in patterns:
                        candidates.extend(str(p) for p in sorted((Path(prefix) / subdir).glob(pattern)))
        found = ctypes.util.find_library('xtb')
        if found:
            candidates.append(found)
    errors = []
    for path in dict.fromkeys(candidates):
        if path.endswith('.a'):
            continue
        try:
            ctypes.CDLL(path)
        except OSError as exc:
            errors.append(str(exc))
        else:
            return str(Path(path).resolve()) if Path(path).exists() else path
    detail = '; '.join(errors[-2:])
    raise XTBNativeNotFoundError(
        f'could not load shared xTB library {candidate or "libxtb"!r}'
        + (f': {detail}' if detail else '')
    )


@lru_cache(maxsize=64)
def _file_digest(path: str, size: int, mtime_ns: int, ctime_ns: int) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def _fingerprint(path: Path) -> str:
    stat = path.stat()
    return _file_digest(str(path.resolve()), stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def _engine_identity(backend: str, path: str) -> dict:
    location = Path(path)
    digest = _fingerprint(location) if location.is_file() else None
    version = None
    # The C API exposes an ABI version, not the xTB program version. Record a
    # package version only when this file belongs to an identifiable conda install.
    if backend == 'native' and location.is_file():
        for parent in location.parents:
            records = parent / 'conda-meta'
            if records.is_dir():
                for record in records.glob('xtb-*.json'):
                    data = json.loads(record.read_text())
                    if data.get('name') == 'xtb' and str(location.relative_to(parent)) in data.get('files', []):
                        version = data['version']
                break
    identity = {'backend': backend, 'path': path, 'sha256': digest, 'version': version}
    if backend == 'native':
        identity['native_gfn0'] = _native_gfn0_build(path)
        identity['orbital_basis'] = _native_basis_build(path)
    return identity


def _parameter_context() -> dict:
    names = ('XTBPATH', 'XTBHOME', 'XTB_PREFIX', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')
    environment = {name: os.environ.get(name) for name in names}
    roots = [Path.cwd()]
    roots.extend(Path(p) for p in os.environ.get('XTBPATH', '').split(os.pathsep) if p)
    files = {str(p.resolve()): _fingerprint(p) for root in roots for p in root.glob('param_gfn*.txt') if p.is_file()}
    return {'environment': environment, 'parameter_files': files}


def _parameter_file(method: str, engine_path: str) -> dict | None:
    if method == 'GFN-FF':
        return None
    filename = f'param_{method.lower()}.txt'
    roots = [Path(p) for p in os.environ.get('XTBPATH', '').split(os.pathsep) if p]
    location = Path(engine_path)
    roots.extend(parent / 'share' / 'xtb' for parent in location.parents)
    for root in roots:
        path = root / filename
        if path.is_file():
            return {'path': str(path.resolve()), 'sha256': _fingerprint(path)}
    if method == 'GFN0-xTB':
        raise XTBRunError(f'cannot locate {filename}; set XTBPATH to the xTB parameter directory')
    return None  # GFN1/2 have built-in parameters as well.
