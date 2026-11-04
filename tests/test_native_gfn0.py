"""Native GFN0 needs evidence for the selected binary, not just a version."""

from copy import deepcopy
import hashlib
import json
from types import SimpleNamespace

import pytest

from rdkit_eth.options import XTBOptions
from rdkit_eth.runtime import _D4_LOCAL_REFERENCES_PATCH, _native_gfn0_build
from rdkit_eth.xtb import _resolve_options


@pytest.fixture
def verified_install(tmp_path):
    library = tmp_path / 'lib/libxtb.so.6.7.1'
    library.parent.mkdir()
    library.write_bytes(b'patched test binary')
    (library.parent / 'libxtb.so').symlink_to(library.name)
    manifest = tmp_path / 'share/xtb/rdkit-eth/build.json'
    manifest.parent.mkdir(parents=True)
    record = {'schema_version': 1, 'capabilities': {'native_gfn0': True},
              'source': {'patch_sha256': _D4_LOCAL_REFERENCES_PATCH, 'commit': 'test-commit'},
              'files': {'lib/libxtb.so': hashlib.sha256(library.read_bytes()).hexdigest()}}
    manifest.write_text(json.dumps(record))
    return library, manifest, record


def test_verified_binary_is_recognized_through_its_symlink(verified_install):
    library, manifest, _ = verified_install
    evidence = _native_gfn0_build(str(library))
    assert evidence['manifest'] == str(manifest)
    assert evidence['source_commit'] == 'test-commit'
    assert _native_gfn0_build(str(library.parent / 'libxtb.so')) == evidence


@pytest.mark.parametrize('damage', ['missing', 'invalid-json', 'non-object', 'schema',
                                  'capability', 'patch', 'hash', 'binary', 'malformed-fields'])
def test_stale_or_unverified_builds_are_not_enabled(verified_install, damage):
    library, manifest, record = verified_install
    record = deepcopy(record)
    if damage == 'missing':
        manifest.unlink()
    elif damage == 'invalid-json':
        manifest.write_text('{')
    elif damage == 'non-object':
        manifest.write_text('[]')
    elif damage == 'binary':
        assert _native_gfn0_build(str(library))
        library.write_bytes(b'unpatched replacement')
    else:
        if damage == 'schema':
            record['schema_version'] = 99
        elif damage == 'capability':
            record['capabilities']['native_gfn0'] = False
        elif damage == 'patch':
            record['source']['patch_sha256'] = 'unknown-patch'
        elif damage == 'hash':
            record['files']['lib/libxtb.so'] = 'stale-hash'
        elif damage == 'malformed-fields':
            record['capabilities'] = []
        manifest.write_text(json.dumps(record))
    assert _native_gfn0_build(str(library)) is None


def test_gfn0_routing_changes_with_library_verification(verified_install, monkeypatch):
    library, manifest, _ = verified_install
    monkeypatch.setattr('rdkit_eth.xtb._find_xtb_library', lambda _: str(library))
    monkeypatch.setattr('rdkit_eth.xtb._NativeXTB', lambda path: SimpleNamespace(library=path))
    options = XTBOptions(method='gfn0', backend='native')
    assert _resolve_options(options).backend == 'native'
    assert _resolve_options(XTBOptions(method='gfn0', backend='native')).backend == 'native'
    manifest.unlink()
    with pytest.raises(ValueError, match='verified patched'):
        _resolve_options(XTBOptions(method='gfn0', backend='native'))
    # GFN1/2 do not require a GFN0 capability record.
    assert _resolve_options(XTBOptions(method='gfn2', backend='native')).backend == 'native'


@pytest.mark.parametrize('damage', [None, 'commit', 'fix-commit', 'kind', 'file-hash', 'malformed'])
def test_upstream_fix_requires_reviewed_identity(verified_install, damage):
    from rdkit_eth.runtime import _UPSTREAM_GFN0_COMMIT, _UPSTREAM_D4_SHA256
    library, manifest, record = verified_install
    record['source'] = {'commit': _UPSTREAM_GFN0_COMMIT,
                        'gfn0_fix': {'kind': 'upstream', 'commit': _UPSTREAM_GFN0_COMMIT,
                                     'file_sha256': _UPSTREAM_D4_SHA256}}
    if damage == 'commit': record['source']['commit'] = 'another-commit'
    elif damage == 'fix-commit': record['source']['gfn0_fix']['commit'] = 'another-commit'
    elif damage == 'kind': record['source']['gfn0_fix']['kind'] = 'version-only'
    elif damage == 'file-hash': record['source']['gfn0_fix']['file_sha256'] = 'unknown'
    elif damage == 'malformed': record['source']['gfn0_fix'] = []
    manifest.write_text(json.dumps(record))
    result = _native_gfn0_build(str(library))
    if damage is None:
        assert result['fix_kind'] == 'upstream'
        assert result['patch_sha256'] is None
    else:
        assert result is None
