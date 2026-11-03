"""Source preparation must preserve edits and never silently change revisions."""

import hashlib
import importlib.util
from pathlib import Path
import subprocess

import pytest

spec = importlib.util.spec_from_file_location(
    'xtb_builder', Path(__file__).resolve().parents[1] / 'tools/build_patched_xtb.py')
builder = importlib.util.module_from_spec(spec)
spec.loader.exec_module(builder)


@pytest.fixture
def source_case(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    subprocess.run(['git', 'init', '-q', source], check=True)
    (source / 'model.txt').write_text('global references\n')
    subprocess.run(['git', '-C', source, 'add', 'model.txt'], check=True)
    subprocess.run(['git', '-C', source, '-c', 'user.name=Build test',
                    '-c', 'user.email=build-test@example.invalid',
                    '-c', 'commit.gpgsign=false', 'commit', '-qm', 'baseline'], check=True)
    commit = subprocess.check_output(['git', '-C', source, 'rev-parse', 'HEAD'], text=True).strip()
    patch = tmp_path / 'model.patch'
    patch.write_text('--- a/model.txt\n+++ b/model.txt\n@@ -1 +1 @@\n'
                     '-global references\n+local references\n')
    lock = {'commit': commit, 'patched_file': 'model.txt',
            'patched_file_sha256': hashlib.sha256(b'local references\n').hexdigest()}
    return source, lock, patch


def test_patch_can_be_applied_and_reused_without_resetting(source_case, tmp_path):
    source, lock, patch = source_case
    for _ in range(2):
        assert builder.source_tree(tmp_path, source, lock, patch) == source
        assert (source / 'model.txt').read_text() == 'local references\n'
        assert builder.git_output(source, 'rev-parse', 'HEAD') == lock['commit']


def test_wrong_revision_is_rejected_without_checkout(source_case, tmp_path):
    source, lock, patch = source_case
    original_head = lock['commit']
    lock['commit'] = '0' * 40
    with pytest.raises(ValueError, match='Source HEAD'):
        builder.source_tree(tmp_path, source, lock, patch)
    assert builder.git_output(source, 'rev-parse', 'HEAD') == original_head
    assert (source / 'model.txt').read_text() == 'global references\n'


def test_unrelated_work_is_preserved(source_case, tmp_path):
    source, lock, patch = source_case
    (source / 'new-research.txt').write_text('keep this work\n')
    with pytest.raises(ValueError, match='unrelated changes'):
        builder.source_tree(tmp_path, source, lock, patch)
    assert (source / 'new-research.txt').read_text() == 'keep this work\n'
    assert (source / 'model.txt').read_text() == 'global references\n'


def test_modified_target_is_not_overwritten(source_case, tmp_path):
    source, lock, patch = source_case
    (source / 'model.txt').write_text('experimental references\n')
    with pytest.raises(ValueError, match='reviewed D4 patch'):
        builder.source_tree(tmp_path, source, lock, patch)
    assert (source / 'model.txt').read_text() == 'experimental references\n'


def test_basis_patch_is_idempotent_and_modified_source_is_preserved(source_case, tmp_path):
    source, lock, patch = source_case
    extra = tmp_path / 'basis.patch'
    # A tracked baseline models the real extension's edits to existing files.
    (source / 'basis.txt').write_text('old\n')
    subprocess.run(['git', '-C', source, 'add', 'basis.txt'], check=True)
    subprocess.run(['git', '-C', source, '-c', 'user.name=Build test',
                    '-c', 'user.email=build-test@example.invalid', '-c', 'commit.gpgsign=false',
                    'commit', '-qm', 'basis baseline'], check=True)
    lock['commit'] = builder.git_output(source, 'rev-parse', 'HEAD')
    extra.write_text('--- a/basis.txt\n+++ b/basis.txt\n@@ -1 +1 @@\n-old\n+snapshot\n')
    lock['basis_extension'] = {'patch': str(extra), 'patch_sha256': builder.digest(extra),
                              'patched_files': {'basis.txt': hashlib.sha256(b'snapshot\n').hexdigest()}}
    for _ in range(2):
        builder.source_tree(tmp_path, source, lock, patch)
        assert (source / 'basis.txt').read_text() == 'snapshot\n'
    (source / 'basis.txt').write_text('experimental\n')
    with pytest.raises(ValueError, match='reviewed basis extension'):
        builder.source_tree(tmp_path, source, lock, patch)
    assert (source / 'basis.txt').read_text() == 'experimental\n'
    extra.write_text('tampered patch\n')
    with pytest.raises(ValueError, match='patch hash differs'):
        builder.source_tree(tmp_path, source, lock, patch)


def test_native_capability_constants_match_reviewed_patch():
    import json
    from rdkit_eth.runtime import _BASIS_EXPORT_PATCH
    root = Path(__file__).resolve().parents[1]
    lock = json.loads((root / 'tools/xtb-source.json').read_text())['basis_extension']
    assert lock['patch_sha256'] == _BASIS_EXPORT_PATCH == builder.digest(root / 'tools' / lock['patch'])
    assert _BASIS_EXPORT_PATCH in (root / 'cpp/rdXTBTools/XTBTools.cpp').read_text()
