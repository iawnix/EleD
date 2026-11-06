"""The unified installer must not silently reuse drifted environments."""
import importlib.util
import json
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location('native_installer', Path(__file__).resolve().parents[1]/'tools/install_native.py')
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


def test_locked_environment_reuse_rejects_changed_packages(tmp_path):
    url = 'https://example.invalid/linux-64/package-1.0-0.conda'
    lock = tmp_path / 'lock.txt'; lock.write_text('@EXPLICIT\n' + url + '#known-sha256\n')
    meta = tmp_path / 'prefix/conda-meta'; meta.mkdir(parents=True)
    record = meta / 'package.json'
    record.write_text(json.dumps({'url':url, 'sha256':'known-sha256'}))
    installer.verify_environment(meta.parent,lock)
    record.write_text(json.dumps({'url':url, 'sha256':'changed'}))
    with pytest.raises(ValueError, match='fresh --directory'):
        installer.verify_environment(meta.parent,lock)
    assert json.loads(record.read_text())['sha256'] == 'changed'


def test_nonempty_unmanaged_directory_is_preserved(tmp_path):
    prefix = tmp_path/'prefix'; prefix.mkdir(); work = prefix/'research'; work.write_text('keep')
    with pytest.raises(ValueError, match='Refusing nonempty'):
        installer.environment('not-run', prefix, tmp_path/'unused-lock', {})
    assert work.read_text() == 'keep'


def test_install_checks_ignore_external_engine_selection(monkeypatch):
    for name in ('XTB_LIBRARY','XTBPATH','LD_PRELOAD','PYTHONPATH','RDXTB_XTB_PREFIX'):
        monkeypatch.setenv(name,'another-installation')
    env = installer.clean_environment()
    assert all(name not in env for name in ('XTB_LIBRARY','XTBPATH','LD_PRELOAD','PYTHONPATH','RDXTB_XTB_PREFIX'))
    assert env['OMP_NUM_THREADS'] == '1'


def test_release_builder_runs_outside_checkout_build_directory():
    """The repository's build/ directory must not shadow the build frontend."""
    # This check is intentionally source-level: the actual wheel build is
    # covered by the locked installed-profile workflow.
    release = Path(__file__).resolve().parents[1] / 'tools/build_release.py'
    text = release.read_text()
    assert 'cwd=output' in text
    assert 'cwd=temporary' in text
