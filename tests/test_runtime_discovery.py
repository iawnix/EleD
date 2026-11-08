"""Native libraries inherited by an isolated venv, with explicit override priority."""
from pathlib import Path

import pytest

from eled import runtime
from eled.errors import XTBNativeNotFoundError


@pytest.fixture
def prefixes(tmp_path, monkeypatch):
    roots = [tmp_path/name for name in ('venv', 'base', 'activated')]
    for root in roots: (root/'lib').mkdir(parents=True)
    monkeypatch.setattr(runtime.sys, 'prefix', str(roots[0]))
    monkeypatch.setattr(runtime.sys, 'base_prefix', str(roots[1]))
    monkeypatch.setenv('CONDA_PREFIX', str(roots[2]))
    monkeypatch.delenv('XTB_PREFIX', raising=False)
    monkeypatch.delenv('XTB_LIBRARY', raising=False)
    monkeypatch.setattr(runtime.ctypes.util, 'find_library', lambda name: None)
    def load(path):
        if not Path(path).is_file(): raise OSError('missing')
        return object()
    monkeypatch.setattr(runtime.ctypes, 'CDLL', load)
    return roots


def test_venv_prefers_own_library_then_base_before_activated_environment(prefixes):
    venv, base, activated = prefixes
    for root in (base, activated): (root/'lib/libxtb.so').touch()
    assert runtime._find_xtb_library() == str(base/'lib/libxtb.so')
    (venv/'lib/libxtb.so').touch()
    assert runtime._find_xtb_library() == str(venv/'lib/libxtb.so')


def test_explicit_library_failure_never_falls_back_to_base(prefixes, monkeypatch):
    venv, base, _ = prefixes; (base/'lib/libxtb.so').touch()
    with pytest.raises(XTBNativeNotFoundError): runtime._find_xtb_library(venv/'missing.so')
    monkeypatch.setenv('XTB_LIBRARY', str(venv/'missing.so'))
    with pytest.raises(XTBNativeNotFoundError): runtime._find_xtb_library()
    assert runtime._find_xtb_library(base/'lib/libxtb.so') == str(base/'lib/libxtb.so')
