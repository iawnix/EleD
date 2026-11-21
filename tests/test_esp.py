"""ESP extension provenance and installed C ABI/CLI numerical contracts."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('verify_esp', ROOT/'tools/diagnostics/verify_esp.py')
verifier = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verifier)


def test_esp_patch_matches_source_lock():
    extension = json.loads((ROOT/'tools/xtb-source.json').read_text())['esp_extension']
    assert extension['api_version'] == 1
    assert hashlib.sha256((ROOT/'tools'/extension['patch']).read_bytes()).hexdigest() == extension['patch_sha256']


@pytest.mark.integration
def test_esp_c_abi_matches_cli_and_preserves_result_lifecycle():
    from eled.runtime import _find_xtb_library
    from eled import XTBNativeNotFoundError
    try:
        library = Path(_find_xtb_library()).resolve()
    except XTBNativeNotFoundError:
        pytest.skip('extended libxtb not installed')
    prefix = library.parent.parent
    manifest = prefix/'share/xtb/eled/build.json'
    if not manifest.exists() or json.loads(manifest.read_text()).get('capabilities', {}).get('esp_api') != 1:
        pytest.skip('ESP extension not installed')
    rows = verifier.verify(prefix, library)
    assert len(rows) == 15
    assert max(row.get('max_cli_error_au', 0) for row in rows) < 2e-7
