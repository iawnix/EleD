import json

from rdkit_eth.diagnostics import main


def test_inventory_does_not_require_xtb(monkeypatch, tmp_path, capsys):
    monkeypatch.setenv("XTB_LIBRARY", str(tmp_path / "missing.so"))
    assert main(["doctor", "--json"]) == 0
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "inventory"
    assert report["compiled_module_import_checked"] is False
    assert report["package_path"].endswith("rdkit_eth/__init__.py")


def test_explicit_broken_library_reports_failure_without_fallback(tmp_path, capsys):
    missing = str(tmp_path / "absent-libxtb.so")
    assert main(["doctor", "--check", "native", "--library", missing, "--json"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert report["status"] == "error"
    assert missing in report["error"]
    assert "check" not in report


def test_removed_calculation_backend_and_executable_flag_are_rejected():
    import pytest
    for args in (["doctor", "--check", "cli"], ["doctor", "--executable", "/bin/xtb"]):
        with pytest.raises(SystemExit) as error:
            main(args)
        assert error.value.code == 2
