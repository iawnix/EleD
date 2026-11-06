"""Regression checks against the xTB CLI optimizer in the same runtime."""
import importlib.util
from pathlib import Path

import pytest

from rdkit_eth.errors import XTBNativeNotFoundError
from rdkit_eth.runtime import _extended_libxtb_build, _find_xtb_library


_TOOL = Path(__file__).resolve().parents[1] / "tools/benchmark_xtb_optimization.py"
_SPEC = importlib.util.spec_from_file_location("xtb_optimization_benchmark", _TOOL)
_MODULE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(_MODULE)


@pytest.fixture(scope="module")
def extended_runtime():
    try:
        library = Path(_find_xtb_library()).resolve()
    except XTBNativeNotFoundError:
        pytest.skip("extended libxtb unavailable")
    if not _extended_libxtb_build(str(library)):
        pytest.skip("reviewed extended libxtb required")
    executable = library.parent.parent / "bin" / "xtb"
    if not executable.is_file():
        pytest.skip("matching xTB CLI is unavailable")
    return library, executable


@pytest.mark.integration
@pytest.mark.parametrize("smiles, seed", [("O", 17), ("CCO", 18)])
def test_python_optimizer_matches_cli_regression(smiles, seed, extended_runtime, tmp_path):
    library, executable = extended_runtime
    result = _MODULE.run_case(
        smiles,
        executable=executable,
        library=str(library),
        method="GFN2-xTB",
        seed=seed,
        cycles=100,
        scf_max_iterations=250,
        accuracy=1.0,
        workdir=tmp_path,
    )
    assert result["python"]["success"]
    assert abs(result["comparison"]["energy_difference_hartree"]) < 1.0e-5
    assert result["comparison"]["aligned_rmsd_angstrom"] < 5.0e-2
