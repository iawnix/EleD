"""Unit checks for the standalone xTB optimization comparison tool."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


_PATH = Path(__file__).resolve().parents[1] / "tools/benchmark_xtb_optimization.py"
_SPEC = importlib.util.spec_from_file_location("xtb_optimization_benchmark", _PATH)
benchmark = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(benchmark)


def test_aligned_rmsd_removes_rigid_translation_and_rotation():
    reference = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.2, 0.7, 0.0]])
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    candidate = reference @ rotation + np.array([3.0, -4.0, 2.0])
    assert benchmark._aligned_rmsd(reference, candidate) == pytest.approx(0.0, abs=1e-12)


def test_read_xyz_extracts_comment_and_coordinates(tmp_path):
    path = tmp_path / "optimized.xyz"
    path.write_text("2\n energy: -1.25 gnorm: 0.0001\nH 0 0 0\nH 0 0 0.74\n")
    symbols, coordinates, comment = benchmark._read_xyz(path)
    assert symbols == ["H", "H"]
    np.testing.assert_allclose(coordinates, [[0, 0, 0], [0, 0, 0.74]])
    assert "gnorm" in comment
