"""Mathematical invariants for AO conversion and nonorthogonal MO comparison."""

import numpy as np
import pytest

from eled.orbital_alignment import (
    compare_orbitals, libxtb_to_tblite_rows, validate_tblite_basis,
)


def test_nonorthogonal_degenerate_rotation_and_phase_are_equivalent():
    s = np.array([[1., .2, .1], [.2, 1., .3], [.1, .3, 1.]])
    c = np.linalg.inv(np.linalg.cholesky(s).T)
    theta = .73
    rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                         [np.sin(theta), np.cos(theta), 0], [0, 0, -1]])
    candidate = c @ rotation
    result = compare_orbitals(c, candidate, s, [-2, -2, 1], [-2, -2, 1], [2, 2, 0], [2, 2, 0])
    assert result["sign_aligned_max_abs"] > .1
    assert result["subspace_aligned_max_abs"] < 1e-14
    assert result["max_projector_distance"] < 1e-14
    assert result["density_orthonormal_max_abs"] < 1e-14
    assert [g["size"] for g in result["groups"]] == [2, 1]


def test_rotating_distinct_or_unequally_occupied_orbitals_is_detected():
    c = np.eye(2)
    rotation = np.array([[.8, -.6], [.6, .8]])
    for energies, occupations in [([-1, 1], [2, 2]), ([-1, -1], [2, 0])]:
        result = compare_orbitals(c, rotation, c, energies, energies, occupations, occupations)
        assert [g["size"] for g in result["groups"]] == [1, 1]
        assert result["max_projector_distance"] == pytest.approx(.6)
        assert result["subspace_aligned_max_abs"] > .5
    assert result["density_orthonormal_max_abs"] > .5


def test_clustering_checks_both_full_energy_spans_and_reports_norm_errors():
    c = np.eye(3)
    r = compare_orbitals(c, c * 2, c, [0, .0007, .0014], [0, .0008, .002], [2]*3, [2]*3)
    assert [g["size"] for g in r["groups"]] == [2, 1]  # no transitive gap chaining
    assert r["candidate_orthonormality_error"] == pytest.approx(3)
    r = compare_orbitals(c, c, c, [0, 0, 0], [0, .002, .004], [2]*3, [2]*3)
    assert [g["size"] for g in r["groups"]] == [1, 1, 1]


@pytest.mark.parametrize("overlap", [np.array([[1., .5], [0., 1.]]), np.diag([1., -1.])])
def test_invalid_metric_is_rejected(overlap):
    with pytest.raises((ValueError, np.linalg.LinAlgError)):
        compare_orbitals(np.eye(2), np.eye(2), overlap, [-1, 1], [-1, 1], [2, 0], [2, 0])


def test_shell_topology_and_integer_maps_are_validated(tmp_path):
    parameters = tmp_path / "param.txt"
    parameters.write_text("$Z= 16\n ao=3s3p3d\n$end\n$Z= 1\n ao=1s2s\n$end\n")
    orbital_map = [0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 4]
    args = ([16, 1], parameters, orbital_map, [0, 0, 0, 1, 1], [0, 1, 2, 0, 0], 11)
    rows, phases, labels = libxtb_to_tblite_rows(*args)
    assert sorted(rows) == list(range(11))
    assert [(l["atom"], l["n"], l["component"]) for l in labels[-2:]] == [(1, 1, "s"), (1, 2, "s")]
    assert list(phases).count(-1) == 1
    with pytest.raises(ValueError, match="parameter shells differ"):
        libxtb_to_tblite_rows([1, 16], *args[1:])
    for bad_map in [orbital_map[:-1], [0]*11, np.asarray(orbital_map, dtype=float)]:
        with pytest.raises(ValueError):
            validate_tblite_basis(bad_map, args[3], args[4], 2, 11)
    with pytest.raises(ValueError, match="input order"):
        validate_tblite_basis(orbital_map, [0, 0, 1, 0, 1], args[4], 2, 11)


def test_f_shells_fail_instead_of_guessing_a_conversion(tmp_path):
    parameters = tmp_path / "param.txt"
    parameters.write_text("$Z= 58\n ao=4f\n$end\n")
    with pytest.raises(ValueError, match="only verified for s/p/d"):
        libxtb_to_tblite_rows([58], parameters, [0]*7, [0], [3], 7)
