"""Geometric invariants for the research-only constraint replay."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture
def replay():
    path = Path(__file__).resolve().parents[1]/'tools/embedding_research/audit_constraints.py'
    spec = importlib.util.spec_from_file_location('constraint_audit', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.audit


def fixture_geometry():
    # Positive affine volume = 4, with the center strictly inside all faces.
    p = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                  [0., 0., 1.], [-1., -1., -1.]])
    c = {'bounds': [], 'chiral_sets': [{'atoms': [0, 1, 2, 3, 4],
         'lower': 5., 'upper': 100., 'flags': 0}], 'tetrahedral_sets': []}
    return p, c


def test_orientation_translation_rotation_and_mirror(replay):
    p, c = fixture_geometry()
    original = replay(p, c)
    assert original['chiral'][0]['volume'] == pytest.approx(4.)
    assert original['failed_predicates'] == []
    rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
    moved = replay(p @ rotation + [20., -13., 7.], c)
    assert moved['failed_predicates'] == []
    assert moved['chiral'][0]['volume'] == pytest.approx(4.)
    mirrored = replay(p * [-1., 1., 1.], c)
    assert mirrored['failed_predicates'] == ['CHIRAL_VOLUME']
    assert mirrored['chiral'][0]['center_inside']


def test_containment_is_not_the_same_as_ligand_handedness(replay):
    p, c = fixture_geometry()
    p[0] = [3., 0., 0.]
    result = replay(p, c)
    assert result['chiral'][0]['volume_ok']
    assert result['chiral'][0]['sign_ok']
    assert result['failed_predicates'] == ['FINAL_CENTER_IN_VOLUME']


def test_even_neighbor_permutation_preserves_affine_orientation(replay):
    p, c = fixture_geometry()
    p[0] = [3., 0., 0.]
    original = replay(p, c)
    c['chiral_sets'][0]['atoms'] = [0, 2, 4, 3, 1]
    permuted = replay(p, c)
    assert permuted['chiral'][0]['volume'] == pytest.approx(original['chiral'][0]['volume'])
    assert permuted['failed_predicates'] == original['failed_predicates']


def test_fourth_embedding_coordinate_does_not_change_3d_predicates(replay):
    p, c = fixture_geometry()
    higher_dimension = np.column_stack([p, [10., -30., 1., 0., 200.]])
    assert replay(higher_dimension, c) == replay(p, c)


def test_distance_failure_requires_union_membership_and_tolerance(replay):
    p, c = fixture_geometry()
    p = np.vstack([p, [100., 0., 0.]])
    # A raw bound violation below the final 10%-of-upper tolerance is allowed.
    c['bounds'] = [[0, 1, .5, .95], [0, 5, .5, 1.]]
    assert replay(p, c)['failed_predicates'] == []
    c['bounds'][0][-1] = .8
    assert replay(p, c)['failed_predicates'] == ['FINAL_CHIRAL_BOUNDS']


def test_nonfinite_coordinates_are_not_interpreted_as_success(replay):
    p, c = fixture_geometry()
    p[0, 0] = np.nan
    with pytest.raises(ValueError, match='nonfinite'):
        replay(p, c)
