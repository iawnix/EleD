"""A gain cannot conceal a lost baseline success or an untested critical stratum."""
import importlib.util
from pathlib import Path

import pytest


@pytest.fixture
def gate():
    path = Path(__file__).resolve().parents[1]/'tools/embedding_research/nonregression_gate.py'
    spec = importlib.util.spec_from_file_location('nonregression_gate', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def report(initial=True, post=True, energy=12.):
    return {'initial_valid': initial, 'post_valid': post,
            'post_energy_kcal_mol': energy if post else None}


def test_many_recoveries_do_not_cancel_a_single_lost_success(gate):
    gains = [gate.compare_pair(report(False, False), report()) for _ in range(100)]
    loss = gate.compare_pair(report(), report(False, False))
    result = gate.verdict(gains+[loss], affected_baseline_successes=1)
    assert result['gained_pairs'] == 100
    assert result['regressed_pairs'] == 1
    assert result['status'] == 'reject_observed_regression'


def test_no_baseline_successes_is_inconclusive_even_with_perfect_recovery(gate):
    result = gate.verdict([gate.compare_pair(report(False, False), report())], 0)
    assert result['status'] == 'inconclusive_no_affected_baseline_success'
    assert not result['production_ready']


def test_post_optimization_loss_is_independent_of_initial_success(gate):
    comparison = gate.compare_pair(report(), report(True, False))
    assert comparison['failures'] == ['lost_post_optimization_success']


def test_energy_worsening_is_not_hidden_by_preserved_stereo(gate):
    comparison = gate.compare_pair(report(energy=12.), report(energy=13.))
    assert comparison['failures'] == ['higher_post_mmff_energy']
    assert comparison['post_energy_delta_kcal_mol'] == 1.


def test_numerical_tolerance_and_finite_coverage_are_explicit(gate):
    comparison = gate.compare_pair(report(energy=12.), report(energy=12.0001))
    assert comparison['failures'] == []
    result = gate.verdict([comparison], 1)
    assert result['status'] == 'no_regression_observed_in_tested_pairs'
    assert result['production_ready'] is False
