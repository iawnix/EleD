import pytest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.eas.modeling import (aggregate_molecule_metrics, molecule_metrics,
                                     paired_bootstrap)
from benchmarks.eas.train import fit_and_evaluate


pytestmark = pytest.mark.benchmark


def _rows_and_scores():
    rows = [
        {'record_id': 'a', 'target': 1}, {'record_id': 'a', 'target': 0},
        {'record_id': 'b', 'target': 1}, {'record_id': 'b', 'target': 0},
        {'record_id': 'c', 'target': 1},
    ]
    return rows, [0.9, 0.1, 0.2, 0.8, 0.5]


def test_molecule_metrics_are_molecule_weighted_and_skip_undefined_auc():
    rows, scores = _rows_and_scores()
    metrics = molecule_metrics(rows, scores)
    assert metrics['a']['macro_roc_auc'] == 1.0
    assert metrics['b']['macro_roc_auc'] == 0.0
    assert metrics['c']['macro_roc_auc'] is None
    aggregate = aggregate_molecule_metrics(metrics)
    assert aggregate['macro_roc_auc'] == pytest.approx(.5)
    assert aggregate['macro_roc_auc_n'] == 2


def test_paired_bootstrap_is_reproducible_and_checks_ids():
    rows, scores = _rows_and_scores()
    baseline = molecule_metrics(rows, scores)
    candidate = molecule_metrics(rows, [0.8, 0.2, 0.7, 0.3, 0.4])
    first = paired_bootstrap(baseline, candidate, replicates=100, seed=7)
    second = paired_bootstrap(baseline, candidate, replicates=100, seed=7)
    assert first == second
    with pytest.raises(ValueError):
        paired_bootstrap(baseline, {'other': candidate['a']}, replicates=5)


def test_logistic_training_keeps_test_molecules_untouched():
    rows = []
    for molecule in range(10):
        split = 'train' if molecule < 8 else 'test'
        fold = molecule % 4
        rows.extend([
            {'record_id': f'm{molecule}', 'atom_index': 0, 'target': 1,
             'split': split, 'cv_fold': fold,
             'features': {'signal': float(molecule) + 1.0, 'missing': None}},
            {'record_id': f'm{molecule}', 'atom_index': 1, 'target': 0,
             'split': split, 'cv_fold': fold,
             'features': {'signal': -float(molecule) - 1.0, 'missing': None}},
        ])
    result = fit_and_evaluate(rows, model_name='logistic_regression', seed=7)
    assert result['n_test_molecules'] == 2
    assert set(result['per_molecule']) == {'m8', 'm9'}
    assert result['metrics']['macro_roc_auc'] is not None
