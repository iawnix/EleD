"""Optional sklearn training/evaluation for the isolated EAS benchmark."""
from __future__ import annotations

from pathlib import Path
import json

import numpy as np

from .modeling import aggregate_molecule_metrics, molecule_metrics, paired_bootstrap
from .pipeline import feature_matrix, load_feature_archive


def _sklearn_components(model_name, seed):
    try:
        from sklearn.ensemble import GradientBoostingClassifier
        from sklearn.impute import SimpleImputer
        from sklearn.linear_model import LogisticRegression
        from sklearn.pipeline import Pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.model_selection import GridSearchCV, PredefinedSplit
    except ImportError as exc:
        raise RuntimeError('install tools/benchmark-requirements.txt for model training') from exc
    if model_name == 'logistic_regression':
        estimator = Pipeline([
            ('imputer', SimpleImputer(strategy='median', add_indicator=True)),
            ('scale', StandardScaler()),
            ('model', LogisticRegression(max_iter=2000, class_weight='balanced',
                                         random_state=seed)),
        ])
        grid = {'model__C': [0.1, 1.0, 10.0]}
    elif model_name == 'gradient_boosting':
        estimator = Pipeline([
            ('imputer', SimpleImputer(strategy='median', add_indicator=True)),
            ('model', GradientBoostingClassifier(random_state=seed)),
        ])
        grid = {'model__n_estimators': [100, 200], 'model__max_depth': [1, 2]}
    else:
        raise ValueError('unknown model_name')
    return estimator, grid, GridSearchCV, PredefinedSplit


def _matrix_arrays(rows, names=None):
    names = sorted(names or {key for row in rows for key in row['features']})
    x = np.asarray([[row['features'].get(name) for name in names] for row in rows], dtype=float)
    y = np.asarray([row['target'] for row in rows], dtype=int)
    return x, y, names


def _cv_splits(rows):
    folds = sorted({row['cv_fold'] for row in rows})
    splits = []
    for fold in folds:
        train = np.asarray([i for i, row in enumerate(rows) if row['cv_fold'] != fold], dtype=int)
        valid = np.asarray([i for i, row in enumerate(rows) if row['cv_fold'] == fold], dtype=int)
        if len(train) and len(valid) and len(set(rows[i]['target'] for i in train)) == 2 and len(set(rows[i]['target'] for i in valid)) == 2:
            splits.append((train, valid))
    return splits


def fit_and_evaluate(rows, *, model_name='logistic_regression', seed=20261106):
    """Fit on the primary training partition and evaluate the untouched test set."""
    train_rows = [row for row in rows if row['split'] == 'train']
    test_rows = [row for row in rows if row['split'] == 'test']
    if not train_rows or not test_rows:
        raise ValueError('both train and test rows are required')
    x_train, y_train, names = _matrix_arrays(train_rows)
    # A feature that is entirely missing in the training partition cannot be
    # imputed or learned. Drop it based on training data and record the choice;
    # partially observed columns retain both imputed values and missing masks.
    observed = np.isfinite(x_train).any(axis=0)
    dropped = [name for name, keep in zip(names, observed) if not keep]
    names = [name for name, keep in zip(names, observed) if keep]
    x_train = x_train[:, observed]
    if not names:
        raise ValueError('all feature columns are missing in the training partition')
    x_test, y_test, _ = _matrix_arrays(test_rows, names)
    if len(set(y_train)) != 2 or len(set(y_test)) != 2:
        raise ValueError('train and test partitions must contain both target classes')
    estimator, grid, GridSearchCV, PredefinedSplit = _sklearn_components(model_name, seed)
    cv = _cv_splits(train_rows)
    if cv:
        # PredefinedSplit uses -1 for samples excluded from a validation fold;
        # GridSearchCV receives explicit folds below to preserve molecule group
        # assignments recorded by the feature protocol.
        search = GridSearchCV(estimator, grid, scoring='roc_auc', cv=cv,
                              refit=True, error_score='raise')
        search.fit(x_train, y_train)
        fitted, best = search.best_estimator_, search.best_params_
    else:
        fitted = estimator.fit(x_train, y_train)
        best = {'tuning': 'skipped_insufficient_valid_cv_folds'}
    scores = fitted.predict_proba(x_test)[:, 1]
    per_molecule = molecule_metrics(test_rows, scores)
    return {
        'model': model_name, 'seed': int(seed), 'feature_names': names,
        'n_train_atoms': len(train_rows), 'n_test_atoms': len(test_rows),
        'n_train_molecules': len({row['record_id'] for row in train_rows}),
        'n_test_molecules': len({row['record_id'] for row in test_rows}),
        'dropped_all_missing_features': dropped,
        'best_params': best, 'per_molecule': per_molecule,
        'metrics': aggregate_molecule_metrics(per_molecule),
    }


def run_benchmark(archive, *, groups=('B', 'C', 'D', 'E', 'F', 'G'),
                  model_names=('logistic_regression', 'gradient_boosting'),
                  seeds=(20261106, 20261107, 20261108), bootstrap_replicates=2000):
    """Run the isolated benchmark and return a JSON-serializable result."""
    if isinstance(archive, (str, Path)):
        archive = load_feature_archive(archive)
    results = []
    for group in groups:
        rows = feature_matrix(archive, group)
        for model_name in model_names:
            for seed in seeds:
                results.append({'group': group, **fit_and_evaluate(
                    rows, model_name=model_name, seed=seed)})
    baseline = [row for row in results if row['group'] == 'B']
    for row in results:
        matching = [item for item in baseline if item['model'] == row['model'] and item['seed'] == row['seed']]
        row['bootstrap_vs_group_B'] = (paired_bootstrap(matching[0]['per_molecule'], row['per_molecule'],
                                                        replicates=bootstrap_replicates, seed=row['seed'])
                                       if matching and row['group'] != 'B' else None)
    return {'protocol_id': archive['protocol_id'], 'protocol_sha256': archive['protocol_sha256'],
            'results': results}

