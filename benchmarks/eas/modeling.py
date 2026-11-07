"""Model-independent ranking metrics and paired molecule bootstrap."""
from __future__ import annotations

from collections import defaultdict
import math
from numbers import Integral

import numpy as np


def molecule_metrics(rows, scores):
    """Compute per-molecule EAS ranking metrics from atom scores.

    Molecules without both a positive and negative candidate are returned with
    ``None`` metrics and excluded from aggregate means. This is explicit because
    ROC-AUC and a meaningful ranking cannot be defined for such a molecule.
    """
    if len(rows) != len(scores):
        raise ValueError('rows and scores must have the same length')
    grouped = defaultdict(list)
    for row, score in zip(rows, scores):
        if not math.isfinite(float(score)):
            raise ValueError('scores must be finite')
        grouped[row['record_id']].append((int(row['target']), float(score)))
    output = {}
    for identifier, values in grouped.items():
        positives = sum(target for target, _ in values)
        if positives == 0 or positives == len(values):
            output[identifier] = {'macro_roc_auc': None, 'mrr': None,
                                  'top_1': None, 'top_2': None,
                                  'n_candidates': len(values), 'n_positive': positives}
            continue
        ordered = sorted(values, key=lambda item: (-item[1], -item[0]))
        ranks = [index + 1 for index, (target, _) in enumerate(ordered) if target]
        # AUC is the fraction of positive/negative pairs in the correct order;
        # ties receive half credit.
        pos_scores = [score for target, score in values if target]
        neg_scores = [score for target, score in values if not target]
        pair_score = sum(1.0 if p > n else 0.5 if p == n else 0.0
                         for p in pos_scores for n in neg_scores)
        auc = pair_score / (len(pos_scores) * len(neg_scores))
        output[identifier] = {
            'macro_roc_auc': float(auc), 'mrr': float(1.0 / ranks[0]),
            'top_1': float(any(target for target, _ in ordered[:1])),
            'top_2': float(any(target for target, _ in ordered[:2])),
            'n_candidates': len(values), 'n_positive': positives,
        }
    return output


def aggregate_molecule_metrics(per_molecule):
    """Average defined per-molecule metrics with equal molecule weight."""
    names = ('macro_roc_auc', 'mrr', 'top_1', 'top_2')
    result = {}
    for name in names:
        values = [row[name] for row in per_molecule.values() if row.get(name) is not None]
        result[name] = float(np.mean(values)) if values else None
        result[f'{name}_n'] = len(values)
    result['n_molecules'] = len(per_molecule)
    return result


def paired_bootstrap(baseline, candidate, *, replicates=2000, seed=20261106):
    """Return paired candidate-minus-baseline molecule bootstrap intervals."""
    if isinstance(replicates, bool) or not isinstance(replicates, Integral) or replicates <= 0:
        raise ValueError('replicates must be a positive integer')
    if set(baseline) != set(candidate):
        raise ValueError('baseline and candidate must contain identical molecule IDs')
    names = ('macro_roc_auc', 'mrr', 'top_1', 'top_2')
    rng = np.random.default_rng(seed)
    identifiers = list(baseline)
    deltas = {name: [] for name in names}
    for _ in range(int(replicates)):
        selected = rng.integers(0, len(identifiers), size=len(identifiers))
        for name in names:
            values = []
            for index in selected:
                b, c = baseline[identifiers[index]].get(name), candidate[identifiers[index]].get(name)
                if b is not None and c is not None:
                    values.append(float(c) - float(b))
            deltas[name].append(float(np.mean(values)) if values else np.nan)
    output = {}
    for name in names:
        values = np.asarray(deltas[name], dtype=float)
        valid = values[np.isfinite(values)]
        output[name] = {
            'delta': float(np.mean(valid)) if len(valid) else None,
            'ci95': ([float(np.percentile(valid, 2.5)), float(np.percentile(valid, 97.5))]
                     if len(valid) else None),
            'n_bootstrap': int(len(valid)),
        }
    return output

