"""Paired gates: gains never offset losses of a baseline success."""


def compare_pair(baseline, candidate, energy_tolerance=1e-3):
    failures = []
    if baseline['initial_valid'] and not candidate['initial_valid']:
        failures.append('lost_initial_success')
    if baseline['post_valid'] and not candidate['post_valid']:
        failures.append('lost_post_optimization_success')
    delta = None
    if baseline['post_valid'] and candidate['post_valid']:
        delta = candidate['post_energy_kcal_mol'] - baseline['post_energy_kcal_mol']
        if delta > energy_tolerance:
            failures.append('higher_post_mmff_energy')
    return {'failures': failures, 'post_energy_delta_kcal_mol': delta,
            'gain': not baseline['initial_valid'] and candidate['initial_valid']}


def verdict(comparisons, affected_baseline_successes):
    failures = sum(bool(c['failures']) for c in comparisons)
    if failures:
        status = 'reject_observed_regression'
    elif affected_baseline_successes == 0:
        status = 'inconclusive_no_affected_baseline_success'
    else:
        status = 'no_regression_observed_in_tested_pairs'
    return {'status': status, 'regressed_pairs': failures,
            'affected_baseline_successes': affected_baseline_successes,
            'gained_pairs': sum(c['gain'] for c in comparisons),
            'production_ready': False,
            'limits': 'finite chemistry/seed coverage; diagnostic timings are not a production performance gate'}
