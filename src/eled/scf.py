"""Interpret optional engine SCC diagnostics without inventing convergence data."""
import math
from .errors import XTBEnvironmentError


def scf_record(method, iterations=None, values=None):
    if method == 'GFN0-xTB':
        return {'available': False, 'reason': 'not_applicable_non_scf_method'}
    if iterations is None or iterations < 0:
        return {'available': False, 'reason': 'not_exposed_by_selected_libxtb'}
    if (iterations < 1 or values is None or len(values) != 4 or
            any(not math.isfinite(v) or v < 0 for v in values) or min(values[2:]) <= 0):
        raise XTBEnvironmentError('invalid SCC diagnostics from engine')
    return {'available': True, 'api_version': 1, 'iterations': iterations,
            'energy_change_hartree': values[0], 'mixer_rms': values[1],
            'energy_threshold_hartree': values[2], 'mixer_threshold': values[3],
            'criteria_satisfied': values[0] < values[2] and values[1] < values[3],
            'mixer_definition': 'RMS of SCC mixed variables; GFN2 includes shell charges and multipoles',
            'scope': 'final successful SCC iteration; not a physical accuracy estimate'}
