import json
from pathlib import Path
import hashlib
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from benchmarks.eas.pipeline import (atom_rows, feature_matrix,
                                     validate_feature_archive)


def _protocol_hash():
    payload = json.loads(Path('protocol.json').read_text())
    return hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(',', ':'),
                                     ensure_ascii=False).encode()).hexdigest()


def _archive():
    atom_stats = {
        '0': {
            'atoms.q_e': {'mean': -0.2, 'std': 0.01, 'min': -0.22, 'max': -0.18},
            'atoms.wbo_sum': {'mean': 1.2},
            'shell_populations.mul_s_e': {'mean': 1.0},
            'shell_populations.mul_p_e': {'mean': 2.0},
            'shell_populations.low_s_e': {'mean': 1.1},
            'shell_populations.low_p_e': {'mean': 1.9},
        },
        '1': {
            'atoms.q_e': {'mean': 0.2, 'std': 0.02},
            'atoms.wbo_sum': {'mean': 1.1},
            'shell_populations.mul_s_e': {'mean': 1.2},
            'shell_populations.mul_p_e': {'mean': 1.8},
            'shell_populations.low_s_e': {'mean': 1.3},
            'shell_populations.low_p_e': {'mean': 1.7},
        },
    }
    analysis = {'status': 'complete', 'atom_statistics': atom_stats, 'reports': []}
    return {
        'protocol_id': 'eas_regioselectivity_v1',
        'protocol_sha256': _protocol_hash(),
        'records': [{
            'id': 'm1', 'smiles': 'C=C', 'scaffold': 'C=C',
            'reactive_atoms': [1], 'split': 'train', 'cv_fold': 0,
            'atom_features': [
                {'atomic_number': 6, 'formal_charge': 0, 'degree': 1,
                 'total_hydrogens': 2, 'is_aromatic': 0, 'is_in_ring': 0,
                 'hybridization': 'SP2'},
                {'atomic_number': 6, 'formal_charge': 0, 'degree': 1,
                 'total_hydrogens': 2, 'is_aromatic': 0, 'is_in_ring': 0,
                 'hybridization': 'SP2'},
            ],
            'analyses': {
                'gfn1_300K': analysis,
                'gfn2_300K': analysis,
                'gfn2_1000K': analysis,
                'gfn2_reactivity_300K': [],
            },
        }],
    }


def test_feature_archive_validation_and_atom_rows():
    archive = _archive()
    check = validate_feature_archive(archive)
    assert check['n_records'] == 1
    rows = atom_rows(archive)
    assert [row['target'] for row in rows] == [0, 1]
    assert [row['atom_index'] for row in rows] == [0, 1]


def test_feature_groups_keep_missing_values_explicit():
    matrix = feature_matrix(_archive(), 'D')
    assert matrix[0]['features']['rdkit_atomic_number'] == 6
    assert matrix[0]['features']['gfn2_q_e'] == -0.2
    assert matrix[0]['features']['gfn2_mul_s_e'] == 1.0
    assert matrix[0]['missing']['gfn2_mul_d_e'] is True


def test_conformer_summary_fields_use_nested_statistics():
    import pytest
    matrix = feature_matrix(_archive(), 'G')
    assert matrix[0]['features']['gfn2_q_e_std'] == .01
    assert matrix[0]['features']['gfn2_q_e_range'] == pytest.approx(.04)
    assert matrix[1]['features']['gfn2_q_e_range'] is None


def test_bond_features_read_nested_sections_and_keep_frontier_groups_separate():
    from copy import deepcopy
    import pytest
    archive = _archive()
    pair = {'bond_index': 0, 'begin_atom_index': 0, 'end_atom_index': 1}
    report = {
        'error': None, 'quality': {'status': 'ok'},
        'bonds': [dict(pair, wbo=1.2)],
        'mayer_bonds': {'bonds': [dict(pair, mbo=1.1)]},
        'frontier_bonds': {'bonds': [dict(pair, frontier_group='homo', mul_ovlp_e=.3),
                                    dict(pair, frontier_group='lumo', mul_ovlp_e=-.2)]},
        'pi': {'bonds': [dict(pair, pi_bo=.6)]},
    }
    other = deepcopy(report)
    other['bonds'][0]['wbo'] = 1.4
    other['mayer_bonds']['bonds'][0]['mbo'] = None
    other['pi']['bonds'][0]['pi_bo'] = 0.  # valid zero, not missing
    failed = deepcopy(report); failed['error'] = {'message': 'failed'}
    failed['bonds'][0]['wbo'] = 100.
    quality_failed = deepcopy(failed); quality_failed['error'] = None
    quality_failed['quality']['status'] = 'error'
    archive['records'][0]['analyses']['gfn2_300K']['reports'] = [report, other, failed, quality_failed]
    for group in ('E', 'F'):
        rows = feature_matrix(archive, group)
        for row in rows:
            features = row['features']
            for key, expected in {'wbo': 1.3, 'mbo': 1.1, 'homo_mul_ovlp_e': .3,
                                  'lumo_mul_ovlp_e': -.2, 'pi_bo': .3}.items():
                assert features['gfn2_'+key] == pytest.approx(expected)
                assert not row['missing']['gfn2_'+key]
            assert 'gfn2_mul_ovlp_e' not in features
    report['frontier_bonds'] = other['frontier_bonds'] = None
    assert feature_matrix(archive, 'E')[0]['features']['gfn2_homo_mul_ovlp_e'] is None


def test_real_report_shape_flows_into_benchmark_bond_features():
    import pytest
    from eled import CalcDescriptorReport
    from test_descriptor_regressions import pi_system
    mol, result = pi_system()
    report = CalcDescriptorReport(mol, result)
    archive = _archive()
    archive['records'][0]['analyses']['gfn2_300K']['reports'] = [report]
    row = feature_matrix(archive, 'E')[0]
    for name, section, field, group in [
        ('mbo', 'mayer_bonds', 'mbo', None), ('pi_bo', 'pi', 'pi_bo', None),
        ('homo_mul_ovlp_e', 'frontier_bonds', 'mul_ovlp_e', 'homo'),
        ('lumo_mul_ovlp_e', 'frontier_bonds', 'mul_ovlp_e', 'lumo'),
    ]:
        values = [r[field] for r in report[section]['bonds']
                  if 0 in (r['begin_atom_index'], r['end_atom_index'])
                  and (group is None or r['frontier_group'] == group) and r[field] is not None]
        assert values
        assert row['features']['gfn2_'+name] == pytest.approx(sum(values)/len(values))
