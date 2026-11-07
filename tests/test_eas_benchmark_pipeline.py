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
            'atoms.partial_charge_e.mean': -0.2,
            'atoms.partial_charge_e.std': 0.01,
            'atoms.wiberg_sum.mean': 1.2,
            'shell_populations.mulliken_s_electron_count_e.mean': 1.0,
            'shell_populations.mulliken_p_electron_count_e.mean': 2.0,
            'shell_populations.lowdin_s_electron_count_e.mean': 1.1,
            'shell_populations.lowdin_p_electron_count_e.mean': 1.9,
        },
        '1': {
            'atoms.partial_charge_e.mean': 0.2,
            'atoms.partial_charge_e.std': 0.02,
            'atoms.wiberg_sum.mean': 1.1,
            'shell_populations.mulliken_s_electron_count_e.mean': 1.2,
            'shell_populations.mulliken_p_electron_count_e.mean': 1.8,
            'shell_populations.lowdin_s_electron_count_e.mean': 1.3,
            'shell_populations.lowdin_p_electron_count_e.mean': 1.7,
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
    assert matrix[0]['features']['gfn2_partial_charge_e'] == -0.2
    assert matrix[0]['features']['gfn2_mulliken_s_electron_count_e'] == 1.0
    assert matrix[0]['missing']['gfn2_mulliken_d_electron_count_e'] is True

