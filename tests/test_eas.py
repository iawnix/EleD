import json

import pytest

from rdkit_eth.eas import (normalize_eas_records, assign_eas_splits, assign_eas_folds,
                           prepare_eas_dataset, frozen_protocol_sha256)


def test_normalize_eas_records_keeps_atom_labels_and_rejects_bad_rows(tmp_path):
    rows, excluded = normalize_eas_records([
        {'id': 'a', 'smiles': 'c1ccccc1', 'reactive_atoms': [0]},
        {'id': 'bad', 'smiles': 'not-smiles', 'reactive_atoms': [0]},
        {'id': 'missing', 'smiles': 'CC'},
        {'id': 'out', 'smiles': 'CC', 'reactive_atoms': [3]},
    ], source='test')
    assert len(rows) == 1 and rows[0]['reactive_atoms'] == [0]
    assert len(excluded) == 3
    assert rows[0]['source'] == 'test' and rows[0]['scaffold']
    mapped, excluded = normalize_eas_records([
        {'id': 'mapped', 'smiles': '[CH3:10][OH:20]', 'reactive_atom_maps': [20]},
    ])
    assert not excluded and mapped[0]['reactive_atoms'] == [1]


def test_eas_splits_keep_scaffold_and_cluster_components_together():
    rows, excluded = normalize_eas_records([
        {'id': str(i), 'smiles': 'c1ccccc1' if i < 2 else 'CCO', 'reactive_atoms': [0]}
        for i in range(6)
    ])
    assert not excluded
    splits = assign_eas_splits(rows, seed=7)
    assert len(set(splits[r['id']] for r in rows if r['scaffold'] == rows[0]['scaffold'])) == 1
    assert set(splits.values()) <= {'train', 'valid', 'test'}
    assert splits == assign_eas_splits(rows, seed=7)
    folds = assign_eas_folds(rows, seed=7, n_splits=5)
    assert set(folds) == {r['id'] for r in rows}
    assert set(folds.values()) <= set(range(5))
    assert folds == assign_eas_folds(rows, seed=7, n_splits=5)


def test_prepare_eas_dataset_jsonl_manifest(tmp_path):
    source = tmp_path / 'eas.jsonl'
    source.write_text(json.dumps({'id': 'x', 'smiles': 'c1ccccc1', 'target_atoms': [1]}) + '\n')
    report = prepare_eas_dataset(source, source='RegioML-test')
    assert report['protocol_id'] == 'eas_regioselectivity_v1'
    assert len(report['protocol_sha256']) == 64
    assert report['protocol_sha256'] == frozen_protocol_sha256()
    assert report['records'][0]['cv_fold'] in range(5)
    assert report['n_selected'] == 1 and report['n_excluded'] == 0
    assert report['records'][0]['split'] in {'train', 'valid', 'test'}
