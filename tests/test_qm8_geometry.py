"""Matched comparison, generated geometry provenance and clustered uncertainty."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


pytestmark = pytest.mark.benchmark


@pytest.fixture
def study():
    spec = importlib.util.spec_from_file_location('qm8_geometry', Path(__file__).resolve().parents[1]/'tools/benchmark_qm8_geometry.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_smiles_only_geometry_is_deterministic_and_records_selection(study):
    a, pa = study.prepare_geometry('OCCO', 'example')
    b, pb = study.prepare_geometry('OCCO', 'example')
    assert pa == pb and pa['source_coordinates_used'] is False
    assert a.GetNumConformers() == b.GetNumConformers() == 1
    assert np.array_equal(a.GetConformer().GetPositions(), b.GetConformer().GetPositions())
    chosen = min((r['energy_kcal_mol'], r['conf_id']) for r in pa['optimization_results'] if r['status'] == 0)
    assert pa['selected_embedded_id'] == chosen[1]
    assert pa['atomic_numbers'] == [atom.GetAtomicNum() for atom in a.GetAtoms()]
    assert np.array_equal(pa['positions_angstrom'], a.GetConformer().GetPositions())


@pytest.mark.parametrize('failure', ['embedding', 'parameters', 'convergence'])
def test_geometry_failures_are_explicit_without_fallback(study, monkeypatch, failure):
    if failure == 'embedding':
        monkeypatch.setattr(study.AllChem, 'EmbedMultipleConfs', lambda *a, **k: [])
        expected = 'embedding_failed'
    elif failure == 'parameters':
        monkeypatch.setattr(study.AllChem, 'MMFFHasAllMoleculeParams', lambda *a: False)
        expected = 'MMFF_parameters_unavailable'
    else:
        monkeypatch.setattr(study.AllChem, 'MMFFOptimizeMoleculeConfs',
                            lambda mol, **k: [(1, 1.)] * mol.GetNumConformers())
        expected = 'no_converged_MMFF_conformer'
    mol, record = study.prepare_geometry('CCO', 'example')
    assert mol is None and record['failure'] == expected


def test_paired_bootstrap_direction_identity_join_and_constant_effect(study):
    reference = {'b': (1., 3.), 'a': (2., 4.), 'c': (3., 5.)}
    candidate = {'c': (3., 4.), 'b': (1., 2.), 'a': (2., 3.)}
    groups = {'a': 'a', 'b': 'b', 'c': 'c'}
    result = study.paired_interval(reference, candidate, groups=groups, seed=7, replicates=200)
    assert result['delta_mae_eV'] == -1.
    assert result['interval_95_eV'] == [-1., -1.]
    assert study.paired_interval(reference, reference, groups=groups, seed=7)['interval_95_eV'] == [0., 0.]


def test_scaffold_cluster_resampling_and_single_group_limitation(study):
    ref = {str(i): (0., 2.) for i in range(10)}
    cand = {str(i): (0., 1. if i < 5 else 3.) for i in range(10)}
    groups = {str(i): '' if i < 5 else 'ring' for i in range(10)}
    clustered = study.paired_interval(ref, cand, groups=groups, seed=13, replicates=2000)
    assert clustered['n_groups'] == 2 and clustered['interval_95_eV'] == [-1., 1.]
    iid = study.paired_interval(ref, cand, groups={i:i for i in ref}, seed=13, replicates=2000)
    assert iid['interval_95_eV'][0] > clustered['interval_95_eV'][0]
    single = study.paired_interval(ref, cand, groups={i:'' for i in ref}, seed=13)
    assert single['interval_95_eV'] is None and single['warning']


@pytest.mark.parametrize('candidate,groups', [({'x': (0., 1.)}, {'a':'a'}),
    ({'a': (1., 0.)}, {'a':'a'}), ({'a': (0., float('nan'))}, {'a':'a'}),
    ({'a': (0., 0.)}, {'b':'b'})])
def test_paired_bootstrap_rejects_mismatched_or_invalid_data(study, candidate, groups):
    with pytest.raises(ValueError):
        study.paired_interval({'a': (0., 1.)}, candidate, groups=groups, seed=1)


def test_frozen_evidence_integrity(study, tmp_path):
    root = Path(__file__).resolve().parents[1]/'docs/validation/0.9.0'
    manifest, rows, splits, hashes = study.read_evidence(root)
    assert len(rows) == 1024 and len(splits) == 6 and manifest['protocol'] == 'qm8-pilot-v1'
    (tmp_path/'SHA256SUMS.json').write_text((root/'SHA256SUMS.json').read_text())
    first = next(iter(hashes)); (tmp_path/first).write_bytes(b'corrupt')
    with pytest.raises(ValueError, match='checksum'): study.read_evidence(tmp_path)


def test_frozen_splits_do_not_resample_or_train_on_test_targets(study, tmp_path, monkeypatch):
    pytest.importorskip('sklearn')
    rows = [{'id': str(i), 'identity': str(i), 'scaffold': str(i), 'rdkit': [float(i)],
             'basic': [float(i*i)], 'orbital': [float(i+1)], 'target_eV': float(i)} for i in range(8)]
    split = [{'kind': 'random', 'seed': 17, 'train': [str(i) for i in range(6)], 'test': ['6', '7']}]
    def forbidden(*a): raise AssertionError('must use frozen partition')
    monkeypatch.setattr(study.pilot, 'split_indices', forbidden)
    first = tmp_path/'a'; second = tmp_path/'b'; first.mkdir(); second.mkdir()
    study.pilot.fit_ablation(rows, first, splits=split)
    rows[6]['target_eV'] += 100.; rows[7]['target_eV'] += 100.
    study.pilot.fit_ablation(rows, second, splits=split)
    a = json.loads((first/'predictions.json').read_text()); b = json.loads((second/'predictions.json').read_text())
    assert [r['prediction_eV'] for r in a] == [r['prediction_eV'] for r in b]
    assert json.loads((first/'splits.json').read_text()) == split


@pytest.mark.parametrize('mutation', ['overlap', 'missing', 'duplicate', 'group_leak', 'empty', 'duplicate_split'])
def test_invalid_frozen_splits_rejected(study, mutation):
    rows = [{'id': str(i), 'identity': str(i), 'scaffold': str(i)} for i in range(4)]
    splits = [{'kind': 'scaffold', 'seed': 1, 'train': ['0', '1'], 'test': ['2', '3']}]
    if mutation == 'overlap': splits[0]['train'].append('2')
    elif mutation == 'missing': splits[0]['test'].pop()
    elif mutation == 'duplicate': splits[0]['train'].append('0')
    elif mutation == 'group_leak': rows[2]['scaffold'] = rows[0]['scaffold'] = ''
    elif mutation == 'empty': splits[0]['test'] = []
    else: splits.append(deepcopy(splits[0]))
    with pytest.raises(ValueError): study.pilot.validate_splits(rows, splits)
