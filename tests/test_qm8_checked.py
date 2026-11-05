"""Checked benchmark must re-audit geometry and retain split memberships."""
import importlib.util
from pathlib import Path
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem


@pytest.fixture
def study():
    spec = importlib.util.spec_from_file_location('checked', Path(__file__).resolve().parents[1]/'tools/benchmark_qm8_checked.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_reflection_is_rejected_with_unchanged_input_tags(study):
    smiles = 'N[C@@H](C)C(=O)O'
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles)); assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    coords = mol.GetConformer().GetPositions()
    record = dict(smiles=smiles, atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()], positions_angstrom=coords.tolist())
    assert study.audit_record(record)[1]['ok']
    coords[:, 0] *= -1; record['positions_angstrom'] = coords.tolist()
    assert study.audit_record(record)[1]['reason'] == 'stereo_mismatch'


def test_rejected_training_and_test_rows_are_both_removed(study):
    rows = [dict(id=str(i), identity=str(i), scaffold=str(i)) for i in (1, 3)]
    original = [dict(kind='random', seed=17, train=['1','2'], test=['3','4'])]
    restricted = study.restrict_splits(rows, original)
    assert restricted[0]['train'] == ['1'] and restricted[0]['test'] == ['3']
    assert original[0]['train'] == ['1','2']
