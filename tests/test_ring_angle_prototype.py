"""Optional native tests for the research-only bounds patch, not the package API."""
import json
import os
from pathlib import Path
import subprocess

import pytest
from rdkit import Chem


@pytest.fixture(scope='module')
def probes():
    names = ('ELED_RING_ANGLE_CONTROL', 'ELED_RING_ANGLE_CANDIDATE')
    paths = [os.environ.get(n) for n in names]
    if not any(paths):
        pytest.skip('requires separately compiled research bounds probes')
    if not all(paths) or not all(Path(p).is_file() for p in paths):
        pytest.fail('both research probe paths must point to existing executables')
    return paths


def run(binary, smiles, force_field='UFF'):
    result = subprocess.run([binary, smiles, 'ETKDGv4', 'dg', '42', '20', force_field],
                            check=True, capture_output=True, text=True, timeout=30)
    return [json.loads(line) for line in result.stdout.splitlines()]


@pytest.mark.parametrize('smiles,pair', [
    ('N#C[C@@H]1C(=O)N2C[C@H]12', (3, 6)),
    ('O=C1N2C[C@H]2[C@@]12CO2', (1, 3)),
    ('O=C1[C@H]2C=C[C@@H]3[C@H]2N13', (1, 5)),
])
def test_local_angle_changes_bounds_without_changing_stereo_constraints(probes, smiles, pair):
    original, patched = [run(p, smiles)[0] for p in probes]
    before = {tuple(b[:2]): b[2:] for b in original['bounds']}
    after = {tuple(b[:2]): b[2:] for b in patched['bounds']}
    assert after[pair][0] < before[pair][0] - .2
    assert original['chiral_sets'] == patched['chiral_sets']
    assert original['tetrahedral_sets'] == patched['tetrahedral_sets']
    assert original['torsion_atoms'] == patched['torsion_atoms']
    # The source patch precedes downstream bound construction and smoothing.
    # Changes are not confined to overwriting the singled-out matrix entry.
    assert sum(before[p] != after[p] for p in before) > 1


@pytest.mark.parametrize('smiles', [
    'CCO', 'CC(=O)N(C)C', 'O=C1CCN1', 'O=C1CN1', 'c1ccncc1',
])
def test_ordinary_molecules_and_single_small_rings_are_unchanged(probes, smiles):
    assert run(probes[0], smiles) == run(probes[1], smiles)


def test_existing_global_mmff_path_is_unchanged(probes):
    smiles = 'N#C[C@@H]1C(=O)N2C[C@H]12'
    assert run(probes[0], smiles, 'MMFF') == run(probes[1], smiles, 'MMFF')


def test_missing_mmff_parameters_preserve_original_rule(probes):
    from rdkit.Chem import AllChem
    smiles = 'B[C@@H]1C(=O)N2C[C@H]12'
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert not AllChem.MMFFHasAllMoleculeParams(mol)
    assert run(probes[0], smiles) == run(probes[1], smiles)


def test_planarity_ablation_retains_bounds_and_chirality_checks(probes):
    combined = os.environ.get('ELED_RING_ANGLE_COMBINED')
    if not combined:
        pytest.skip('requires the additional diagnostic planarity-ablation build')
    from eled.geometry import _check_geometry, _stereo_constraints
    smiles = 'N#C[C@@H]1C(=O)N2C[C@H]12'
    angle_only = run(probes[1], smiles)
    both = run(combined, smiles)
    assert both[0] == angle_only[0]  # no distance or signed chiral-volume relaxation
    assert both[-1]['conformer_id'] == 0
    accepted = next(s for s in both if s.get('stage') == 'ACCEPTED')
    assert all(c['center_inside'] for c in accepted['chiral_checks'])
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for i, xyz in enumerate(accepted['positions']):
        conf.SetAtomPosition(i, xyz[:3])
    mol.AddConformer(conf)
    assert _check_geometry(mol, 0, _stereo_constraints(mol))['ok']
