"""Dataset identity joins and group separation without downloading QM8."""
import hashlib
import importlib.util
import io
from pathlib import Path
import tarfile

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem


pytestmark = pytest.mark.benchmark


@pytest.fixture
def benchmark():
    spec = importlib.util.spec_from_file_location('benchmark_qm8', Path(__file__).resolve().parents[1]/'tools/benchmark_qm8.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


def test_id_join_dedup_and_duplicate_property_names(benchmark, tmp_path, monkeypatch):
    blocks = []
    for identifier, smiles in [(8, 'O'), (2, 'C'), (13, 'O')]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        assert AllChem.EmbedMolecule(mol, randomSeed=4) == 0
        mol.SetProp('_Name', f'gdb {identifier}\tunused-source-properties')
        blocks.append(Chem.MolToMolBlock(mol)+'$$$$\n')
    # Targets are deliberately reordered relative to SDF; PBE0 names repeat.
    properties = 'gdb9_index,E1-CC2,E1-PBE0,E1-PBE0\n2,0.2,1,2\n13,0.9,3,4\n8,0.3,5,6\n'
    path = tmp_path/'gdb8.tar.gz'
    with tarfile.open(path, 'w:gz') as archive:
        for name, text in [('qm8.sdf', ''.join(blocks)), ('qm8.sdf.csv', properties)]:
            data = text.encode(); info = tarfile.TarInfo(name); info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
    monkeypatch.setattr(benchmark, 'SHA256', hashlib.sha256(path.read_bytes()).hexdigest())
    candidates, excluded, n = benchmark.read_source(tmp_path)
    assert n == 3 and len(candidates) == 2
    rows = {row['id']: row for row, mol in candidates}
    assert rows['8']['target_eV'] == pytest.approx(.3 * benchmark.HARTREE_TO_EV)
    assert rows['2']['target_eV'] == pytest.approx(.2 * benchmark.HARTREE_TO_EV)
    assert excluded == [{'id': '13', 'reason': 'duplicate_nonstereo_identity'}]
    with path.open('ab') as handle: handle.write(b'corruption')
    with pytest.raises(ValueError, match='checksum'): benchmark.read_source(tmp_path)


@pytest.mark.parametrize('kind', ['random', 'scaffold'])
def test_frozen_split_no_identity_or_scaffold_leakage(benchmark, kind):
    pytest.importorskip('sklearn')
    rows = [{'identity': f'molecule-{i}', 'scaffold': '' if i < 8 else f'ring-{i//3}'} for i in range(32)]
    for seed in benchmark.SEEDS:
        train, test = benchmark.split_indices(rows, kind, seed)
        again = benchmark.split_indices(rows, kind, seed)
        assert np.array_equal(train, again[0]) and np.array_equal(test, again[1])
        assert set(train).isdisjoint(test) and set(train) | set(test) == set(range(32))
        if kind == 'scaffold':
            assert {rows[i]['scaffold'] for i in train}.isdisjoint({rows[i]['scaffold'] for i in test})
            assert set(range(8)).issubset(train) or set(range(8)).issubset(test)
