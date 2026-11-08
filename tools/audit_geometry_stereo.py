#!/usr/bin/env python3
"""Post-hoc coordinate/stereo audit of the immutable 0.9.1 generated geometries."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from rdkit import Chem
from eled.geometry import _check_geometry, _stereo_constraints
from eled.persistence import _implementation_identity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=Path(__file__).resolve().parents[1]/'docs/validation/0.9.1/geometry-preparation.jsonl.gz')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists(): raise ValueError('output already exists')
    rejected, checked = [], 0
    for row in map(json.loads, gzip.open(args.input, 'rt')):
        if row['failure']: continue
        mol = Chem.AddHs(Chem.MolFromSmiles(row['smiles']))
        if [a.GetAtomicNum() for a in mol.GetAtoms()] != row['atomic_numbers']:
            raise ValueError('atom ordering differs from the recorded geometry')
        conf = Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
        if len(row['positions_angstrom']) != mol.GetNumAtoms(): raise ValueError('coordinate count mismatch')
        for i, xyz in enumerate(row['positions_angstrom']): conf.SetAtomPosition(i, xyz)
        mol.AddConformer(conf)
        check = _check_geometry(mol, 0, _stereo_constraints(mol)); checked += 1
        if not check['ok']: rejected.append({'id': row['id'], 'smiles': row['smiles'], 'check': check})
    result = {'scope': 'post-hoc coordinate/stereo audit; old predictions and geometries are not changed',
        'input_sha256': hashlib.sha256(args.input.read_bytes()).hexdigest(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'implementation': _implementation_identity(),
        'n_checked': checked, 'n_rejected': len(rejected), 'rejected': rejected}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(checked, 'checked;', len(rejected), 'rejected')


if __name__ == '__main__': main()
