"""Export descriptor tables; optionally inspect explicit temperature sensitivity.

Requires the standard extended libxtb installation. Generated geometries are
unoptimized protocol examples, not a prediction benchmark or reference dataset.
"""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (XTBOptions, RunDescriptorBatch, RunReactivity,
                       ExportDescriptorReports)


def molecule(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if AllChem.EmbedMolecule(mol, randomSeed=1) != 0:
        raise RuntimeError(f'Embedding failed: {smiles}')
    return mol


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--sensitivity', action='store_true')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    smiles = ['O', 'CCO', '[CH3]']
    mols = [molecule(s) for s in smiles]
    reports = RunDescriptorBatch(mols, record_ids=smiles)
    for kind in ('csv', 'jsonl'):
        ExportDescriptorReports(reports, args.output / f'single-points.{kind}')
    water = RunReactivity(mols[0], record_id='water', population_scheme='lowdin')
    ExportDescriptorReports([water], args.output / 'water-reactivity.jsonl')
    if args.sensitivity:
        # Each temperature is an explicit independent protocol, never a retry.
        # Retain failing rows and the exact atom coordinates for reproducibility.
        cases = [('O', 'gfn2'), ('CCO', 'gfn1'), ('[NH4+]', 'gfn2')]
        rows, inputs = [], []
        for smi, method in cases:
            mol = molecule(smi)
            inputs.append({'smiles': smi, 'method': method,
                'atomic_numbers': [a.GetAtomicNum() for a in mol.GetAtoms()],
                'coordinates_angstrom': mol.GetConformer().GetPositions().tolist()})
            for temperature in (0., 300., 1000.):
                rows.append(RunReactivity(mol,
                    options=XTBOptions(method=method, accuracy=.1,
                                       electronic_temperature=temperature),
                    record_id=f'{smi}/{method}/{temperature}K', return_errors=True))
        (args.output / 'sensitivity-inputs.json').write_text(json.dumps(inputs, indent=2) + '\n')
        for kind in ('csv', 'jsonl'):
            ExportDescriptorReports(rows, args.output / f'sensitivity.{kind}')
        print(json.dumps([{'id': r['record_id'], 'quality': r['quality']['status'],
                           'values': r['values'], 'error': r['error']} for r in rows], indent=2))
    print(f'Wrote reports to {args.output.resolve()}')


if __name__ == '__main__':
    main()
