"""Prepare checked explicit-H conformers, export SDF/report, optionally run xTB."""
import argparse
import json
from pathlib import Path

from rdkit import Chem
from rdkit_eth import XTBCalculator, prepare_conformers


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/prepared-geometry'))
    parser.add_argument('--xtb', action='store_true')
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    prepared = prepare_conformers('C[C@H](O)C(=O)O', num_conformers=4, random_seed=42)
    (args.output/'preparation.json').write_text(json.dumps(prepared.report, indent=2, allow_nan=False)+'\n')
    if not prepared.success: raise RuntimeError(prepared.report['error'])
    with Chem.SDWriter(str(args.output/'conformers.sdf')) as writer:
        for conf in prepared.molecule.GetConformers(): writer.write(prepared.molecule, confId=conf.GetId())
    print('Accepted:', prepared.report['n_conformers'], 'status:', prepared.report['status'])
    if args.xtb:
        result = XTBCalculator().calculate(
            prepared.molecule, conf_id=prepared.report['recommended_conformer_id'])
        print('GFN2 single-point energy (Hartree):', result.total_energy_hartree)


if __name__ == '__main__': main()
