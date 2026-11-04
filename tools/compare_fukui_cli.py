#!/usr/bin/env python3
"""Explicit developer comparison against xTB --vfukui, never a runtime backend.

Uses fixed coordinates, matched settings and explicit spin choices compatible
with the CLI policy. Keeps full precision adapter results and original CLI logs.
The native CLI prints only three decimals: agreement is bounded by that precision.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import XTBOptions, run_reactivity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--library', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True, help='new empty directory')
    args = parser.parse_args()
    executable = args.executable.resolve(strict=True)
    library = args.library.resolve(strict=True)
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()): parser.error('choose an empty output directory')
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1',
               XTBPATH=str(library.parent.parent / 'share/xtb'))
    os.environ.update({key:env[key] for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','XTBPATH')})
    rows = []
    for smi in ('O', 'CCO', '[CH3]', '[NH4+]'):
        mol = Chem.AddHs(Chem.MolFromSmiles(smi))
        if AllChem.EmbedMolecule(mol, randomSeed=1): raise RuntimeError('Embedding failed')
        charge = Chem.GetFormalCharge(mol)
        spin = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
        for method in (1, 2):
            case = output / f'case-{len(rows):02d}'
            case.mkdir()
            xyz = str(mol.GetNumAtoms())+'\nfixed unoptimized geometry\n'+''.join(
                f'{a.GetSymbol()} {x:.17g} {y:.17g} {z:.17g}\n'
                for a,(x,y,z) in zip(mol.GetAtoms(),mol.GetConformer().GetPositions()))
            (case / 'input.xyz').write_text(xyz)
            command = [str(executable),'input.xyz','--gfn',str(method),'--chrg',str(charge),
                       '--uhf',str(spin),'--etemp','300','--acc','.1','--iterations','250','--vfukui']
            completed = subprocess.run(command, cwd=case, env=env, capture_output=True, text=True, timeout=60)
            (case / 'stdout.log').write_text(completed.stdout)
            (case / 'stderr.log').write_text(completed.stderr)
            row = {'smiles':smi,'method':method,'command':command,'directory':case.name,
                   'returncode':completed.returncode,'cli_output_decimal_places':3,
                   'absolute_tolerance':0.00051,'reports':{}}
            tail = completed.stdout.split('Fukui functions:')[-1].split('\n\n',1)[0]
            matches = re.findall(r'^\s*(\d+)\s*([A-Za-z]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$', tail, re.M)
            valid = (completed.returncode == 0 and 'Fukui functions:' in completed.stdout
                and len(matches) == mol.GetNumAtoms()
                and all(int(item[0])==i+1 and item[1]==mol.GetAtomWithIdx(i).GetSymbol()
                        for i,item in enumerate(matches)))
            if valid:
                cli = np.array([[float(v) for v in item[2:]] for item in matches])
                row['cli_fukui'] = cli.tolist()
            for initialization in ('reference','independent'):
                report = run_reactivity(mol, options=XTBOptions(method=f'gfn{method}',
                    xtb_library=str(library), accuracy=.1, electronic_temperature=300.),
                    ionized_uhf=spin-1 if spin else 1, attached_uhf=spin-1 if spin else 1,
                    initialization=initialization, return_errors=True)
                row['reports'][initialization] = report
            reference = row['reports']['reference']
            if valid and reference['error'] is None:
                computed = np.array([[a[k] for k in ('f_plus','f_minus','f_radical')] for a in reference['atoms']])
                row['max_abs_cli_difference'] = float(np.max(np.abs(computed-cli)))
                row['matched_printed_precision'] = row['max_abs_cli_difference'] <= row['absolute_tolerance']
            else:
                row['matched_printed_precision'] = False
            independent = row['reports']['independent']
            if reference['error'] is None and independent['error'] is None:
                row['initialization_difference'] = {
                    'max_abs_state_energy_hartree': max(abs(reference['states'][k]['energy_hartree']-
                        independent['states'][k]['energy_hartree']) for k in reference['states']),
                    'max_abs_fukui': max(abs(a[k]-b[k]) for a,b in zip(reference['atoms'],independent['atoms'])
                                        for k in ('f_plus','f_minus','f_radical'))}
            rows.append(row)
    evidence = {'executable':str(executable),'executable_sha256':hashlib.sha256(executable.read_bytes()).hexdigest(),
                'library':str(library),'library_sha256':hashlib.sha256(library.read_bytes()).hexdigest(),
                'scope':'four fixed unoptimized molecules, GFN1/2, 300 K; CLI rounded to 3 decimals',
                'matched':sum(r['matched_printed_precision'] for r in rows),'total':len(rows),'cases':rows}
    (output / 'comparison.json').write_text(json.dumps(evidence, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k:v for k,v in evidence.items() if k!='cases'},indent=2))
    return int(evidence['matched'] != evidence['total'])


if __name__ == '__main__':
    raise SystemExit(main())
