#!/usr/bin/env python3
"""Compare two native builds on frozen geometries, in separate Python processes.

Reports numerical changes rather than asserting equality across upstream model
changes. Both libraries must carry reviewed GFN0 and basis capability manifests.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]


def calculate(library):
    import numpy as np
    from rdkit import Chem
    from eled import XTBCalculator, calc_orbital_populations
    from eled.runtime import _fingerprint
    cases = json.loads((ROOT / 'tests/data/orbital_alignment_cases.json').read_text())['cases']
    rows = []
    for case in cases:
        mol = Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        conf = Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
        for i, xyz in enumerate(case['positions_angstrom']): conf.SetAtomPosition(i, xyz)
        mol.AddConformer(conf)
        for method in ('GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'):
            result = XTBCalculator(method=method, xtb_library=str(library)).calculate(
                mol, keep_overlap=True)
            population = calc_orbital_populations(result, scheme='mulliken')
            rows.append({'name': case['name'], 'method': method, 'energy_eh': result.total_energy_hartree,
                         'orbital_energies_ev': result.orbital_energies_ev,
                         'charges_e': result.partial_charges, 'dipole_e_bohr': result.dipole_e_bohr,
                         'wiberg': result.wiberg_bond_orders,
                         'orthonormality_error': population['orthonormality_error'],
                         'electron_error': float(abs(sum(population['electron_populations_e']) - result.number_of_electrons))})
    return {'library': str(library.resolve()), 'sha256': _fingerprint(library), 'rows': rows}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--reference', type=Path)
    p.add_argument('--candidate', type=Path)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--worker', type=Path, help=argparse.SUPPRESS)
    args = p.parse_args()
    if args.worker:
        args.output.write_text(json.dumps(calculate(args.worker), indent=2) + '\n')
        return
    if args.reference is None or args.candidate is None: p.error('provide both library paths')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    builds = []
    for label, lib in [('reference', args.reference), ('candidate', args.candidate)]:
        lib = lib.resolve(strict=True)
        output = args.output.with_name(args.output.stem + '-' + label + '.json')
        env = os.environ.copy()
        for name in ('XTB_LIBRARY', 'XTBPATH', 'XTBHOME', 'XTB_PREFIX', 'LD_LIBRARY_PATH', 'LD_PRELOAD'):
            env.pop(name, None)
        env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
        subprocess.run([sys.executable, str(Path(__file__).resolve()), '--worker', str(lib),
                        '--output', str(output)], env=env, check=True)
        builds.append(json.loads(output.read_text()))
    import numpy as np
    records = []
    for old, new in zip(builds[0]['rows'], builds[1]['rows'], strict=True):
        assert (old['name'], old['method']) == (new['name'], new['method'])
        records.append({'name': old['name'], 'method': old['method'],
                        'delta_energy_eh': new['energy_eh'] - old['energy_eh'],
                        **{'max_abs_' + field: float(np.max(np.abs(np.array(new[field])-old[field])))
                           for field in ('orbital_energies_ev', 'charges_e', 'dipole_e_bohr', 'wiberg')}})
    summary = {}
    for method in ('GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'):
        rows = [r for r in records if r['method'] == method]
        summary[method] = {'cases': len(rows), **{key: max(abs(r[key]) for r in rows)
                                                 for key in rows[0] if key not in ('name','method')}}
    report = {'reference': {k:v for k,v in builds[0].items() if k != 'rows'},
              'candidate': {k:v for k,v in builds[1].items() if k != 'rows'},
              'comparison': 'frozen geometry; separate processes; candidate minus reference',
              'summary_max_abs': summary, 'cases': records}
    args.output.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__': main()
