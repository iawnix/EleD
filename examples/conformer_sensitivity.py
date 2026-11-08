"""ETKDG/MMFF geometry preparation followed by native xTB sensitivity analysis."""
import argparse
import gzip
import json
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from eled import XTBDiskCache, RunConformerAnalysis


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=Path('build/conformer-example'))
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    cache = XTBDiskCache(args.output/'single-points.sqlite')
    reports, summary = [], []
    for name, smiles in [('ethanol', 'CCO'), ('ethylene_glycol', 'OCCO'),
                         ('butane', 'CCCC'), ('phenethylamine', 'NCCc1ccccc1'), ('benzene', 'c1ccccc1')]:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        params = AllChem.ETKDGv3(); params.randomSeed = 20261104
        params.pruneRmsThresh = .2; params.numThreads = 1
        ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=12, params=params))
        optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=500)
        retained = [i for i, (status, energy) in zip(ids, optimized) if status == 0]
        if not retained: raise RuntimeError(f'no converged MMFF conformers: {name}')
        report = RunConformerAnalysis(mol, conf_ids=retained, record_id=name, cache=cache)
        report['geometry_preparation'] = {'rdkit': rdBase.rdkitVersion, 'smiles': smiles,
            'embedding': 'ETKDGv3', 'seed': 20261104, 'requested': 12, 'prune_rms_angstrom': .2,
            'optimization': 'MMFF94; maxIters=500; no post-optimization deduplication',
            'optimization_results': [{'conf_id': i, 'status': s, 'energy_kcal_mol': e}
                                     for i, (s, e) in zip(ids, optimized)],
            'retained_converged_ids': retained,
            'atomic_numbers': [a.GetAtomicNum() for a in mol.GetAtoms()],
            'positions_angstrom': {str(i): mol.GetConformer(i).GetPositions().tolist() for i in retained}}
        reports.append(report)
        summary.append({'name': name, 'n_conformers': report['n_conformers'], 'status': report['status'],
                        'gap_eV': report['statistics']['xtb_gap_ev'],
                        'homo_participation_ratio': report['statistics']['frontier_homo_participation_ratio'],
                        'centroid_separation_angstrom': report['statistics']['frontier_homo_lumo_atom_centroid_separation_angstrom']})
    with gzip.open(args.output/'reports.jsonl.gz', 'wt') as handle:
        for report in reports: handle.write(json.dumps(report, allow_nan=False)+'\n')
    text = json.dumps(summary, indent=2, allow_nan=False)+'\n'
    (args.output/'summary.json').write_text(text); print(text)


if __name__ == '__main__': main()
