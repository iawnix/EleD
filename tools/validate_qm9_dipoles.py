#!/usr/bin/env python3
"""Target-blind QM9 nine-heavy-atom dipole validation against B3LYP/6-31G(2df,p).

This is a physical observable comparison, not a fitted model. Nine-heavy-atom
molecules cannot overlap the <=8-heavy-atom QM8 pilot. Both share the GDB family.
"""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import heapq
import io
import json
from pathlib import Path
import tarfile
import time

import numpy as np
from rdkit import Chem, rdBase
from eled import XTBCalculator, calc_xtb_descriptors, XTBEnvironmentError, XTBRunError
from eled.input import _xyz_from_mol, _geometry_hash
from eled.persistence import _implementation_identity

SHA256 = '45255048ac6d83ea4b923ecdf7d6fb6dc62bfec5e80fbc5bcfd93a62157a31db'
PROTOCOL = 'qm9-nine-heavy-dipole-v1'


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--size', type=int, default=256)
    args = parser.parse_args()
    if args.size < 2: raise ValueError('size must be >=2')
    if hashlib.sha256(args.archive.read_bytes()).hexdigest() != SHA256: raise ValueError('archive checksum mismatch')
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()): raise ValueError('output must be empty')
    selected, seen, exclusions = [], set(), Counter()
    with tarfile.open(args.archive) as archive:
        labels = {r['mol_id']: float(r['mu']) for r in csv.DictReader(io.TextIOWrapper(archive.extractfile('gdb9.sdf.csv')))}
        with rdBase.BlockLogs():
            for mol in Chem.ForwardSDMolSupplier(archive.extractfile('gdb9.sdf'), removeHs=False):
                if mol is None: exclusions['sanitization_failed'] += 1; continue
                if mol.GetNumHeavyAtoms() != 9: exclusions['not_nine_heavy_atoms'] += 1; continue
                if Chem.GetFormalCharge(mol) or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
                    exclusions['not_neutral_closed_shell'] += 1; continue
                identity = Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=False)
                if identity in seen: exclusions['duplicate_identity'] += 1; continue
                try: _xyz_from_mol(mol, -1)
                except ValueError:
                    exclusions['unsupported_source_geometry'] += 1; continue
                seen.add(identity)
                rank = int(hashlib.sha256((PROTOCOL+':'+identity).encode()).hexdigest(),16)
                item = (-rank, identity, mol)
                if len(selected) < args.size: heapq.heappush(selected, item)
                elif rank < -selected[0][0]: heapq.heapreplace(selected, item)
    chosen = sorted(selected, reverse=True)
    metadata = []
    for rank, identity, mol in chosen:
        identifier = mol.GetProp('_Name').replace(' ', '_')
        # SDF titles and CSV IDs are verified explicitly, never joined by row order.
        if identifier not in labels: raise ValueError('missing target for '+identifier)
        metadata.append({'id': identifier, 'identity': identity, 'reference_dipole_debye': labels[identifier],
            'geometry_sha256': _geometry_hash(_xyz_from_mol(mol,-1)),
            'atomic_numbers': [a.GetAtomicNum() for a in mol.GetAtoms()],
            'positions_angstrom': mol.GetConformer().GetPositions().tolist()})
    manifest = {'protocol': PROTOCOL, 'archive_sha256': SHA256,
        'source_url': 'https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/gdb9.tar.gz',
        'reference_doi': '10.1038/sdata.2014.22', 'reference': 'QM9 mu, B3LYP/6-31G(2df,p), Debye',
        'geometry': 'provided SDF ground-state coordinates, no optimization',
        'eligible_unique_identities': len(seen), 'exclusions': dict(exclusions), 'selected': metadata,
        'methods': ['GFN1-xTB','GFN2-xTB'], 'settings': {'accuracy':1.,'electronic_temperature':300.},
        'selection': 'first N SHA256(protocol:nonstereo identity), no target-dependent selection',
        'implementation': _implementation_identity(),
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'limitation': 'new nonoverlapping molecules and target; same GDB family, not broad chemical-space validation'}
    write(args.output/'manifest.json',manifest)
    records=[];started=time.perf_counter()
    with gzip.open(args.output/'records.jsonl.gz','wt') as out:
        for (_, _, mol), original in zip(chosen, metadata):
            for method in manifest['methods']:
                row={'id':original['id'],'method':method,'reference_debye':original['reference_dipole_debye']}
                try:
                    result=XTBCalculator(method=method).calculate(mol)
                    row.update(predicted_debye=calc_xtb_descriptors(result)['xtb_dipole_norm_debye'],
                               provenance=result.metadata,error=None)
                except XTBEnvironmentError: raise
                except XTBRunError as exc: row.update(predicted_debye=None,error=str(exc))
                records.append(row);out.write(json.dumps(row,allow_nan=False)+'\n');out.flush()
    common={r['id'] for r in metadata if all(x['error'] is None for x in records if x['id']==r['id'])}
    metrics={}; absolute={}
    for method in manifest['methods']:
        rows=[r for r in records if r['method']==method and r['id'] in common]
        error=np.array([r['predicted_debye']-r['reference_debye'] for r in rows])
        absolute[method]=abs(error)
        metrics[method]={'n':len(rows),'mae_debye':float(np.mean(abs(error))),
                         'rmse_debye':float(np.sqrt(np.mean(error**2))),'bias_debye':float(np.mean(error))}
    diff=absolute['GFN2-xTB']-absolute['GFN1-xTB'];rng=np.random.default_rng(20261105)
    draws=rng.integers(0,len(diff),size=(5000,len(diff)))
    summary={'n_selected':len(metadata),'n_common':len(common),'failed_records':[r for r in records if r['error']],
        'metrics':metrics,'delta_mae_gfn2_minus_gfn1_debye':float(diff.mean()),
        'paired_molecule_bootstrap_95_debye':np.quantile(diff[draws].mean(axis=1),[.025,.975]).tolist(),
        'bootstrap_scope':'selected molecules, not scaffold-independent chemical generalization',
        'calculation_seconds':time.perf_counter()-started}
    write(args.output/'summary.json',summary)
    write(args.output/'SHA256SUMS.json',{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(args.output.iterdir()) if p.is_file()})
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
