#!/usr/bin/env python
"""Generate the frozen EAS feature triplet for a prepared manifest.

The manifest must come from ``PrepareEASDataset.py``. Failures are retained
per molecule and are never converted into zero-valued atom features.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from rdkit import Chem
from eled import (XTBOptions, PrepareConformers, RunConformerAnalysis,
                       RunReactivity)
from eled.eas import FrozenProtocolSHA256, LoadFrozenProtocol


def _atom_features(smiles):
    """Serialize stable RDKit graph features for an external study archive."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('SMILES could not be parsed for atom feature export')
    values = []
    for atom in mol.GetAtoms():
        values.append({
            'atomic_number': atom.GetAtomicNum(),
            'formal_charge': atom.GetFormalCharge(),
            'degree': atom.GetDegree(),
            'total_hydrogens': atom.GetTotalNumHs(),
            'is_aromatic': int(atom.GetIsAromatic()),
            'is_in_ring': int(atom.IsInRing()),
            'hybridization': str(atom.GetHybridization()),
        })
    return values

def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('manifest', type=Path)
    parser.add_argument('-o', '--output', type=Path, required=True)
    parser.add_argument('--xtb-library', type=Path)
    parser.add_argument('--num-conformers', type=int, default=None)
    parser.add_argument('--seed', type=int, default=None)
    args = parser.parse_args(argv)
    protocol = LoadFrozenProtocol()
    requested_conformers = protocol['geometry']['conformers_per_molecule']
    if args.num_conformers is not None and args.num_conformers != requested_conformers:
        raise SystemExit(f'--num-conformers must equal frozen protocol value {requested_conformers}')
    args.num_conformers = requested_conformers
    protocol_seed = protocol['geometry']['seed']
    if args.seed is not None and args.seed != protocol_seed:
        raise SystemExit(f'--seed must equal frozen protocol value {protocol_seed}')
    args.seed = protocol_seed
    manifest = json.loads(args.manifest.read_text(encoding='utf-8'))
    if (manifest.get('protocol_id') != protocol['protocol_id'] or
            manifest.get('protocol_sha256') != FrozenProtocolSHA256(protocol)):
        raise SystemExit('manifest was not produced by the frozen EAS protocol')
    records = []
    for row in manifest.get('records', []):
        digest = hashlib.sha256(f"{args.seed}:{row['id']}".encode()).hexdigest()
        seed = int(digest[:8], 16) % 2147483646 + 1
        item = {'id': row['id'], 'smiles': row['smiles'],
                'scaffold': row['scaffold'], 'reactive_atoms': row['reactive_atoms'],
                'split': row['split'], 'cv_fold': row['cv_fold'],
                'charge': row.get('charge', 0), 'uhf': row.get('uhf', 0),
                'atom_features': _atom_features(row['smiles'])}
        try:
            prepared = PrepareConformers(row['smiles'], numConformers=args.num_conformers,
                method=protocol['geometry']['embedding']['method'], fallbackMethods=(),
                attemptsPerStrategy=protocol['geometry']['embedding']['attempts_per_strategy'],
                retryRandomCoords=protocol['geometry']['embedding']['retry_random_coords'],
                maxEmbeddingIterations=protocol['geometry']['embedding']['max_embedding_iterations'],
                optimize='xTB', xtbMethod='GFN2-xTB', xtbLibrary=args.xtb_library,
                optimizationIterations=protocol['geometry']['optimization']['max_iterations'],
                randomSeed=seed, charge=row.get('charge'), numUnpairedElectrons=row.get('uhf'))
            item['preparation'] = prepared.report
            if not prepared.success:
                item['status'] = 'preparation_failed'; records.append(item); continue
            common = dict(xtbLibrary=args.xtb_library, keepOverlap=True,
                          electronicTemperature=protocol['calculation']['electronic_temperature_K'],
                          charge=row.get('charge'), numUnpairedElectrons=row.get('uhf'))
            primary = RunConformerAnalysis(prepared.molecule,
                options=XTBOptions(method='GFN2-xTB', **common), recordId=row['id'],
                returnErrors=True)
            gfn1 = RunConformerAnalysis(prepared.molecule,
                options=XTBOptions(method='GFN1-xTB', **common), recordId=row['id'],
                returnErrors=True)
            hot_common = {**common, 'electronicTemperature':
                          protocol['calculation']['diagnostic_electronic_temperature_K']}
            hot = RunConformerAnalysis(prepared.molecule,
                options=XTBOptions(method='GFN2-xTB', **hot_common),
                recordId=row['id'], returnErrors=True)
            reactivity = []
            valid_conformer_ids = [int(report['record_id']) for report in primary.get('reports', [])
                                   if report.get('quality', {}).get('status') != 'error'
                                   and report.get('error') is None]
            for conformer_id in valid_conformer_ids:
                reactivity.append(RunReactivity(
                    prepared.molecule,
                    options=XTBOptions(method='GFN2-xTB', confId=conformer_id, **common),
                    ionizedUnpaired=row.get('ionized_uhf'),
                    attachedUnpaired=row.get('attached_uhf'),
                    recordId=f"{row['id']}:{conformer_id}", returnErrors=True))
            item['analysis'] = primary
            item['analyses'] = {'gfn2_300K': primary, 'gfn1_300K': gfn1,
                                'gfn2_1000K': hot, 'gfn2_reactivity_300K': reactivity}
            def _geometry_hashes(report):
                return tuple(sorted(i.get('geometry_sha256') for i in report.get('inputs', [])
                                    if i.get('geometry_sha256')))
            prepared_hashes = tuple(sorted(c.get('geometry_sha256')
                                           for c in prepared.report.get('conformers', [])
                                           if c.get('geometry_sha256')))
            same_geometry = bool(prepared_hashes and _geometry_hashes(primary) == prepared_hashes and
                                 _geometry_hashes(primary) == _geometry_hashes(gfn1) == _geometry_hashes(hot))
            item['protocol_checks'] = {
                'same_optimized_geometry': same_geometry,
                'comparison_method': gfn1.get('calculation_provenance', {}).get('method'),
                'diagnostic_temperature_K': 1000.0,
                'reactivity_conformers': len(reactivity),
            }
            if not same_geometry:
                raise ValueError('GFN1/GFN2/temperature analyses did not retain identical optimized geometries')
            statuses = [primary.get('status'), gfn1.get('status'), hot.get('status')]
            reactivity_ok = all(r.get('error') is None and r.get('quality', {}).get('status') != 'error'
                                for r in reactivity)
            item['status'] = 'complete' if all(s == 'complete' for s in statuses) and reactivity_ok else (
                'failed' if primary.get('status') == 'failed' else 'partial')
        except Exception as exc:
            item['status'] = 'failed'
            item['error'] = {'type': type(exc).__name__, 'message': str(exc)}
        records.append(item)
    output = {'protocol_id': manifest['protocol_id'],
              'protocol_sha256': FrozenProtocolSHA256(protocol),
              'feature_method': 'GFN2-xTB', 'electronic_temperature_K': 300.0,
              'diagnostic_electronic_temperature_K': 1000.0,
              'n_records': len(records), 'records': records}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True, allow_nan=False) + '\n', encoding='utf-8')
    print(json.dumps({'n_records': len(records),
                      'n_complete': sum(r.get('status') == 'complete' for r in records),
                      'n_failed': sum(r.get('status') not in ('complete','partial') for r in records)}, indent=2))
    return 0

if __name__ == '__main__': raise SystemExit(main())
