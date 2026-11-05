#!/usr/bin/env python3
"""Check preservation of baseline successes on fresh graph-selected QM8 cases.

All 40 matching non-development identities, all 45 previous preset/mode/seed
combinations, plus 40 nonmatching controls under ETKDGv4/DG and five seeds.
Inputs are frozen before generation. Source coordinates never seed generation.
"""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from benchmark_qm8 import read_source, SHA256
from evaluate_ring_angle import PRESETS, MODES, SEEDS, assess, digest, invoke, target_centers
from run_diagnosis import attach_and_check
from nonregression_gate import compare_pair, verdict

ROOT = Path(__file__).resolve().parents[2]


def select_cohort(data):
    candidates, excluded, source_count = read_source(data)
    selected = json.load(gzip.open(ROOT/'docs/validation/0.9.0/qm8-manifest.json.gz', 'rt'))['selected']
    known_identities = {r['identity'] for r in selected}
    known_scaffolds = {r['scaffold'] for r in selected}
    affected, controls = [], []
    for row, _ in candidates:
        if row['identity'] in known_identities:
            continue
        match = bool(target_centers(row['smiles']))
        out = {k: row[k] for k in ('id', 'smiles', 'identity', 'scaffold')}
        out.update(group='affected' if match else 'control', scaffold_seen=row['scaffold'] in known_scaffolds)
        if match:
            affected.append(out)
        elif len(controls) < 40:
            controls.append(out)
    return {'source_count': source_count, 'eligible_unique': len(candidates),
            'excluded_count': len(excluded), 'rows': affected + controls,
            'selection': 'all matching non-development identities; first 40 nonmatching in existing stable hash order',
            'known_identity_count': len(known_identities),
            'affected_new_scaffold_molecules': sum(not r['scaffold_seen'] for r in affected)}


def measure(smiles, events):
    base = assess(smiles, events, optimize=False)
    initial_valid = bool(base['coordinate_check'] and base['coordinate_check']['ok'])
    report = {'probe': base, 'initial_valid': initial_valid, 'post_valid': False,
              'initial_energy_kcal_mol': None, 'post_energy_kcal_mol': None,
              'optimization_seconds': None, 'accepted_attempt': None}
    if not initial_valid:
        return report
    report['accepted_attempt'] = next(e['attempt'] for e in events if e.get('stage') == 'ACCEPTED')
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for i, xyz in enumerate(base['positions_angstrom']):
        conf.SetAtomPosition(i, xyz)
    mol.AddConformer(conf)
    params = AllChem.MMFFGetMoleculeProperties(mol)
    if params is None:
        report['post_status'] = 'missing_mmff_parameters'
        return report
    field = AllChem.MMFFGetMoleculeForceField(mol, params)
    report['initial_energy_kcal_mol'] = float(field.CalcEnergy())
    start = time.perf_counter()
    status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    report['optimization_seconds'] = time.perf_counter()-start
    report['post_status'] = status
    report['post_coordinate_check'] = attach_and_check(mol, mol.GetConformer().GetPositions())
    post_field = AllChem.MMFFGetMoleculeForceField(mol, AllChem.MMFFGetMoleculeProperties(mol))
    report['post_energy_kcal_mol'] = float(post_field.CalcEnergy())
    if not np.isfinite([report['initial_energy_kcal_mol'], report['post_energy_kcal_mol']]).all():
        raise ValueError('nonfinite MMFF energy')
    report['post_valid'] = status == 0 and report['post_coordinate_check']['ok']
    report['post_positions_angstrom'] = mol.GetConformer().GetPositions().tolist()
    coords = mol.GetConformer().GetPositions()
    ratios = []
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        length = np.linalg.norm(coords[i]-coords[j])
        denominator = sum(Chem.GetPeriodicTable().GetRcovalent(mol.GetAtomWithIdx(k).GetAtomicNum()) for k in (i, j))
        ratios.append(float(length/denominator))
    report['post_bond_covalent_radius_ratios'] = {'min': min(ratios), 'max': max(ratios)}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, required=True)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--planarity-only', type=Path, required=True)
    parser.add_argument('--combined', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != '2026.09.1':
        raise ValueError('requires RDKit 2026.09.1')
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise ValueError('output must be empty')
    cohort = select_cohort(args.data)
    (args.output/'cohort.json').write_text(json.dumps(cohort, indent=2)+'\n')
    binaries = {'baseline': args.baseline, 'planarity_only': args.planarity_only, 'combined': args.combined}
    builds = {}
    for name, binary in binaries.items():
        builds[name] = json.loads(binary.with_name('build-manifest.json').read_text())
        if digest(binary) != builds[name]['binary_sha256']:
            raise ValueError('binary changed: '+name)
    manifest = {'protocol': 'ring-nonregression-v1', 'builds': builds,
        'source_archive_sha256': SHA256, 'cohort_sha256': digest(args.output/'cohort.json'),
        'presets': PRESETS, 'modes': MODES, 'seeds': SEEDS, 'max_embedding_attempts': 100,
        'controls': 'ETKDGv4/DG only, all five seeds', 'post_mmff_max_iterations': 500,
        'energy_tolerance_kcal_mol': .001,
        'gates': 'zero lost initial successes; zero lost post-MMFF successes; no >.001 kcal/mol post-MMFF energy increase on paired valid results',
        'timing_policy': 'rotated execution order; subprocess/logging wall time reported, not claimed to be algorithm-only time',
        'sources': {p.name: digest(p) for p in [Path(__file__), Path(__file__).with_name('nonregression_gate.py'),
            Path(__file__).with_name('evaluate_ring_angle.py'), Path(__file__).with_name('run_diagnosis.py'),
            Path(__file__).with_name('audit_constraints.py'), Path(__file__).parent.parent/'benchmark_qm8.py']},
        'source_coordinates_used_for_generation': False}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    comparisons = defaultdict(list)
    counts = defaultdict(Counter)
    affected_successes = 0
    failures = []
    timings = defaultdict(list)
    trial = 0
    with gzip.open(args.output/'records.jsonl.gz', 'wt') as handle:
        for row in cohort['rows']:
            presets = PRESETS if row['group'] == 'affected' else ('ETKDGv4',)
            modes = MODES if row['group'] == 'affected' else ('dg',)
            for preset in presets:
                for mode in modes:
                    for seed in SEEDS:
                        labels = list(binaries)
                        order = labels[trial % 3:] + labels[:trial % 3]
                        reports, seconds, stderr = {}, {}, {}
                        events_by_label = {}
                        for label in order:
                            events, seconds[label], stderr[label] = invoke(binaries[label], row['smiles'], preset, mode, seed, 'UFF')
                            reports[label] = measure(row['smiles'], events)
                            events_by_label[label] = events
                            key = row['group']+'/'+label
                            counts[key]['runs'] += 1
                            counts[key]['initial_valid'] += reports[label]['initial_valid']
                            counts[key]['post_valid'] += reports[label]['post_valid']
                            timings[key].append(seconds[label])
                        original = reports['baseline']
                        if row['group'] == 'affected':
                            affected_successes += original['initial_valid']
                        record = {'id': row['id'], 'group': row['group'], 'smiles': row['smiles'],
                            'preset': preset, 'mode': mode, 'seed': seed, 'execution_order': order,
                            'reports': reports, 'seconds': seconds, 'stderr': stderr, 'comparisons': {}}
                        for label in ('planarity_only', 'combined'):
                            comparison = compare_pair(original, reports[label])
                            if row['group'] == 'control' and events_by_label[label] != events_by_label['baseline']:
                                comparison['failures'].append('changed_non_target_events')
                            comparisons[label].append(comparison)
                            record['comparisons'][label] = comparison
                            if comparison['failures']:
                                failures.append({'id': row['id'], 'group': row['group'], 'preset': preset,
                                    'mode': mode, 'seed': seed, 'variant': label, **comparison})
                                record.setdefault('regression_events', {})[label] = events_by_label[label]
                                record.setdefault('regression_events', {})['baseline'] = events_by_label['baseline']
                        handle.write(json.dumps(record, allow_nan=False)+'\n')
                        handle.flush()
                        trial += 1
            print(row['id'], row['group'], 'complete', trial, 'paired trials;', len(failures), 'flagged comparisons', flush=True)
    summary = {'counts': counts, 'paired_trials': trial, 'failures': failures,
        'verdicts': {label: verdict(values, affected_successes) for label, values in comparisons.items()},
        'diagnostic_wall_seconds': {key: {'total': sum(values), 'median': float(np.median(values)),
                                        'p95': float(np.percentile(values, 95))} for key, values in timings.items()}}
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary['verdicts'], indent=2), flush=True)


if __name__ == '__main__':
    main()
