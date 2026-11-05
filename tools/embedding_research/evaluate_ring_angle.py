#!/usr/bin/env python3
"""Paired evaluation of the local ring-angle prototype; no reference initialization.

All frozen 1,024 SMILES: ETKDGv4/DG/UFF, seed 42, 100 attempts.
Topological target subset: 3 presets x 3 modes x 5 seeds, baseline/local/MMFF.
This is development-set evidence, not independent generalization validation.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from audit_constraints import audit
from run_diagnosis import attach_and_check, validate_replay

ROOT = Path(__file__).resolve().parents[2]
SEEDS = (17, 29, 42, 43, 101)
PRESETS = ('ETKDGv3', 'ETKDGv4', 'srETKDGv3')
MODES = ('dg', 'random', 'internal')


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def target_centers(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    return [a.GetIdx() for a in mol.GetAtoms() if
            a.GetAtomicNum() == 7 and a.GetFormalCharge() == 0 and
            a.GetDegree() == 3 and not a.GetIsAromatic() and
            a.GetHybridization() == Chem.HybridizationType.SP2 and
            a.IsInRingSize(3) and a.IsInRingSize(4)]


def invoke(binary, smiles, preset, mode, seed, ff):
    start = time.perf_counter()
    proc = subprocess.run([str(binary.resolve()), smiles, preset, mode, str(seed), '100', ff],
                          text=True, capture_output=True, check=True, timeout=120)
    events = [json.loads(line) for line in proc.stdout.splitlines()]
    if events[0]['kind'] != 'constraints' or events[-1]['kind'] != 'result':
        raise ValueError('unexpected probe output')
    return events, time.perf_counter()-start, proc.stderr


def assess(smiles, events, optimize=False):
    constraints, *snapshots, result = events
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    if result['atomic_numbers'] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
        raise ValueError('probe atom order differs')
    for snapshot in snapshots:
        replay = audit(snapshot['positions'], constraints)
        validate_replay(snapshot, replay)
    accepted = [s for s in snapshots if s['stage'] == 'ACCEPTED']
    if len(accepted) != int(result['conformer_id'] >= 0):
        raise ValueError('accepted snapshot count differs')
    report = {'result': result, 'coordinate_check': None, 'post_mmff': None,
              'bounds_sha256': hashlib.sha256(json.dumps(constraints['bounds']).encode()).hexdigest(),
              'snapshot_count': len(snapshots)}
    if not accepted:
        return report
    xyz = np.asarray(accepted[0]['positions'])[:, :3]
    report['positions_angstrom'] = xyz.tolist()
    report['coordinate_check'] = attach_and_check(mol, xyz)
    if optimize and report['coordinate_check']['ok']:
        conf = Chem.Conformer(mol.GetNumAtoms())
        conf.Set3D(True)
        for i, point in enumerate(xyz):
            conf.SetAtomPosition(i, point)
        mol.AddConformer(conf)
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            report['post_mmff'] = {'status': 'missing_parameters'}
        else:
            status = AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            report['post_mmff'] = {'status': status,
                'coordinate_check': attach_and_check(mol, mol.GetConformer().GetPositions())}
    return report


def tally(counter, report):
    counter['runs'] += 1
    counter['returned'] += report['result']['conformer_id'] >= 0
    counter['optimizer_exceptions'] += report['result']['optimizer_exception']
    counter['initial_valid'] += bool(report['coordinate_check'] and report['coordinate_check']['ok'])
    post = report['post_mmff']
    counter['post_mmff_valid_converged'] += bool(post and post['status'] == 0 and post['coordinate_check']['ok'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=Path, required=True)
    parser.add_argument('--control', type=Path, required=True)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if rdBase.rdkitVersion != '2026.09.1':
        raise ValueError('requires RDKit 2026.09.1')
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise ValueError('output directory must be empty')
    selected_path = ROOT/'docs/validation/0.9.0/qm8-manifest.json.gz'
    rows = json.load(gzip.open(selected_path, 'rt'))['selected']
    targets = [r for r in rows if target_centers(r['smiles'])]
    builds = {}
    for name in ('baseline', 'control', 'candidate'):
        binary = getattr(args, name)
        builds[name] = json.loads(binary.with_name('build-manifest.json').read_text())
        if builds[name]['binary_sha256'] != digest(binary):
            raise ValueError('binary identity mismatch: '+name)
    if builds['control']['patch_sha256'] is not None or not builds['candidate']['patch_sha256']:
        raise ValueError('expected an unmodified control and patched candidate')
    manifest = {'protocol': 'local-ring-angle-experiment-v1', 'builds': builds,
        'cohort_sha256': digest(selected_path), 'ids': [r['id'] for r in rows],
        'target_ids': [r['id'] for r in targets], 'target_selection': 'graph predicate in target_centers()',
        'screen': {'preset': 'ETKDGv4', 'mode': 'dg', 'seed': 42, 'force_field': 'UFF'},
        'sweep': {'presets': PRESETS, 'modes': MODES, 'seeds': SEEDS}, 'max_iterations': 100,
        'post_mmff': 'target sweep only; same MMFF94, 500 iterations',
        'script_sha256': digest(Path(__file__)),
        'dependencies': {n: digest(Path(__file__).with_name(n)) for n in ('audit_constraints.py', 'run_diagnosis.py')},
        'scope': 'development experiment; no new molecular holdout or time-matched benchmark',
        'source_coordinate_policy': 'source references read only after all generation finishes'}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    summary = {'screen': {'baseline': Counter(), 'candidate': Counter()},
               'control_parity': 0, 'unchanged_non_targets': 0, 'changed_bounds_ids': [],
               'screen_gains': [], 'screen_losses': [], 'sweep': {}, 'target_source_audits': {}}
    target_set = {r['id'] for r in targets}
    target_constraints = {}
    with gzip.open(args.output/'screen.jsonl.gz', 'wt') as handle:
        for n, row in enumerate(rows):
            smi = row['smiles']
            baseline, tb, eb = invoke(args.baseline, smi, 'ETKDGv4', 'dg', 42, 'UFF')
            control, tc, ec = invoke(args.control, smi, 'ETKDGv4', 'dg', 42, 'UFF')
            candidate, tp, ep = invoke(args.candidate, smi, 'ETKDGv4', 'dg', 42, 'UFF')
            if baseline != control:
                raise AssertionError('recompiled bounds changed baseline: '+row['id'])
            summary['control_parity'] += 1
            if row['id'] not in target_set:
                if candidate != control:
                    raise AssertionError('patch changed a non-target: '+row['id'])
                summary['unchanged_non_targets'] += 1
            if candidate[0]['bounds'] != baseline[0]['bounds']:
                summary['changed_bounds_ids'].append(row['id'])
            reports = {name: assess(smi, ev) for name, ev in (('baseline', baseline), ('candidate', candidate))}
            for name, report in reports.items():
                tally(summary['screen'][name], report)
            valid = [bool(reports[k]['coordinate_check'] and reports[k]['coordinate_check']['ok']) for k in ('baseline', 'candidate')]
            if valid == [False, True]: summary['screen_gains'].append(row['id'])
            if valid == [True, False]: summary['screen_losses'].append(row['id'])
            record = {'id': row['id'], 'smiles': smi, 'reports': reports,
                      'control_equal': True, 'candidate_equal': candidate == control,
                      'seconds': {'baseline': tb, 'control': tc, 'candidate': tp},
                      'stderr': {'baseline': eb, 'control': ec, 'candidate': ep}}
            if row['id'] in target_set:
                record['events'] = {'baseline': baseline, 'candidate': candidate}
            handle.write(json.dumps(record, allow_nan=False)+'\n')
            handle.flush()
            if (n+1) % 100 == 0: print('screen', n+1, '/', len(rows), flush=True)
    with gzip.open(args.output/'sweep.jsonl.gz', 'wt') as handle:
        for row in targets:
            for preset in PRESETS:
                for mode in MODES:
                    for seed in SEEDS:
                        for label, binary, ff in (('baseline', args.control, 'UFF'),
                                                  ('candidate', args.candidate, 'UFF'),
                                                  ('global_mmff', args.control, 'MMFF')):
                            events, elapsed, stderr = invoke(binary, row['smiles'], preset, mode, seed, ff)
                            report = assess(row['smiles'], events, optimize=True)
                            key = '/'.join((row['id'], preset, mode, label))
                            tally(summary['sweep'].setdefault(key, Counter()), report)
                            target_constraints[row['id'], preset, mode, label] = events[0]
                            record = {'id': row['id'], 'smiles': row['smiles'], 'preset': preset,
                                'mode': mode, 'seed': seed, 'variant': label, 'report': report,
                                'events': events, 'seconds': elapsed, 'stderr': stderr}
                            handle.write(json.dumps(record, allow_nan=False)+'\n')
                            handle.flush()
            print('sweep', row['id'], 'complete', flush=True)
    # Reference coordinates are diagnostic only, never available to invoke().
    ref_path = ROOT/'docs/research/embedding-constraints-v1/references.json.gz'
    refs = json.load(gzip.open(ref_path, 'rt'))
    for key, constraints in target_constraints.items():
        summary['target_source_audits']['/'.join(key)] = audit(refs[key[0]]['positions_angstrom'], constraints)
    summary['diagnostic_references_sha256'] = digest(ref_path)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps({k: v for k, v in summary.items() if k not in ('sweep', 'target_source_audits')}, indent=2), flush=True)


if __name__ == '__main__':
    main()
