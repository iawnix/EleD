#!/usr/bin/env python3
"""Matched QM8 source-vs-SMILES geometry study with conditional paired intervals.

Uses the frozen 0.9 pilot cohort, feature sets, model and split IDs. No model
tuning, source coordinates in embedding, xTB optimization or excited-state xTB.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import importlib.util
import json
from numbers import Integral
from pathlib import Path
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, Descriptors, rdFingerprintGenerator

from rdkit_eth import XTBOptions, XTBEnvironmentError, XTBRunError, run_xtb, calc_descriptor_report
from rdkit_eth.input import _geometry_hash, _xyz_from_mol
from rdkit_eth.persistence import _implementation_identity

# Resolve the sibling explicitly so the tool also works with isolated Python -I.
_spec = importlib.util.spec_from_file_location('qm8_pilot', Path(__file__).with_name('benchmark_qm8.py'))
pilot = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(pilot)

PROTOCOL = 'qm8-geometry-v1'
BOOTSTRAP_REPLICATES = 5000


def prepare_geometry(smiles, identifier):
    """Create coordinates only from SMILES; retain all preparation outcomes."""
    seed = int(hashlib.sha256((PROTOCOL + ':' + identifier).encode()).hexdigest()[:8], 16) % 2147483646 + 1
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('invalid frozen SMILES')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3(); params.randomSeed = seed; params.numThreads = 1
    params.pruneRmsThresh = .2
    info = {'smiles': smiles, 'seed': seed, 'embedding': 'ETKDGv3', 'requested_conformers': 8,
            'prune_rms_angstrom': .2, 'optimizer': 'MMFF94', 'max_iterations': 500,
            'selection': 'lowest MMFF energy among converged conformers; ID breaks ties',
            'source_coordinates_used': False}
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=8, params=params))
    info['embedded_ids'] = ids
    if not ids: return None, {**info, 'failure': 'embedding_failed'}
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        return None, {**info, 'failure': 'MMFF_parameters_unavailable'}
    outcomes = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=500, mmffVariant='MMFF94')
    if len(outcomes) != len(ids): raise RuntimeError('MMFF result count differs from embedded conformers')
    info['optimization_results'] = [{'conf_id': i, 'status': int(status), 'energy_kcal_mol': float(energy)}
                                    for i, (status, energy) in zip(ids, outcomes)]
    if any(not np.isfinite(r['energy_kcal_mol']) for r in info['optimization_results']):
        raise ValueError('nonfinite MMFF energy')
    converged = [(energy, i) for i, (status, energy) in zip(ids, outcomes) if status == 0]
    if not converged: return None, {**info, 'failure': 'no_converged_MMFF_conformer'}
    _, chosen = min(converged)
    conf = Chem.Conformer(mol.GetConformer(chosen))
    mol.RemoveAllConformers(); mol.AddConformer(conf, assignId=True)
    info.update(selected_embedded_id=chosen, conformer_id=mol.GetConformer().GetId(),
                atomic_numbers=[a.GetAtomicNum() for a in mol.GetAtoms()],
                positions_angstrom=mol.GetConformer().GetPositions().tolist(), failure=None)
    _xyz_from_mol(mol, -1)
    return mol, info


def paired_interval(reference, candidate, *, groups, seed, replicates=BOOTSTRAP_REPLICATES):
    """Paired candidate-minus-reference MAE with test-group percentile bootstrap.

    Each input is an ID -> (target, prediction) mapping. Resamples whole test
    groups; each drawn group retains its molecules. Models are NOT retrained.
    """
    if not reference or set(reference) != set(candidate) or set(reference) != set(groups):
        raise ValueError('paired predictions and groups must have identical nonempty IDs')
    if isinstance(replicates, bool) or not isinstance(replicates, Integral) or replicates < 1:
        raise ValueError('replicates must be a positive integer')
    ids = sorted(reference)
    differences = []
    for identifier in ids:
        ya, pa = reference[identifier]; yb, pb = candidate[identifier]
        if not all(np.isfinite(v) for v in (ya, pa, yb, pb)) or ya != yb:
            raise ValueError('nonfinite predictions or mismatched targets')
        differences.append(abs(pb-yb)-abs(pa-ya))
    labels = sorted(set(groups.values()))
    sums = np.zeros(len(labels)); counts = np.zeros(len(labels), dtype=int)
    index = {label: i for i, label in enumerate(labels)}
    for identifier, difference in zip(ids, differences):
        group = index[groups[identifier]]; sums[group] += difference; counts[group] += 1
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(labels), size=(replicates, len(labels)))
    boot = sums[draws].sum(axis=1) / counts[draws].sum(axis=1)
    interval = np.quantile(boot, [.025, .975]).tolist() if len(labels) > 1 else None
    return {'delta_mae_eV': float(np.mean(differences)), 'interval_95_eV': interval,
            'n_test': len(ids), 'n_groups': len(labels), 'bootstrap_replicates': replicates,
            'bootstrap_seed': seed, 'direction': 'candidate minus reference; negative favors candidate',
            'resampling': 'paired test groups with replacement; molecule-weighted MAE',
            'scope': 'conditional on these fitted models; no retraining or multiplicity adjustment',
            'warning': 'few test groups' if len(labels) < 20 else None}


def read_evidence(root):
    required = ('qm8-manifest.json.gz', 'qm8-features.jsonl.gz', 'qm8-splits.json.gz')
    checksums = json.loads((root/'SHA256SUMS.json').read_text())
    for name in required:
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != checksums[name]:
            raise ValueError('frozen evidence checksum mismatch: ' + name)
    def read(name):
        with gzip.open(root/name, 'rt') as handle: return json.load(handle)
    manifest = read(required[0]); splits = read(required[2])
    with gzip.open(root/required[1], 'rt') as handle: rows = [json.loads(line) for line in handle]
    if manifest['protocol'] != 'qm8-pilot-v1' or manifest['sha256'] != pilot.SHA256:
        raise ValueError('unsupported source protocol')
    if (manifest['options'] != XTBOptions(keep_overlap=True).to_dict()
            or manifest['basic_features'] != list(pilot.BASIC) or manifest['orbital_features'] != list(pilot.ORBITAL)
            or manifest['rdkit_features'] != list(pilot.TWO_D) + ['Morgan radius2, 1024bits, no chirality']):
        raise ValueError('feature or calculation protocol differs from frozen pilot')
    selected = {row['id']: row for row in manifest['selected']}
    if len(selected) != len(manifest['selected']) or set(selected) != {row['id'] for row in rows}:
        raise ValueError('geometry study requires the complete frozen pilot cohort')
    for row in rows:
        if any(row[k] != selected[row['id']][k] for k in selected[row['id']]):
            raise ValueError('feature metadata differs from manifest')
    pilot.validate_splits(rows, splits)
    return manifest, rows, splits, {name: checksums[name] for name in required}


def compare_predictions(rows, predictions):
    by_id = {r['id']: r for r in rows}
    indexed = {}
    for geometry, items in predictions.items():
        for r in items:
            key = geometry, r['split'], r['seed'], r['features']
            group = indexed.setdefault(key, {})
            if r['id'] in group: raise ValueError('duplicate prediction ID')
            group[r['id']] = (r['target_eV'], r['prediction_eV'])
    comparisons = []
    for split in ('random', 'scaffold'):
        for seed in pilot.SEEDS:
            contrasts = [(g, 'rdkit', g, 'rdkit_basic') for g in predictions]
            contrasts += [(g, 'rdkit_basic', g, 'rdkit_basic_orbital') for g in predictions]
            contrasts += [('source', f, 'generated', f) for f in ('rdkit_basic', 'rdkit_basic_orbital')]
            for ga, fa, gb, fb in contrasts:
                a = indexed[(ga, split, seed, fa)]; b = indexed[(gb, split, seed, fb)]
                grouping = 'identity' if split == 'random' else 'scaffold'
                groups = {i: by_id[i][grouping] for i in a}
                comparisons.append({'split': split, 'seed': seed, 'reference': [ga, fa], 'candidate': [gb, fb],
                                    **paired_interval(a, b, groups=groups, seed=9000+seed)})
            a = indexed[('source', split, seed, 'rdkit')]; b = indexed[('generated', split, seed, 'rdkit')]
            if set(a) != set(b) or any(abs(a[i][1]-b[i][1]) > 1e-12 for i in a):
                raise ValueError('2D baseline changed between geometry protocols')
    return comparisons


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference', type=Path, default=Path(__file__).resolve().parents[1]/'docs/validation/0.9.0')
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()): raise ValueError('output directory must be empty')
    frozen, source, splits, checksums = read_evidence(args.reference)
    manifest = {'protocol': PROTOCOL, 'reference_evidence_sha256': checksums,
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'pilot_script_sha256': hashlib.sha256(Path(pilot.__file__).read_bytes()).hexdigest(),
        'implementation': _implementation_identity(), 'rdkit': rdBase.rdkitVersion,
        'selected_ids': [r['id'] for r in source], 'target': frozen['target'],
        'options': frozen['options'], 'embedding': 'ETKDGv3; 8 requested conformers; RMS prune .2 angstrom',
        'geometry_selection': 'lowest converged MMFF94 energy, 500 iterations, no source coordinates',
        'split_policy': 'restrict original frozen train/test IDs to common successful cohort; no resampling',
        'failure_policy': 'retain all failures, no replacement; reference refit on the same common cohort',
        'bootstrap_replicates': BOOTSTRAP_REPLICATES, 'model_tuning': False}
    pilot.write_json(args.output/'manifest.json', manifest)
    rows, failures, preparations = [], [], []
    fpgen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024, includeChirality=False)
    started = time.perf_counter()
    with gzip.open(args.output/'generated-features.jsonl.gz', 'wt') as output:
        for n, original in enumerate(source):
            identifier = original['id']; t = time.perf_counter()
            mol, info = prepare_geometry(original['smiles'], identifier)
            info.update(id=identifier, elapsed_seconds=time.perf_counter()-t); preparations.append(info)
            if mol is None:
                failures.append({'id': identifier, 'stage': 'geometry', 'reason': info['failure']}); continue
            # Confirm the canonical reconstruction has identical 2D features.
            heavy = Chem.RemoveHs(mol)
            two_d = [float(getattr(Descriptors, name)(heavy)) for name in pilot.TWO_D] + fpgen.GetFingerprintAsNumPy(heavy).astype(float).tolist()
            if not np.allclose(two_d, original['rdkit'], atol=1e-10, rtol=1e-12):
                raise ValueError('2D features changed after canonical reconstruction: ' + identifier)
            t = time.perf_counter()
            try: result = run_xtb(mol, **frozen['options'])
            except XTBEnvironmentError: raise
            except XTBRunError as exc:
                failures.append({'id': identifier, 'stage': 'xtb', 'reason': str(exc)}); continue
            elapsed = time.perf_counter()-t
            if (result.metadata['engine']['sha256'] != original['provenance']['engine']['sha256']
                    or (result.metadata['parameter_file'] or {}).get('sha256') !=
                       (original['provenance']['parameter_file'] or {}).get('sha256')):
                raise ValueError('engine/parameter identity differs from reference evidence')
            t = time.perf_counter(); report = calc_descriptor_report(mol, result)
            analysis_seconds = time.perf_counter()-t
            basic = [report['values'][k] for k in pilot.BASIC]
            orbital = [report['values'][k] for k in pilot.ORBITAL]
            if report['quality']['status'] == 'error' or any(v is None or not np.isfinite(v) for v in basic+orbital):
                failures.append({'id': identifier, 'stage': 'features', 'reason': 'quality_error_or_missing_features',
                                 'quality': report['quality']}); continue
            row = {k: original[k] for k in ('id', 'smiles', 'identity', 'scaffold', 'target_eV', 'rdkit')}
            row.update(basic=basic, orbital=orbital, provenance=report['provenance'], quality=report['quality'],
                source_geometry_sha256=original['geometry_sha256'], geometry_sha256=_geometry_hash(_xyz_from_mol(mol, -1)),
                geometry_seconds=info['elapsed_seconds'], xtb_seconds=elapsed, analysis_seconds=analysis_seconds)
            output.write(json.dumps(row, allow_nan=False)+'\n'); output.flush(); rows.append(row)
            if (n+1) % 100 == 0: print(f'generated/calculated {n+1}/{len(source)}', flush=True)
    with gzip.open(args.output/'geometry-preparation.jsonl.gz', 'wt') as handle:
        for row in preparations: handle.write(json.dumps(row, allow_nan=False)+'\n')
    pilot.write_json(args.output/'failures.json', failures)
    common = {r['id'] for r in rows}
    matched_source = [r for r in source if r['id'] in common]
    fixed = [{**s, 'train': [i for i in s['train'] if i in common],
                    'test': [i for i in s['test'] if i in common]} for s in splits]
    pilot.validate_splits(rows, fixed)
    summary = {'n_selected': len(source), 'n_common': len(rows), 'n_failed': len(failures),
               'failure_counts': dict(Counter(r['reason'] for r in failures)),
               'featurization_wall_seconds': time.perf_counter()-started,
               'geometry_seconds_all_attempts': sum(p['elapsed_seconds'] for p in preparations),
               'xtb_seconds_successful': sum(r['xtb_seconds'] for r in rows),
               'analysis_seconds_successful': sum(r['analysis_seconds'] for r in rows)}
    predictions = {}
    for label, features in [('source', matched_source), ('generated', rows)]:
        folder = args.output/label; folder.mkdir()
        summary[label] = pilot.fit_ablation(features, folder, splits=fixed)
        predictions[label] = json.loads((folder/'predictions.json').read_text())
    pilot.write_json(args.output/'paired-comparisons.json', compare_predictions(rows, predictions))
    pilot.write_json(args.output/'summary.json', summary)


if __name__ == '__main__': main()
