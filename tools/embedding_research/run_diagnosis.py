#!/usr/bin/env python3
"""Frozen source-level diagnosis, with installed-library parity and source audit.

Source coordinates are only loaded after all generation subprocesses finish.
This is a diagnostic study of selected hard cases, not a general benchmark.
"""
import argparse
from collections import Counter
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdDistGeom
from eled.geometry import _check_geometry, _stereo_constraints

from audit_constraints import audit

ROOT = Path(__file__).resolve().parents[2]
SOURCE_SHA256 = 'e668b58f8dfba329bc20cd7e0eb3944e8a73a5d9db37160f686a8b31b20548a0'


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def attach_and_check(mol, positions):
    copy = Chem.Mol(mol)
    copy.RemoveAllConformers()
    conf = Chem.Conformer(copy.GetNumAtoms())
    conf.Set3D(True)
    if len(positions) != copy.GetNumAtoms():
        raise ValueError('coordinate count mismatch')
    for i, xyz in enumerate(positions):
        conf.SetAtomPosition(i, [float(v) for v in xyz[:3]])
    copy.AddConformer(conf)
    return _check_geometry(copy, 0, _stereo_constraints(copy))


def validate_replay(snapshot, result):
    for actual, replay in zip(snapshot['chiral_checks'], result['chiral'], strict=True):
        if not np.isclose(actual['volume'], replay['volume'], rtol=1e-12, atol=1e-12):
            raise AssertionError('signed volume differs from upstream')
        if actual['center_inside'] != replay['center_inside']:
            raise AssertionError('chiral containment differs from upstream')
    for actual, replay in zip(snapshot['tetrahedral_checks'], result['tetrahedral'], strict=True):
        if any(actual[k] != replay[k] for k in ('volume_ok', 'center_inside')):
            raise AssertionError('tetrahedral predicate differs from upstream')
    if snapshot['stage'] in ('FINAL_CHIRAL_BOUNDS', 'FINAL_CENTER_IN_VOLUME', 'CHECK_TETRAHEDRAL_CENTERS'):
        if snapshot['stage'] not in result['failed_predicates']:
            raise AssertionError('reported failure is absent from replay')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--probe', type=Path, required=True)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--force-field', choices=('UFF', 'MMFF'), default='UFF',
                        help='existing RDKit bounds parameterization, not a new method')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    if any(args.output.iterdir()):
        raise ValueError('output directory must be empty')
    if rdBase.rdkitVersion != '2026.09.1':
        raise ValueError('this diagnosis requires RDKit 2026.09.1')
    evidence = ROOT/'docs/validation'
    cohort_path = evidence/'0.10.0/preparation-manifest.json'
    prior_path = evidence/'0.10.0/preparation-summary.json'
    selected_path = evidence/'0.9.0/qm8-manifest.json.gz'
    ids = json.loads(cohort_path.read_text())['ids']
    prior = json.loads(prior_path.read_text())
    remaining = set(prior['remaining_ids'])
    recovered = set(prior['recovered_ids'])
    selected = {r['id']: r for r in json.load(gzip.open(selected_path, 'rt'))['selected']}
    rows = [selected[i] for i in ids]
    presets = ('ETKDGv3', 'ETKDGv4', 'srETKDGv3')
    modes = ('dg', 'internal')
    seed, budget = 42, 100
    build_manifest = args.probe.with_name('build-manifest.json')
    build_info = json.loads(build_manifest.read_text())
    if digest(args.probe) != build_info['binary_sha256']:
        raise ValueError('probe binary differs from build manifest')
    manifest = {'protocol': 'embedding-constraint-diagnosis-v1', 'ids': ids,
        'presets': presets, 'modes': modes, 'seed': seed, 'max_iterations': budget,
        'embed_force_field': args.force_field,
        'rdkit': rdBase.rdkitVersion, 'numpy': np.__version__, 'probe_build': build_info,
        'inputs': {str(p.relative_to(ROOT)): digest(p) for p in (cohort_path, prior_path, selected_path)},
        'tools': {p.name: digest(p) for p in Path(__file__).parent.glob('*.py')},
        'source_archive_sha256': SOURCE_SHA256,
        'scope': '81 selected cases; one conformer; no MMFF/xTB; max 3 snapshots per failure type; all counters retained',
        'reference_policy': 'source coordinates loaded only after generation; never supplied to embedding'}
    (args.output/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
    summary = {'groups': {}, 'parity_checks': 0, 'snapshot_replays': 0,
               'source_geometry_rejections': [], 'source_predicate_failures': {},
               'source_failure_ids': {}, 'source_sign_mismatch_ids': []}
    records = []
    for row in rows:
        identifier, smiles = row['id'], row['smiles']
        group = 'remaining_19' if identifier in remaining else 'recovered_22' if identifier in recovered else 'control_40'
        for preset in presets:
            for mode in modes:
                start = time.perf_counter()
                proc = subprocess.run([str(args.probe.resolve()), smiles, preset, mode, str(seed), str(budget), args.force_field],
                                      check=True, text=True, capture_output=True, timeout=120)
                events = [json.loads(line) for line in proc.stdout.splitlines()]
                elapsed = time.perf_counter()-start
                constraints, *snapshots, result = events
                assert constraints['kind'] == 'constraints' and result['kind'] == 'result'
                mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
                if result['atomic_numbers'] != [a.GetAtomicNum() for a in mol.GetAtoms()]:
                    raise AssertionError('probe atom order mismatch')
                p = getattr(rdDistGeom, preset)()
                p.randomSeed, p.maxIterations, p.numThreads, p.trackFailures = seed, budget, 1, True
                p.embedForceField = getattr(rdDistGeom, args.force_field)
                if mode == 'internal':
                    p.initialEmbeddingMode = rdDistGeom.INTERNAL_COORDINATE_EMBEDDING
                parity_start = time.perf_counter()
                optimizer_exception = False
                try:
                    cid = rdDistGeom.EmbedMolecule(mol, p)
                except RuntimeError as exc:
                    if 'bad direction in linearSearch' not in str(exc):
                        raise
                    optimizer_exception, cid = True, -1
                parity_seconds = time.perf_counter()-parity_start
                if optimizer_exception != result['optimizer_exception']:
                    raise AssertionError('probe/library exception differs')
                if cid != result['conformer_id'] or list(p.GetFailureCounts()) != result['failure_counts']:
                    raise AssertionError(f'probe/library outcome differs: {identifier}, {preset}, {mode}')
                accepted = [s for s in snapshots if s['stage'] == 'ACCEPTED']
                check = None
                if cid >= 0:
                    assert len(accepted) == 1
                    coords = np.asarray(accepted[0]['positions'])[:, :3]
                    if not np.allclose(mol.GetConformer(cid).GetPositions(), coords, rtol=0, atol=1e-10):
                        raise AssertionError('probe/library coordinates differ')
                    check = attach_and_check(mol, coords)
                summary['parity_checks'] += 1
                for snapshot in snapshots:
                    replay = audit(snapshot['positions'], constraints)
                    validate_replay(snapshot, replay)
                    snapshot['audit'] = replay
                    summary['snapshot_replays'] += 1
                label = preset+'/'+mode
                totals = summary['groups'].setdefault(group, {}).setdefault(label,
                    {'total': 0, 'returned': 0, 'stereo_valid': 0, 'optimizer_exceptions': 0,
                     'installed_embedding_seconds': 0., 'failure_counts': {}})
                totals['total'] += 1
                totals['returned'] += cid >= 0
                totals['stereo_valid'] += bool(check and check['ok'])
                totals['optimizer_exceptions'] += optimizer_exception
                totals['installed_embedding_seconds'] += parity_seconds
                names = {int(v): str(v) for v in rdDistGeom.EmbedFailureCauses.values.values()}
                counts = Counter(totals['failure_counts'])
                counts.update({names.get(i, f'UNKNOWN_SLOT_{i}'): n for i, n in enumerate(result['failure_counts']) if n})
                totals['failure_counts'] = dict(counts)
                records.append({'id': identifier, 'smiles': smiles, 'group': group,
                    'preset': preset, 'mode': mode, 'embed_force_field': args.force_field,
                    'probe_seconds': elapsed, 'parity_verified': True,
                    'stderr': proc.stderr, 'constraints': constraints, 'snapshots': snapshots,
                    'result': result, 'coordinate_check': check})
                with gzip.open(args.output/'generation-checkpoint.jsonl.gz', 'at') as checkpoint:
                    checkpoint.write(json.dumps(records[-1], allow_nan=False)+'\n')
        print(identifier, group, 'generation/parity complete', flush=True)

    # Strict separation: reference coordinates become available only here.
    if digest(args.archive) != SOURCE_SHA256:
        raise ValueError('QM8 archive differs from frozen input')
    sources = {}
    with tarfile.open(args.archive) as archive:
        for mol in Chem.ForwardSDMolSupplier(io.BytesIO(archive.extractfile('qm8.sdf').read()), removeHs=False):
            if mol is None:
                continue
            identifier = mol.GetProp('_Name').split()[1]
            if identifier not in selected or identifier not in ids:
                continue
            query = Chem.AddHs(Chem.MolFromSmiles(selected[identifier]['smiles']))
            if mol.GetNumAtoms() != query.GetNumAtoms():
                raise ValueError('source/query atom count differs')
            match = mol.GetSubstructMatch(query, useChirality=True)
            if len(match) != query.GetNumAtoms():
                raise ValueError('source stereochemical mapping failed: '+identifier)
            coords = mol.GetConformer().GetPositions()[list(match)]
            check = attach_and_check(query, coords)
            sources[identifier] = {'mapping_query_to_source': list(match),
                'positions_angstrom': coords.tolist(), 'coordinate_check': check}
            if not check['ok']:
                summary['source_geometry_rejections'].append(identifier)
    if set(sources) != set(ids):
        raise ValueError('source geometries missing')
    unique_failures = {}
    sign_failures = set()
    for record in records:
        identifier, group = record['id'], record['group']
        source = sources[identifier]
        result = audit(source['positions_angstrom'], record['constraints'])
        record['source_audit'] = result
        label = group+'/'+record['preset']+'/'+record['mode']
        counts = Counter(summary['source_predicate_failures'].get(label, {}))
        counts.update(result['failed_predicates'])
        summary['source_predicate_failures'][label] = dict(counts)
        for cause in result['failed_predicates']:
            unique_failures.setdefault(group+'/'+cause, set()).add(identifier)
        if any(not c['sign_ok'] for c in result['chiral']):
            sign_failures.add(identifier)
    summary['source_failure_ids'] = {k: sorted(v, key=int) for k, v in unique_failures.items()}
    summary['source_sign_mismatch_ids'] = sorted(sign_failures, key=int)
    with gzip.open(args.output/'records.jsonl.gz', 'wt') as handle:
        for record in records:
            handle.write(json.dumps(record, allow_nan=False)+'\n')
    with gzip.open(args.output/'references.json.gz', 'wt') as handle:
        json.dump(sources, handle, allow_nan=False)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2, allow_nan=False)+'\n')
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
