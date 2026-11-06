"""Actual spawned processes: timeout, crash isolation, cleanup and native parity."""
from contextlib import closing
import multiprocessing
import os
import time

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import (run_xtb_batch, run_descriptor_batch, XTBOptions, XTBCache,
                       XTBTimeoutError, XTBWorkerError, ReactivityJob, run_reactivity_batch)
from rdkit_eth.workers import isolated_map, validate_execution
from rdkit_eth.runtime import _extended_libxtb_build, _find_xtb_library
from rdkit_eth.errors import XTBNativeNotFoundError


def task(value):
    if value == 'crash': os._exit(23)
    if value == 'slow': time.sleep(10)
    if value == 'error': raise ValueError('explicit test failure')
    return value


def molecule(smiles='O'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=17) == 0
    return mol


def delayed_reactivity(payload):
    from rdkit_eth import run_reactivity
    binary, kwargs = payload
    mol = Chem.Mol(binary)
    if mol.GetNumAtoms() > 3: time.sleep(10)
    return run_reactivity(mol, **kwargs)


@pytest.fixture
def extended_runtime():
    try:
        library = _find_xtb_library()
    except XTBNativeNotFoundError:
        pytest.skip('extended libxtb unavailable')
    if not _extended_libxtb_build(library):
        pytest.skip('reviewed extended libxtb required')
    return library


@pytest.mark.parametrize('args', [(0, None, False), (True, None, False), (1, 0, False),
                                 (1, float('nan'), False), (1, True, False), (1, None, 1)])
def test_invalid_execution_options(args):
    with pytest.raises(ValueError): validate_execution(*args)


def test_crashed_worker_does_not_kill_other_jobs():
    values = dict(isolated_map(task, ['crash', 3, 'error', 7], max_workers=2, timeout=15))
    assert isinstance(values[0], XTBWorkerError)
    assert '23' in str(values[0])
    assert values[1] == 3 and values[3] == 7
    assert isinstance(values[2], ValueError)
    assert not multiprocessing.active_children()


def test_timeout_and_early_close_clean_up_children():
    values = dict(isolated_map(task, ['slow'], max_workers=1, timeout=1.))
    assert isinstance(values[0], XTBTimeoutError)
    with closing(isolated_map(task, [3, 'slow', 'slow'], max_workers=2, timeout=None)) as results:
        assert next(results) == (0, 3)
    assert not multiprocessing.active_children()


@pytest.mark.integration
def test_native_parallel_parity_cache_order_and_reports(extended_runtime):
    mols = [molecule('O'), molecule('CO')]
    serial = run_xtb_batch(mols)
    cache = XTBCache()
    parallel = run_xtb_batch(mols, max_workers=2, timeout=30, cache=cache)
    for a, b in zip(serial, parallel):
        assert a.total_energy_hartree == pytest.approx(b.total_energy_hartree, abs=1e-10)
        assert a.metadata['geometry_sha256'] == b.metadata['geometry_sha256']
    assert len(cache) == 2
    # A cache hit needs no worker, even with a budget too small for startup.
    hits = run_xtb_batch(mols, timeout=1e-6, cache=cache)
    assert hits[0].total_energy_hartree == parallel[0].total_energy_hartree
    reports = run_descriptor_batch(mols, max_workers=2, timeout=30,
                                  options=XTBOptions(keep_overlap=False))
    assert all(r['error'] is None for r in reports)
    failed = run_descriptor_batch(mols, max_workers=2, timeout=1e-6,
                                 options=XTBOptions(keep_overlap=False))
    assert all(r['error']['type'] == 'XTBTimeoutError' for r in failed)


@pytest.mark.integration
def test_reactivity_timeout_checkpoint_retry_and_resume(tmp_path, extended_runtime):
    jobs = [ReactivityJob('water', molecule())]
    path = tmp_path/'checkpoint.sqlite'
    failed = run_reactivity_batch(jobs, checkpoint=path, timeout=1e-6)
    assert failed[0]['error']['type'] == 'XTBTimeoutError'
    resumed = run_reactivity_batch(jobs, checkpoint=path, max_workers=2, timeout=30)
    assert resumed[0]['execution']['source'] == 'checkpoint'
    retried = run_reactivity_batch(jobs, checkpoint=path, max_workers=2, timeout=30, retry_failed=True)
    assert retried[0]['error'] is None and retried[0]['execution']['attempt'] == 2
    serial = run_reactivity_batch(jobs)
    assert retried[0]['values'] == pytest.approx(serial[0]['values'])


@pytest.mark.integration
def test_completed_rows_commit_before_a_slow_neighbor_finishes(tmp_path, monkeypatch, extended_runtime):
    from concurrent.futures import ThreadPoolExecutor
    from rdkit_eth import ReactivityCheckpoint
    import rdkit_eth.workers as workers
    monkeypatch.setattr(workers, 'reactivity', delayed_reactivity)
    checkpoint = ReactivityCheckpoint(tmp_path/'completion.sqlite')
    # Slow input comes first; completion commits must not wait for input ordering.
    jobs = [ReactivityJob('slow', molecule('CCO')), ReactivityJob('water', molecule())]
    with ThreadPoolExecutor(max_workers=1) as executor:
        future = executor.submit(run_reactivity_batch, jobs, checkpoint=checkpoint,
                                 max_workers=2, timeout=4.)
        deadline = time.monotonic()+8.
        while not checkpoint.history('water') and not future.done() and time.monotonic()<deadline:
            time.sleep(.02)
        assert checkpoint.history('water')[0]['report']['error'] is None
        assert not future.done()
        rows = future.result(timeout=8.)
    assert [r['record_id'] for r in rows] == ['slow', 'water']
    assert rows[0]['error']['type'] == 'XTBTimeoutError'
