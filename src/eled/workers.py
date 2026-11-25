"""Spawned, disposable libxtb workers; no xTB executable or shared native handles."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor, as_completed
import math
import multiprocessing
from numbers import Real
from pathlib import Path
import pickle
import tempfile
import threading
import time
import traceback

from .errors import XTBEnvironmentError, XTBTimeoutError, XTBWorkerError


def validate_execution(maxWorkers, timeout, isolate):
    if isinstance(maxWorkers, bool) or not isinstance(maxWorkers, int) or maxWorkers < 1:
        raise ValueError('maxWorkers must be a positive integer')
    if not isinstance(isolate, bool):
        raise ValueError('isolate must be Boolean')
    if timeout is not None and (
        isinstance(timeout, bool) or not isinstance(timeout, Real) or not math.isfinite(timeout) or timeout <= 0
    ):
        raise ValueError('timeout must be a finite positive number of seconds or None')
    return isolate or maxWorkers > 1 or timeout is not None


def _entry(function, payload, output):
    try:
        value = function(payload)
    except Exception as exc:
        value = exc
        if hasattr(value, 'add_note'):
            value.add_note(traceback.format_exc())
    try:
        data = pickle.dumps(value, protocol=pickle.HIGHEST_PROTOCOL)
    except Exception as exc:
        data = pickle.dumps(XTBEnvironmentError('worker result serialization failed: ' + str(exc)))
    Path(output).write_bytes(data)


def _one(function, payload, timeout, cancelled):
    if cancelled.is_set():
        return XTBWorkerError('batch cancelled')
    context = multiprocessing.get_context('spawn')
    with tempfile.TemporaryDirectory(prefix='eled-worker-') as directory:
        output = Path(directory) / 'result.pickle'
        process = context.Process(target=_entry, args=(function, payload, str(output)))
        started = time.monotonic()
        try:
            process.start()
            while process.is_alive():
                if cancelled.is_set():
                    return XTBWorkerError('batch cancelled')
                if timeout is not None and time.monotonic() - started >= timeout:
                    return XTBTimeoutError(
                        f'isolated calculation exceeded {timeout:g} seconds (including worker startup)'
                    )
                process.join(0.05)
            if process.exitcode != 0 or not output.is_file():
                return XTBWorkerError(
                    f'isolated calculation exited without a result (exit code {process.exitcode})'
                )
            return pickle.loads(output.read_bytes())
        finally:
            if process.pid is not None:
                if process.is_alive():
                    process.terminate()
                    process.join(1.0)
                if process.is_alive():
                    process.kill()
                process.join()
                process.close()


def isolated_map(function, payloads, *, maxWorkers, timeout):
    """Yield (input index, result/exception) on completion; stop children on abort."""
    cancelled = threading.Event()
    executor = ThreadPoolExecutor(max_workers=maxWorkers)
    futures = {}
    try:
        for index, payload in enumerate(payloads):
            futures[executor.submit(_one, function, payload, timeout, cancelled)] = index
        for future in as_completed(futures):
            yield futures[future], future.result()
    finally:
        cancelled.set()
        for future in futures:
            future.cancel()
        executor.shutdown(wait=True, cancel_futures=True)


def single_point(payload):
    from rdkit import Chem
    from .calculator import XTBCalculator
    mol, options = payload
    return XTBCalculator(**options).Calculate(Chem.Mol(mol))


def reactivity(payload):
    from rdkit import Chem
    from .reactivity import RunReactivity
    mol, kwargs = payload
    return RunReactivity(Chem.Mol(mol), **kwargs)


def molecule_bytes(mol):
    from rdkit import Chem
    return mol.ToBinary(Chem.PropertyPickleOptions.AllProps | Chem.PropertyPickleOptions.CoordsAsDouble)
