"""Ordered batch execution with provenance-aware, in-memory result caching."""

from __future__ import annotations
from concurrent.futures import Future
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
import sqlite3
from threading import RLock
from typing import Any, Callable, Iterable

from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .result import XTBResult
from .runtime import _engine_identity, _parameter_context, _parameter_file
from .xtb import _resolve_options, run_xtb
from .errors import XTBCacheError, XTBEnvironmentError, XTBNativeNotFoundError


def _prepare(mol: Any, options: XTBOptions) -> tuple[XTBOptions, dict]:
    xyz = _xyz_from_mol(mol, options.conf_id)
    charge, uhf = _effective_charge_uhf(mol, options.charge, options.uhf)
    metadata = _input_metadata(mol, options.conf_id, charge, uhf, xyz)
    resolved = _resolve_options(replace(options, charge=charge, uhf=uhf,
                                       conf_id=metadata['conformer_id']))
    return resolved, metadata


def _key(options: XTBOptions, metadata: dict) -> str:
    path = options.xtb_library
    engine = _engine_identity(options.backend, path)
    if engine['sha256'] is None:
        raise ValueError('caching needs an absolute xTB library path so its binary can be fingerprinted')
    values = options.to_dict()
    values.pop('xtb_library', None)
    from .persistence import _implementation_identity
    payload = {'schema': 5, 'input': metadata, 'options': values,
               'implementation': _implementation_identity(),
               'engine': engine,
               'parameters': _parameter_context(),
               'parameter_file': _parameter_file(options.method, path)}
    return hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()


def xtb_cache_key(mol: Any, **options: Any) -> str:
    """Identify indexed input, settings, resolved engine version and binary hash.

    Requires the requested backend to be installed. Bare library sonames are
    usable for calculations, but caching requires a fingerprintable file path.
    """
    config = XTBOptions(**options)
    resolved, metadata = _prepare(mol, config)
    return _key(resolved, metadata)


class XTBCache:
    """In-memory cache; concurrent identical requests share a single calculation.

    Only successful results are stored. Copies protect the cached metadata and
    raw JSON from caller mutation. Clear the cache after changing external
    parameter/dependency installations during a running Python session.
    """

    def __init__(self):
        self._values: dict[str, XTBResult] = {}
        self._pending: dict[str, Future] = {}
        self._lock = RLock()

    def get(self, key: str) -> XTBResult | None:
        with self._lock:
            return deepcopy(self._values.get(key))

    def set(self, key: str, value: XTBResult) -> None:
        with self._lock:
            self._values[key] = deepcopy(value)

    def get_or_compute(self, key: str, calculate: Callable[[], XTBResult]) -> XTBResult:
        with self._lock:
            if key in self._values:
                return deepcopy(self._values[key])
            pending = self._pending.get(key)
            owner = pending is None
            if owner:
                pending = self._pending[key] = Future()
        if not owner:
            return deepcopy(pending.result())
        try:
            value = calculate()
            self.set(key, value)
            pending.set_result(value)
            return deepcopy(value)
        except BaseException as exc:
            pending.set_exception(exc)
            raise
        finally:
            with self._lock:
                self._pending.pop(key, None)

    def clear(self) -> None:
        with self._lock:
            if self._pending:
                raise RuntimeError('cannot clear cache while calculations are running')
            self._values.clear()

    def __len__(self) -> int:
        with self._lock:
            return len(self._values)


def run_xtb_batch(
    molecules: Iterable[Any], *, cache: XTBCache | None = None,
    max_workers: int = 1, return_exceptions: bool = False,
    options: XTBOptions | None = None, timeout: float | None = None,
    isolate: bool = False, **overrides: Any,
) -> list[XTBResult | Exception]:
    """Calculate in input order, optionally retaining per-molecule errors.

    ``options`` can hold serializable defaults; keyword overrides take priority.
    max_workers>1, timeout or isolate=True selects disposable spawned Python
    workers calling libxtb. Timeout includes worker startup, not queue time.
    The default remains in-process. Use a __main__ guard for spawned workers.
    """
    from .workers import validate_execution, isolated_map, single_point, molecule_bytes
    isolated = validate_execution(max_workers, timeout, isolate)
    config = replace(options or XTBOptions(), **overrides)
    if isolated:
        from contextlib import closing
        from rdkit import Chem
        config = _resolve_options(config)  # installation errors must not become failed molecular rows
        values, tasks, pending = [], [], []
        for mol in molecules:
            index = len(values); values.append(None)
            try:
                snapshot = Chem.Mol(mol) if mol is not None else None
                resolved, metadata = _prepare(snapshot, config)
                key = _key(resolved, metadata) if cache is not None else None
                cached = cache.get(key) if cache is not None else None
                if cached is not None:
                    values[index] = cached
                else:
                    tasks.append((molecule_bytes(snapshot), resolved.to_dict())); pending.append((index, key))
            except (XTBCacheError, sqlite3.Error, XTBEnvironmentError, XTBNativeNotFoundError): raise
            except (ValueError, TypeError) as exc:
                if not return_exceptions: raise
                values[index] = exc
        with closing(isolated_map(single_point, tasks, max_workers=max_workers, timeout=timeout)) as completed:
            for task_index, value in completed:
                index, key = pending[task_index]
                if isinstance(value, Exception):
                    if isinstance(value, (XTBEnvironmentError, XTBNativeNotFoundError)) or not return_exceptions:
                        raise value
                elif cache is not None: cache.set(key, value)
                values[index] = value
        return values

    def calculate(mol):
        try:
            if cache is None:
                return run_xtb(mol, **config.to_dict())
            resolved, metadata = _prepare(mol, config)
            key = _key(resolved, metadata)
            return cache.get_or_compute(key, lambda: run_xtb(mol, **resolved.to_dict()))
        except (XTBCacheError, sqlite3.Error, XTBEnvironmentError, XTBNativeNotFoundError):
            raise
        except Exception as exc:
            if return_exceptions:
                return exc
            raise

    return [calculate(mol) for mol in molecules]
