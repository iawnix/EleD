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

from .charges import _charge_model_for_method
from .errors import XTBCacheError, XTBEnvironmentError, XTBNativeNotFoundError
from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .result import XTBResults
from .runtime import _engine_identity, _parameter_context, _parameter_file
from .xtb import _resolve_options, _run_xtb


def _prepare(mol: Any, options: XTBOptions) -> tuple[XTBOptions, dict]:
    xyz = _xyz_from_mol(mol, options.confId)
    charge, unpaired = _effective_charge_uhf(mol, options.charge, options.numUnpairedElectrons)
    metadata = _input_metadata(mol, options.confId, charge, unpaired, xyz)
    resolved = _resolve_options(
        replace(options, charge=charge, numUnpairedElectrons=unpaired, confId=metadata['conformer_id'])
    )
    return resolved, metadata


def _key(options: XTBOptions, metadata: dict) -> str:
    path = options.xtbLibrary
    engine = _engine_identity('native', path)
    if engine['sha256'] is None:
        raise ValueError('caching needs an absolute xTB library path so its binary can be fingerprinted')
    values = options.ToDict()
    values.pop('xtbLibrary', None)
    from .persistence import _implementation_identity
    payload = {
        'schema': 6,
        'input': metadata,
        'options': values,
        'implementation': _implementation_identity(),
        'engine': engine,
        'charge_model': _charge_model_for_method(options.method),
        'parameters': _parameter_context(),
        'parameter_file': _parameter_file(options.method, path),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True, allow_nan=False).encode()).hexdigest()


def XTBCacheKey(mol: Any, **options: Any) -> str:
    """Identify indexed input, settings, resolved engine version and binary hash."""
    config = XTBOptions(**options)
    resolved, metadata = _prepare(mol, config)
    return _key(resolved, metadata)


class XTBCache:
    """In-memory cache; concurrent identical requests share a single calculation."""

    def __init__(self):
        self._values: dict[str, XTBResults] = {}
        self._pending: dict[str, Future] = {}
        self._lock = RLock()

    def Get(self, key: str) -> XTBResults | None:
        with self._lock:
            return deepcopy(self._values.get(key))

    def Set(self, key: str, value: XTBResults) -> None:
        with self._lock:
            self._values[key] = deepcopy(value)

    def GetOrCompute(self, key: str, calculate: Callable[[], XTBResults]) -> XTBResults:
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
            self.Set(key, value)
            pending.set_result(value)
            return deepcopy(value)
        except BaseException as exc:
            pending.set_exception(exc)
            raise
        finally:
            with self._lock:
                self._pending.pop(key, None)

    def Clear(self) -> None:
        with self._lock:
            if self._pending:
                raise RuntimeError('cannot clear cache while calculations are running')
            self._values.clear()

    def __len__(self) -> int:
        with self._lock:
            return len(self._values)


def RunXTBBatch(
    molecules: Iterable[Any],
    *,
    cache: Any = None,
    maxWorkers: int = 1,
    returnExceptions: bool = False,
    options: XTBOptions | None = None,
    timeout: float | None = None,
    isolate: bool = False,
    **overrides: Any,
) -> list[XTBResults | Exception]:
    """Calculate in input order, optionally retaining per-molecule errors."""
    from .workers import isolated_map, molecule_bytes, single_point, validate_execution
    isolated = validate_execution(maxWorkers, timeout, isolate)
    config = replace(options or XTBOptions(), **overrides)
    if isolated:
        from contextlib import closing
        from rdkit import Chem
        config = _resolve_options(config)
        values, tasks, pending = [], [], []
        for mol in molecules:
            index = len(values)
            values.append(None)
            try:
                snapshot = Chem.Mol(mol) if mol is not None else None
                resolved, metadata = _prepare(snapshot, config)
                key = _key(resolved, metadata) if cache is not None else None
                cached = cache.Get(key) if cache is not None else None
                if cached is not None:
                    values[index] = cached
                else:
                    tasks.append((molecule_bytes(snapshot), resolved.ToDict()))
                    pending.append((index, key))
            except (XTBCacheError, sqlite3.Error, XTBEnvironmentError, XTBNativeNotFoundError):
                raise
            except (ValueError, TypeError) as exc:
                if not returnExceptions:
                    raise
                values[index] = exc
        with closing(isolated_map(single_point, tasks, maxWorkers=maxWorkers, timeout=timeout)) as completed:
            for task_index, value in completed:
                index, key = pending[task_index]
                if isinstance(value, Exception):
                    if isinstance(value, (XTBEnvironmentError, XTBNativeNotFoundError)) or not returnExceptions:
                        raise value
                elif cache is not None:
                    cache.Set(key, value)
                values[index] = value
        return values

    def calculate(mol):
        try:
            if cache is None:
                return _run_xtb(mol, **config.ToDict())
            resolved, metadata = _prepare(mol, config)
            key = _key(resolved, metadata)
            return cache.GetOrCompute(key, lambda: _run_xtb(mol, **resolved.ToDict()))
        except (XTBCacheError, sqlite3.Error, XTBEnvironmentError, XTBNativeNotFoundError):
            raise
        except Exception as exc:
            if returnExceptions:
                return exc
            raise

    return [calculate(mol) for mol in molecules]
