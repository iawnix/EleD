"""Ordered batch execution with provenance-aware, in-memory result caching."""

from __future__ import annotations
from concurrent.futures import Future
from copy import deepcopy
from dataclasses import replace
import hashlib
import json
from threading import RLock
from typing import Any, Callable, Iterable

from .input import _effective_charge_uhf, _input_metadata, _xyz_from_mol
from .options import XTBOptions
from .result import XTBResult
from .runtime import _engine_identity, _parameter_context, _parameter_file
from .xtb import _resolve_options, run_xtb


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
    payload = {'schema': 4, 'input': metadata, 'options': values,
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
    options: XTBOptions | None = None, **overrides: Any,
) -> list[XTBResult | Exception]:
    """Calculate sequentially in input order, optionally retaining per-molecule errors.

    ``options`` can hold serializable defaults; keyword overrides take priority.
    Only max_workers=1 is supported. Native calls in this package are serialized;
    the former CLI subprocess workers have been removed.
    """
    if isinstance(max_workers, bool) or not isinstance(max_workers, int) or max_workers < 1:
        raise ValueError('max_workers must be a positive integer')
    config = replace(options or XTBOptions(), **overrides)
    if max_workers > 1:
        raise ValueError('parallel batch execution is unavailable; use max_workers=1')

    def calculate(mol):
        try:
            if cache is None:
                return run_xtb(mol, **config.to_dict())
            resolved, metadata = _prepare(mol, config)
            key = _key(resolved, metadata)
            return cache.get_or_compute(key, lambda: run_xtb(mol, **resolved.to_dict()))
        except Exception as exc:
            if return_exceptions:
                return exc
            raise

    return [calculate(mol) for mol in molecules]
