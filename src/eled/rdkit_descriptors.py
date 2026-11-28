"""Single-step RDKit descriptor wrappers and ``rdMolDescriptors.Properties`` registration.

For multi-view analyses or datasets, prefer :class:`XTBCalculator` so one
calculation feeds molecular, atomic, bond, density and report views directly.
This module provides RDKit-style single-call convenience and property-registry
integration built on the same validated :class:`XTBCalculator` contract.
"""
from __future__ import annotations

from collections import OrderedDict
import math
import os
from threading import RLock
from typing import Any

from rdkit.Chem import Descriptors, rdMolDescriptors

from .calculator import XTBCalculator
from .electronic_descriptors import CalcXTBDescriptors
from .input import _effective_charge_uhf, _geometry_hash, _ordered_structure_hash, _xyz_from_mol
from .options import XTBOptions
from .result import XTBResults

# Keep Python functor instances alive for the lifetime of the process so C++
# rdMolDescriptors.Properties never holds dangling Python object references.
_REGISTERED_FUNCTORS: dict[str, Descriptors.PropertyFunctor] = {}
_REGISTRY_LOCK = RLock()


def _numeric_descriptor_names(frontierDepth: int, method: str) -> tuple[str, ...]:
    template = CalcXTBDescriptors(XTBResults(method), frontierDepth=frontierDepth)
    excluded = {'xtb_method'}
    names = [name for name in template if name not in excluded]
    if method in ('GFN1-xTB', 'GFN2-xTB'):
        for label in ('sum', 'min', 'max', 'mean', 'sd', 'abs'):
            names.append(f'xtb_cm5_q_{label}_e')
    return tuple(names)


def CalcMolXTBDescriptors(
    mol: Any,
    confId: int = -1,
    *,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    numUnpairedElectrons: int | None = None,
    xtbLibrary: str | os.PathLike | None = None,
    frontierDepth: int = 2,
    missingVal: Any = None,
    silent: bool = False,
    cache: Any = None,
) -> dict[str, Any]:
    """Calculate scalar xTB descriptors in one call, mirroring ``Descriptors.CalcMolDescriptors``.

    Unlike ``Descriptors.CalcMolDescriptors``, ``silent`` defaults to ``False``
    so invalid 3D conformers, implicit hydrogens or SCF failures raise
    explicitly unless suppressed by the caller.
    """
    if not isinstance(silent, bool):
        raise ValueError("silent must be a bool")
    options = XTBOptions(
        confId=confId,
        method=method,
        charge=charge,
        numUnpairedElectrons=numUnpairedElectrons,
        xtbLibrary=xtbLibrary,
    )
    try:
        calculator = XTBCalculator(**options.ToDict())
        if cache is not None:
            from .batch import XTBCacheKey
            key = XTBCacheKey(mol, **options.ToDict())
            result = cache.GetOrCompute(key, lambda: calculator.Calculate(mol))
        else:
            result = calculator.Calculate(mol)
        values = CalcXTBDescriptors(result, frontierDepth=frontierDepth)
    except Exception:
        if not silent:
            raise
        template = CalcXTBDescriptors(XTBResults(options.method), frontierDepth=frontierDepth)
        if options.method in ('GFN1-xTB', 'GFN2-xTB'):
            for label in ('sum', 'min', 'max', 'mean', 'sd', 'abs'):
                template[f'xtb_cm5_q_{label}_e'] = None
        return {name: missingVal for name in template}
    if missingVal is not None:
        values = {name: (missingVal if val is None else val) for name, val in values.items()}
    return values


class _SharedMolEvaluator:
    """Deduplicate consecutive property evaluations on the same molecule/conformer."""

    def __init__(
        self,
        options: XTBOptions,
        frontierDepth: int,
        cache: Any = None,
        maxEntries: int = 32,
    ) -> None:
        self._calculator = XTBCalculator(**options.ToDict())
        self._options = options
        self._frontier_depth = frontierDepth
        self._external_cache = cache
        self._max_entries = maxEntries
        self._recent: OrderedDict[tuple[Any, ...], dict[str, Any]] = OrderedDict()
        self._lock = RLock()

    def evaluate(self, mol: Any) -> dict[str, Any]:
        xyz = _xyz_from_mol(mol, self._options.confId)
        charge, unpaired = _effective_charge_uhf(
            mol, self._options.charge, self._options.numUnpairedElectrons
        )
        fingerprint = (
            _ordered_structure_hash(mol),
            _geometry_hash(xyz),
            charge,
            unpaired,
            self._options.method,
            self._options.xtbLibrary,
            self._frontier_depth,
        )
        with self._lock:
            if fingerprint in self._recent:
                self._recent.move_to_end(fingerprint)
                return self._recent[fingerprint]
        if self._external_cache is not None:
            from .batch import XTBCacheKey
            key = XTBCacheKey(mol, **self._options.ToDict())
            result = self._external_cache.GetOrCompute(key, lambda: self._calculator.Calculate(mol))
        else:
            result = self._calculator.Calculate(mol)
        values = CalcXTBDescriptors(result, frontierDepth=self._frontier_depth)
        with self._lock:
            self._recent[fingerprint] = values
            self._recent.move_to_end(fingerprint)
            while len(self._recent) > self._max_entries:
                self._recent.popitem(last=False)
        return values


class _XTBPropertyFunctor(Descriptors.PropertyFunctor):
    """RDKit ``PropertyFunctor`` extracting one scalar field from a shared evaluation."""

    def __init__(self, propName: str, fieldName: str, version: str, evaluator: _SharedMolEvaluator) -> None:
        super().__init__(propName, version)
        self._field_name = fieldName
        self._evaluator = evaluator

    def __call__(self, mol: Any) -> float:
        value = self._evaluator.evaluate(mol).get(self._field_name)
        return math.nan if value is None else float(value)


def RegisterRDKitProperties(
    *,
    method: str = "GFN2-xTB",
    confId: int = -1,
    xtbLibrary: str | os.PathLike | None = None,
    frontierDepth: int = 2,
    cache: Any = None,
    prefix: str = "",
) -> tuple[str, ...]:
    """Register numeric xTB scalar descriptors with ``rdMolDescriptors.Properties``.

    Returns the tuple of registered property names in order. All properties in
    a registration group share a bounded in-memory conformer evaluator so
    ``rdMolDescriptors.Properties(names).ComputeProperties(mol)`` executes at
    most one xTB single point per distinct molecule/geometry. Undefined scalar
    values return ``NaN`` as required by RDKit's numeric ``double`` property API.
    """
    from .reports import DESCRIPTOR_VERSION
    if not isinstance(prefix, str):
        raise ValueError("prefix must be a string")
    options = XTBOptions(confId=confId, method=method, xtbLibrary=xtbLibrary)
    fields = _numeric_descriptor_names(frontierDepth, options.method)
    evaluator = _SharedMolEvaluator(options, frontierDepth, cache=cache)
    registered: list[str] = []
    with _REGISTRY_LOCK:
        for field_name in fields:
            prop_name = f"{prefix}{field_name}"
            functor = _XTBPropertyFunctor(prop_name, field_name, DESCRIPTOR_VERSION, evaluator)
            _REGISTERED_FUNCTORS[prop_name] = functor
            rdMolDescriptors.Properties.RegisterProperty(functor)
            registered.append(prop_name)
    return tuple(registered)
