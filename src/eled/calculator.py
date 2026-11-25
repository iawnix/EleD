"""Reusable RDKit-facing calculators.

The calculator facade owns a validated :class:`XTBOptions` object and returns
an immutable :class:`XTBResults`. Descriptor functions stay separate so a
calculation can be reused for several analyses without another engine call.
"""
from __future__ import annotations

from dataclasses import fields, replace
from typing import Any

from .options import XTBOptions
from .result import XTBResults
from .xtb import _run_xtb


_OPTION_NAMES = frozenset(field.name for field in fields(XTBOptions))


class XTBCalculator:
    """Reusable calculator for RDKit molecules and conformers.

    Parameters are validated when the calculator is constructed. Per-call
    overrides are limited to fields of :class:`XTBOptions` (``confId``,
    ``keepOverlap``, ``numUnpairedElectrons``, ``espPointsAngstrom``, etc.).

    Examples
    --------
    >>> calculator = XTBCalculator(method="GFN2-xTB")
    >>> result = calculator.Calculate(mol, keepOverlap=True)
    """

    def __init__(self, **options: Any) -> None:
        unknown = sorted(set(options) - _OPTION_NAMES)
        if unknown:
            names = ', '.join(unknown)
            raise TypeError(f'unknown XTBCalculator option(s): {names}')
        self._options = XTBOptions(**options)

    @property
    def options(self) -> XTBOptions:
        """Return the immutable calculation settings owned by this object."""
        return self._options

    def GetOptions(self) -> XTBOptions:
        """Return the immutable calculation settings owned by this object."""
        return self._options

    def Calculate(self, mol: Any, **overrides: Any) -> XTBResults:
        """Run one single point for an RDKit molecule.

        ``overrides`` accepts only ``XTBOptions`` fields, for example
        ``confId=2``, ``keepOverlap=True`` or ``espPointsAngstrom=points``.
        The input molecule is not modified.
        """
        unknown = sorted(set(overrides) - _OPTION_NAMES)
        if unknown:
            names = ', '.join(unknown)
            raise TypeError(f'unknown XTBCalculator override(s): {names}')
        selected = replace(self._options, **overrides) if overrides else self._options
        return _run_xtb(mol, **selected.ToDict())

    __call__ = Calculate

    def __repr__(self) -> str:
        settings = self._options.ToDict()
        settings.pop('xtbLibrary', None)
        return f'XTBCalculator({settings!r})'
