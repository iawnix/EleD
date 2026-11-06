"""Reusable RDKit-facing calculators.

The calculator facade owns a validated :class:`XTBOptions` object and returns
the same immutable :class:`XTBResult` as ``run_xtb``.  Descriptor functions stay
separate so a calculation can be reused for several analyses without another
engine call.
"""
from __future__ import annotations

from dataclasses import fields, replace
from typing import Any

from .options import XTBOptions
from .result import XTBResult
from .xtb import run_xtb


_OPTION_NAMES = frozenset(field.name for field in fields(XTBOptions))


class XTBCalculator:
    """Reusable calculator for RDKit molecules and conformers.

    Parameters are validated when the calculator is constructed.  Per-call
    overrides are limited to fields of :class:`XTBOptions`; this keeps the
    RDKit-facing facade aligned with the functional ``run_xtb`` API.

    Examples
    --------
    >>> calculator = XTBCalculator(method="GFN2-xTB")
    >>> result = calculator.calculate(mol, keep_overlap=True)
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

    def to_options(self) -> XTBOptions:
        """Return the same immutable settings for batch/cache APIs."""

        return self._options

    def calculate(self, mol: Any, **overrides: Any) -> XTBResult:
        """Run one single point for an RDKit molecule.

        ``overrides`` accepts only ``XTBOptions`` fields, for example
        ``conf_id=2`` or ``keep_overlap=True``.  The input molecule is not
        modified.
        """

        unknown = sorted(set(overrides) - _OPTION_NAMES)
        if unknown:
            names = ', '.join(unknown)
            raise TypeError(f'unknown XTBCalculator override(s): {names}')
        selected = replace(self._options, **overrides) if overrides else self._options
        return run_xtb(mol, **selected.to_dict())

    __call__ = calculate

    def __repr__(self) -> str:
        settings = self._options.to_dict()
        settings.pop('xtb_library', None)
        return f'XTBCalculator({settings!r})'
