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

    def __repr__(self) -> str:
        settings = self._options.ToDict()
        settings.pop('xtbLibrary', None)
        return f'XTBCalculator({settings!r})'


def _calculate(mol: Any, options: XTBOptions, interface: str = 'native') -> XTBResults:
    """Run one single point through the native ctypes or compiled ``rdXTBTools`` interface.

    This is the only place where the two interfaces are dispatched; both return
    the same normalized :class:`XTBResults` (see ``_as_result``).
    """
    if interface == 'native':
        return XTBCalculator(**options.ToDict()).Calculate(mol)
    if interface != 'compiled':
        raise ValueError("interface must be 'native' or 'compiled'")
    from . import rdXTBTools
    from .electronic_descriptors import _as_result
    from .errors import XTBEnvironmentError, XTBRunError
    ok, raw = rdXTBTools.RunMol(
        mol,
        confId=options.confId,
        method=options.method,
        charge=options.charge,
        numUnpairedElectrons=options.numUnpairedElectrons,
        keepOrbitalCoefficients=options.keepOverlap,
        xtbLibrary=options.xtbLibrary or '',
        maxIterations=options.maxIterations,
        accuracy=options.accuracy,
        electronicTemperature=options.electronicTemperature,
        keepOverlap=options.keepOverlap,
        espPointsAngstrom=options.espPointsAngstrom or None,
        solvent=options.solvent,
        solvationModel=options.solvationModel,
        solvationGrid=options.solvationGrid,
        realSpacePointsAngstrom=options.realSpacePointsAngstrom or None,
        realSpaceOccupationThreshold=options.realSpaceOccupationThreshold,
        realSpaceLEAWindowEV=options.realSpaceLEAWindowEV,
    )
    if not ok:
        error_type = XTBRunError if raw.failureKind == 'calculation' else XTBEnvironmentError
        raise error_type(raw.errorMessage)
    return _as_result(raw)
