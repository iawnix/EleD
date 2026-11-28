"""Internal xTB execution entry points.

The public RDKit-facing calculation API is :class:`XTBCalculator`. This module
keeps the engine adapter private so the package can evolve its result and
backend layers without carrying a legacy function API.
"""

from __future__ import annotations

from dataclasses import replace
import os
from typing import Any

from .backends.native import _NativeXTB, run_xtb_native
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .result import XTBResults
from .runtime import _find_xtb_library, _require_extended_libxtb


def _run_xtb(
    mol: Any,
    *,
    confId: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    numUnpairedElectrons: int | None = None,
    xtbLibrary: str | os.PathLike[str] | None = None,
    keepOverlap: bool = False,
    espPointsAngstrom=(),
    realSpacePointsAngstrom=(),
    realSpaceOccupationThreshold: float = 1e-7,
    realSpaceLEAWindowEV: float = 0.0,
    maxIterations: int = 250,
    accuracy: float = 1.0,
    electronicTemperature: float = 300.0,
    solvent: str | None = None,
    solvationModel: str = "alpb",
    solvationGrid: int | None = None,
) -> XTBResults:
    """Run one single point (gas-phase or implicit-solvent) for ``XTBCalculator``."""
    opts = XTBOptions(
        confId=confId,
        method=method,
        charge=charge,
        numUnpairedElectrons=numUnpairedElectrons,
        xtbLibrary=xtbLibrary,
        keepOverlap=keepOverlap,
        espPointsAngstrom=espPointsAngstrom,
        realSpacePointsAngstrom=realSpacePointsAngstrom,
        realSpaceOccupationThreshold=realSpaceOccupationThreshold,
        realSpaceLEAWindowEV=realSpaceLEAWindowEV,
        maxIterations=maxIterations,
        accuracy=accuracy,
        electronicTemperature=electronicTemperature,
        solvent=solvent,
        solvationModel=solvationModel,
        solvationGrid=solvationGrid,
    )
    _xyz_from_mol(mol, confId)
    charge, unpaired = _effective_charge_uhf(mol, charge, numUnpairedElectrons)
    opts = _resolve_options(replace(opts, charge=charge, numUnpairedElectrons=unpaired))
    return run_xtb_native(
        mol,
        confId=opts.confId,
        method=opts.method,
        charge=opts.charge,
        numUnpairedElectrons=opts.numUnpairedElectrons,
        xtbLibrary=opts.xtbLibrary,
        keepOverlap=opts.keepOverlap,
        espPointsAngstrom=opts.espPointsAngstrom,
        realSpacePointsAngstrom=opts.realSpacePointsAngstrom,
        realSpaceOccupationThreshold=opts.realSpaceOccupationThreshold,
        realSpaceLEAWindowEV=opts.realSpaceLEAWindowEV,
        maxIterations=opts.maxIterations,
        accuracy=opts.accuracy,
        electronicTemperature=opts.electronicTemperature,
        solvent=opts.solvent,
        solvationModel=opts.solvationModel,
        solvationGrid=opts.solvationGrid,
    )


def _resolve_options(options: XTBOptions) -> XTBOptions:
    library = _NativeXTB(_find_xtb_library(options.xtbLibrary))
    _require_extended_libxtb(library.library)
    return replace(options, xtbLibrary=library.library)
