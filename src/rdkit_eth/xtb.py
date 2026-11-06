"""Internal xTB execution entry points.

The public RDKit-facing calculation API is :class:`XTBCalculator`. This module
keeps the engine adapter private so the package can evolve its result and
backend layers without carrying a legacy function API.
"""

from __future__ import annotations
import os
from dataclasses import replace
from typing import Any

from .errors import (
    XTBError,
    XTBRunError,
    XTBNativeNotFoundError,
)
from .result import XTBResult, _frontier_energies
from .electronic_descriptors import calc_xtb_descriptors
from .input import _xyz_from_mol, _effective_charge_uhf
from .options import XTBOptions
from .backends.native import run_xtb_native, _NativeXTB
from .runtime import _find_xtb_library, _require_extended_libxtb

def _run_xtb(
    mol: Any,
    *,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    uhf: int | None = None,
    xtb_library: str | os.PathLike[str] | None = None,
    backend: str = "native",
    keep_overlap: bool = False,
    max_iterations: int = 250,
    accuracy: float = 1.0,
    electronic_temperature: float = 300.0,
) -> XTBResult:
    """Run one gas-phase single point for ``XTBCalculator``.

    The selected library must contain the project basis-export, SCC-diagnostics
    and verified GFN0 extensions. No external engine or fallback is used.
    """

    opts = XTBOptions(conf_id=conf_id, method=method, charge=charge, uhf=uhf,
                      xtb_library=xtb_library, backend=backend, keep_overlap=keep_overlap,
                      max_iterations=max_iterations, accuracy=accuracy, electronic_temperature=electronic_temperature)
    _xyz_from_mol(mol, conf_id)
    charge, uhf = _effective_charge_uhf(mol, charge, uhf)
    opts = _resolve_options(replace(opts, charge=charge, uhf=uhf))
    return run_xtb_native(mol, conf_id=opts.conf_id, method=opts.method,
                          charge=opts.charge, uhf=opts.uhf, xtb_library=opts.xtb_library,
                          keep_overlap=opts.keep_overlap, max_iterations=opts.max_iterations,
                          accuracy=opts.accuracy, electronic_temperature=opts.electronic_temperature)


def _resolve_options(options: XTBOptions) -> XTBOptions:
    library = _NativeXTB(_find_xtb_library(options.xtb_library))
    _require_extended_libxtb(library.library)
    return replace(options, xtb_library=library.library)
