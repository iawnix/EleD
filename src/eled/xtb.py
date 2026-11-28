"""Internal xTB execution entry points.

The public RDKit-facing calculation API is :class:`XTBCalculator`. This module
keeps the engine adapter private so the package can evolve its result and
backend layers without carrying a legacy function API.
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

from .backends.native import _NativeXTB, _run_xtb_native
from .options import XTBOptions
from .result import XTBResults
from .runtime import _find_xtb_library, _require_extended_libxtb


def _run_xtb(mol: Any, **options: Any) -> XTBResults:
    """Run one single point (gas-phase or implicit-solvent) for ``XTBCalculator``.

    ``options`` are ``XTBOptions`` fields; the native backend validates them
    once (through ``XTBOptions``), resolves charge/spin and loads the library.
    """
    return _run_xtb_native(mol, **options)


def _resolve_options(options: XTBOptions) -> XTBOptions:
    library = _NativeXTB(_find_xtb_library(options.xtbLibrary))
    _require_extended_libxtb(library.library)
    return replace(options, xtbLibrary=library.library)
