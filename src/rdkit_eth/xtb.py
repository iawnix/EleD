"""Public xTB entry points; backend implementations are kept separate."""

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
from .runtime import _find_xtb_library, _native_gfn0_build

def run_xtb(
    mol: Any,
    *,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    uhf: int | None = None,
    xtb_library: str | os.PathLike[str] | None = None,
    backend: str = "native",
    keep_overlap: bool = False,
) -> XTBResult:
    """Run a gas-phase libxtb single point, optionally retaining native basis data.

    keep_overlap requires the verified project basis extension. GFN0 also
    requires the D4 state fix. No external engine or fallback is used.
    """

    opts = XTBOptions(conf_id=conf_id, method=method, charge=charge, uhf=uhf,
                      xtb_library=xtb_library, backend=backend, keep_overlap=keep_overlap)
    _xyz_from_mol(mol, conf_id)
    charge, uhf = _effective_charge_uhf(mol, charge, uhf)
    opts = _resolve_options(replace(opts, charge=charge, uhf=uhf))
    return run_xtb_native(mol, conf_id=opts.conf_id, method=opts.method,
                          charge=opts.charge, uhf=opts.uhf, xtb_library=opts.xtb_library,
                          keep_overlap=opts.keep_overlap)


def _resolve_options(options: XTBOptions) -> XTBOptions:
    library = _NativeXTB(_find_xtb_library(options.xtb_library))
    if options.method == 'GFN0-xTB' and not _native_gfn0_build(library.library):
        raise ValueError('native GFN0 requires a verified patched xTB build; '
                         'run tools/build_patched_xtb.sh')
    return replace(options, xtb_library=library.library)
