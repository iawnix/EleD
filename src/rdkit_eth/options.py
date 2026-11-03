"""Serializable calculation settings shared by batch and single-molecule APIs."""

from __future__ import annotations
from dataclasses import asdict, dataclass
import os
from .input import _method_flag


@dataclass(frozen=True)
class XTBOptions:
    conf_id: int = -1
    method: str = 'GFN2-xTB'
    charge: int | None = None
    uhf: int | None = None
    xtb_library: str | os.PathLike | None = None
    backend: str = 'native'
    keep_overlap: bool = False

    def __post_init__(self):
        object.__setattr__(self, 'method', _method_flag(self.method)[0])
        object.__setattr__(self, 'backend', self.backend.strip().lower())
        if self.method == 'GFN-FF':
            raise ValueError('GFN-FF has no molecular orbitals and is not supported by this descriptor API')
        if self.backend != 'native':
            raise ValueError("backend must be 'native'; CLI, auto and tblite backends have been removed")
        if not isinstance(self.keep_overlap, bool):
            raise ValueError('keep_overlap must be a bool')
        if self.xtb_library is not None:
            object.__setattr__(self, 'xtb_library', os.fspath(self.xtb_library))

    def to_dict(self) -> dict:
        """Return JSON-compatible options, also suitable for ``run_xtb(mol, **opts)``."""
        return asdict(self)
