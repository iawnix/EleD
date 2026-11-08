"""Serializable calculation settings shared by batch and single-molecule APIs."""

from __future__ import annotations
from dataclasses import asdict, dataclass
import os
import math
from numbers import Integral, Real
from .input import _method_flag


def CalculationSettings(max_iterations=250, accuracy=1.0, electronic_temperature=300.0):
    """Validate values before calling C; record xTB's minimum temperature."""
    if (isinstance(max_iterations, bool) or not isinstance(max_iterations, Integral)
            or not 1 <= max_iterations <= 2147483647):
        raise ValueError('max_iterations must be an integer in [1, 2147483647]')
    for name, value in [('accuracy', accuracy), ('electronic_temperature', electronic_temperature)]:
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f'{name} must be a finite real number')
    if not 1e-4 <= accuracy <= 1e3:
        raise ValueError('accuracy must be in [0.0001, 1000]; avoid silent engine clamping')
    if electronic_temperature < 0:
        raise ValueError('electronic_temperature must be nonnegative (kelvin)')
    return {'max_iterations': int(max_iterations), 'accuracy': float(accuracy),
            'electronic_temperature': float(electronic_temperature),
            'effective_electronic_temperature': max(float(electronic_temperature), 1e-6)}


@dataclass(frozen=True)
class XTBOptions:
    conf_id: int = -1
    method: str = 'GFN2-xTB'
    charge: int | None = None
    uhf: int | None = None
    xtb_library: str | os.PathLike | None = None
    backend: str = 'native'
    keep_overlap: bool = False
    max_iterations: int = 250
    accuracy: float = 1.0
    electronic_temperature: float = 300.0

    def __post_init__(self):
        settings = CalculationSettings(self.max_iterations, self.accuracy, self.electronic_temperature)
        for key in ("max_iterations", "accuracy", "electronic_temperature"):
            object.__setattr__(self, key, settings[key])
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
        """Return JSON-compatible options for ``XTBCalculator`` and batches."""
        return asdict(self)
