"""Serializable calculation settings shared by batch and single-molecule APIs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from numbers import Integral, Real
import os
from typing import Any, Mapping

from .input import _method_flag


def _validate_calculation_settings(maxIterations=250, accuracy=1.0, electronicTemperature=300.0) -> dict[str, Any]:
    """Validate values before calling C; record xTB's minimum temperature."""
    if (isinstance(maxIterations, bool) or not isinstance(maxIterations, Integral)
            or not 1 <= maxIterations <= 2147483647):
        raise ValueError('maxIterations (iteration limit) must be an integer in [1, 2147483647]')
    for label, value in [('accuracy', accuracy), ('electronic temperature', electronicTemperature)]:
        if isinstance(value, bool) or not isinstance(value, Real) or not math.isfinite(value):
            raise ValueError(f'{label} must be a finite real number')
    if not 1e-4 <= accuracy <= 1e3:
        raise ValueError('accuracy must be in [0.0001, 1000]; avoid silent engine clamping')
    if electronicTemperature < 0:
        raise ValueError('electronic temperature must be nonnegative (kelvin)')
    return {
        'max_iterations': int(maxIterations),
        'accuracy': float(accuracy),
        'electronic_temperature': float(electronicTemperature),
        'effective_electronic_temperature': max(float(electronicTemperature), 1e-6),
    }


@dataclass(frozen=True)
class XTBOptions:
    confId: int = -1
    method: str = 'GFN2-xTB'
    charge: int | None = None
    numUnpairedElectrons: int | None = None
    xtbLibrary: str | os.PathLike | None = None
    backend: str = 'native'
    keepOverlap: bool = False
    maxIterations: int = 250
    accuracy: float = 1.0
    electronicTemperature: float = 300.0
    espPointsAngstrom: tuple[tuple[float, float, float], ...] = ()

    def __post_init__(self) -> None:
        from .esp_input import _normalize_esp_points
        if isinstance(self.confId, bool) or not isinstance(self.confId, Integral):
            raise ValueError('confId must be an integer')
        object.__setattr__(self, 'confId', int(self.confId))
        object.__setattr__(self, 'espPointsAngstrom', _normalize_esp_points(self.espPointsAngstrom))
        settings = _validate_calculation_settings(self.maxIterations, self.accuracy, self.electronicTemperature)
        object.__setattr__(self, 'maxIterations', settings['max_iterations'])
        object.__setattr__(self, 'accuracy', settings['accuracy'])
        object.__setattr__(self, 'electronicTemperature', settings['electronic_temperature'])
        object.__setattr__(self, 'method', _method_flag(self.method)[0])
        object.__setattr__(self, 'backend', self.backend.strip().lower())
        if self.method == 'GFN-FF':
            raise ValueError('GFN-FF has no molecular orbitals and is not supported by this descriptor API')
        if self.backend != 'native':
            raise ValueError("backend must be 'native'; CLI, auto and tblite backends have been removed")
        if not isinstance(self.keepOverlap, bool):
            raise ValueError('keepOverlap must be a bool')
        if self.xtbLibrary is not None:
            object.__setattr__(self, 'xtbLibrary', os.fspath(self.xtbLibrary))

    def ToDict(self) -> dict[str, Any]:
        """Return JSON-compatible options for ``XTBCalculator`` and batches."""
        values = asdict(self)
        values['espPointsAngstrom'] = [list(row) for row in self.espPointsAngstrom]
        return values

    @classmethod
    def FromDict(cls, values: Mapping[str, Any]) -> XTBOptions:
        """Reconstruct validated ``XTBOptions`` from :meth:`ToDict` output."""
        return cls(**dict(values))
