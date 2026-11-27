"""Serializable calculation settings shared by batch and single-molecule APIs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from numbers import Integral, Real
import os
from typing import Any, Mapping

from .input import _method_flag


_VALID_SOLVATION_MODELS = frozenset({'alpb', 'gbsa'})
_VALID_SOLVATION_GRIDS = frozenset({
    6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434,
    590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890,
    4334, 4802, 5294, 5810,
})


def _validate_solvation_settings(solvent=None, solvationModel='alpb', solvationGrid=None):
    """Validate optional implicit solvation settings without mutating gas-phase defaults."""
    if solvent is not None:
        if not isinstance(solvent, str) or not solvent.strip():
            raise ValueError('solvent must be a nonempty solvent name string or None')
        solvent = solvent.strip()
    if not isinstance(solvationModel, str) or not solvationModel.strip():
        raise ValueError("solvationModel must be 'alpb' or 'gbsa'")
    model = solvationModel.strip().lower()
    if model not in _VALID_SOLVATION_MODELS:
        raise ValueError("solvationModel must be 'alpb' or 'gbsa'")
    if solvationGrid is not None:
        if isinstance(solvationGrid, bool) or not isinstance(solvationGrid, Integral):
            raise ValueError('solvationGrid must be an integer Lebedev grid size or None')
        grid = int(solvationGrid)
        if grid not in _VALID_SOLVATION_GRIDS:
            raise ValueError('solvationGrid must be a valid xTB Lebedev grid size (e.g. 110, 230, 302, 434)')
    else:
        grid = None
    return solvent, model, grid


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
    realSpacePointsAngstrom: tuple[tuple[float, float, float], ...] = ()
    realSpaceOccupationThreshold: float = 1e-7
    realSpaceLEAWindowEV: float = 0.0
    solvent: str | None = None
    solvationModel: str = 'alpb'
    solvationGrid: int | None = None

    def __post_init__(self) -> None:
        from .esp_input import _normalize_esp_points
        if isinstance(self.confId, bool) or not isinstance(self.confId, Integral):
            raise ValueError('confId must be an integer')
        object.__setattr__(self, 'confId', int(self.confId))
        object.__setattr__(self, 'espPointsAngstrom', _normalize_esp_points(self.espPointsAngstrom))
        object.__setattr__(
            self,
            'realSpacePointsAngstrom',
            _normalize_esp_points(self.realSpacePointsAngstrom),
        )
        for label, attr in (
            ('realSpaceOccupationThreshold', 'realSpaceOccupationThreshold'),
            ('realSpaceLEAWindowEV', 'realSpaceLEAWindowEV'),
        ):
            val = getattr(self, attr)
            if isinstance(val, bool) or not isinstance(val, Real) or not math.isfinite(val) or val < 0:
                raise ValueError(f'{label} must be a finite nonnegative real number')
            object.__setattr__(self, attr, float(val))
        settings = _validate_calculation_settings(self.maxIterations, self.accuracy, self.electronicTemperature)
        object.__setattr__(self, 'maxIterations', settings['max_iterations'])
        object.__setattr__(self, 'accuracy', settings['accuracy'])
        object.__setattr__(self, 'electronicTemperature', settings['electronic_temperature'])
        solvent, solv_model, solv_grid = _validate_solvation_settings(
            self.solvent, self.solvationModel, self.solvationGrid
        )
        object.__setattr__(self, 'solvent', solvent)
        object.__setattr__(self, 'solvationModel', solv_model)
        object.__setattr__(self, 'solvationGrid', solv_grid)
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
        values['realSpacePointsAngstrom'] = [list(row) for row in self.realSpacePointsAngstrom]
        # Preserve existing cache keys and serialization when real-space sampling is unset.
        if (
            not self.realSpacePointsAngstrom
            and self.realSpaceOccupationThreshold == 1e-7
            and self.realSpaceLEAWindowEV == 0.0
        ):
            values.pop('realSpacePointsAngstrom', None)
            values.pop('realSpaceOccupationThreshold', None)
            values.pop('realSpaceLEAWindowEV', None)
        # Preserve gas-phase cache keys and serialization when solvation is unset.
        if self.solvent is None and self.solvationModel == 'alpb' and self.solvationGrid is None:
            values.pop('solvent', None)
            values.pop('solvationModel', None)
            values.pop('solvationGrid', None)
        return values

    @classmethod
    def FromDict(cls, values: Mapping[str, Any]) -> XTBOptions:
        """Reconstruct validated ``XTBOptions`` from :meth:`ToDict` output."""
        return cls(**dict(values))
