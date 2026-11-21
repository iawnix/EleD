"""Shared validation for ESP sampling coordinates (Angstrom)."""
import numpy as np


def ESPPoints(values):
    """Freeze finite real [point, xyz] coordinates; None/empty means no sampling."""
    if values is None:
        return ()
    raw = np.asarray(values)
    if raw.size == 0 and raw.shape in ((0,), (0, 3)):
        return ()
    if raw.dtype.kind not in 'iuf' or raw.ndim != 2 or raw.shape[1] != 3:
        raise ValueError('esp_points_angstrom must be a real [npoints, 3] array')
    points = np.asarray(raw, dtype=float)
    if len(points) > 2147483647 or not np.isfinite(points).all():
        raise ValueError('esp_points_angstrom must be finite and fit the C API point count')
    return tuple(tuple(map(float, row)) for row in points)
