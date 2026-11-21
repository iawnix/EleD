"""Explicit xTB atomic-charge models.

``xtb_getCharges`` does not return CM5 charges.  For GFN1/2 it returns the
converged SCC Mulliken charge (nuclear charge minus the Mulliken population);
GFN0 uses the PEEQ/EEQ charge model.  CM5 is a geometry-dependent correction
to Mulliken charges and is implemented here with the same parameters as the
reviewed xTB source.
"""
from __future__ import annotations

import math
from typing import Sequence

CHARGE_MODEL_SCC = "xtb_scc_mulliken"
CHARGE_MODEL_EEQ = "gfn0_eeq"
CHARGE_MODEL_CM5 = "cm5"

# Mantina/Valero/Cramer/Truhlar covalent radii in Angstrom, Z=1..118.
_RADII_A = (
    0.32, 0.37, 1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62,
    1.60, 1.40, 1.24, 1.14, 1.09, 1.04, 1.00, 1.01, 2.00, 1.74,
    1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, 1.17, 1.22, 1.20,
    1.23, 1.20, 1.20, 1.18, 1.17, 1.16, 2.15, 1.90, 1.76, 1.64,
    1.56, 1.46, 1.38, 1.36, 1.34, 1.30, 1.36, 1.40, 1.42, 1.40,
    1.40, 1.37, 1.36, 1.36, 2.38, 2.06, 1.94, 1.84, 1.90, 1.88,
    1.86, 1.85, 1.83, 1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78,
    1.74, 1.64, 1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30, 1.32,
    1.44, 1.45, 1.50, 1.42, 1.48, 1.46, 2.42, 2.11, 2.01, 1.90,
    1.84, 1.83, 1.80, 1.80, 1.73, 1.68, 1.68, 1.68, 1.65, 1.67,
    1.73, 1.76, 1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28,
    1.21, 1.22, 1.36, 1.43, 1.62, 1.75, 1.65, 1.57,
)

# Atomic CM5 parameters a0, Z=1..118.  The xTB implementation uses the
# pair-specific values below for H/C/N/O; all other pairs use a0[i]-a0[j].
_A0 = (
    0.0056, -0.1543, 0.0000, 0.0333, -0.1030, -0.0446, -0.1072,
    -0.0802, -0.0629, -0.1088, 0.0184, 0.0000, -0.0726, -0.0790,
    -0.0756, -0.0565, -0.0444, -0.0767, 0.0130, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, -0.0512, -0.0557, -0.0533, -0.0399, -0.0313, -0.0541,
    0.0092, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, -0.0361, -0.0393, -0.0376, -0.0281,
    -0.0220, -0.0381, 0.0065, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, -0.0255, -0.0277, -0.0265, -0.0198, -0.0155,
    -0.0269, 0.0046, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
    0.0000, 0.0000, 0.0000, 0.0000, -0.0179, -0.0195, -0.0187,
    -0.0140, -0.0110, -0.0189,
)

_PAIR = {(1, 6): 0.0502, (1, 7): 0.1747, (1, 8): 0.1671,
         (6, 7): 0.0556, (6, 8): 0.0234, (7, 8): -0.0346}


def charge_model_for_method(method: str) -> str:
    """Return the exact model used by xTB's atomic-charge result getter."""
    if method == "GFN0-xTB":
        return CHARGE_MODEL_EEQ
    if method in ("GFN1-xTB", "GFN2-xTB"):
        return CHARGE_MODEL_SCC
    raise ValueError(f"no atomic charge model is defined for {method!r}")


def calc_cm5_charges(atomic_numbers: Sequence[int], coordinates_angstrom: Sequence[Sequence[float]],
                     mulliken_charges: Sequence[float]) -> tuple[float, ...]:
    """Apply the xTB CM5 correction to SCC Mulliken charges.

    Coordinates are Angstrom and charges are in elementary-charge units.  CM5
    is defined for a Mulliken base charge; callers must not pass GFN0 EEQ
    charges.  The pair correction is antisymmetric, so the molecular charge is
    conserved to floating-point precision.
    """
    numbers = tuple(int(z) for z in atomic_numbers)
    coords = tuple(tuple(float(v) for v in row) for row in coordinates_angstrom)
    q = tuple(float(v) for v in mulliken_charges)
    n = len(numbers)
    if len(coords) != n or len(q) != n:
        raise ValueError("CM5 inputs must have matching atom counts")
    if any(z < 1 or z > 118 for z in numbers):
        raise ValueError("CM5 supports atomic numbers 1..118")
    if any(len(row) != 3 or any(not math.isfinite(v) for v in row) for row in coords):
        raise ValueError("CM5 coordinates must be finite 3-vectors")
    if any(not math.isfinite(v) for v in q):
        raise ValueError("Mulliken charges must be finite")
    correction = [0.0] * n
    for i in range(n):
        zi = numbers[i]
        for j in range(i):
            zj = numbers[j]
                dx = coords[i][0] - coords[j][0]
                dy = coords[i][1] - coords[j][1]
                dz = coords[i][2] - coords[j][2]
                distance = math.sqrt(dx * dx + dy * dy + dz * dz)
                if distance == 0.0:
                    raise ValueError("CM5 requires distinct atom coordinates")
                if zi == zj:
                    continue
                key = (min(zi, zj), max(zi, zj))
                if key in _PAIR:
                    base = _PAIR[key]
                    pij = base if zi < zj else -base
                else:
                    pij = _A0[zi - 1] - _A0[zj - 1]
            weight = math.exp(-2.4740 * (distance - _RADII_A[zi - 1] - _RADII_A[zj - 1]))
            correction[i] += weight * pij
            correction[j] -= weight * pij
    return tuple(qi + di for qi, di in zip(q, correction))

