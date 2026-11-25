"""Component directions for the pinned xTB spherical AO convention."""
from __future__ import annotations

import numpy as np


def _native_p_directions(aoShellIndices, shellAngularMomenta) -> tuple[tuple[float, float, float], ...]:
    """Capture p_x/p_y/p_z directions at the native extraction boundary.

    The pinned xTB src/xtb/basisset.f90:set_p_function enumerates lao=2,3,4;
    src/intpack.f90 defines these as px,py,pz. Persist these vectors with AO
    rows so a later permutation or AO sign change can transform them too.
    Non-p rows contain zero vectors.
    """
    shells = np.asarray(aoShellIndices)
    directions = np.zeros((len(shells), 3))
    for sh, angularMomentum in enumerate(shellAngularMomenta):
        if angularMomentum != 1:
            continue
        rows = np.flatnonzero(shells == sh)
        if len(rows) != 3:
            raise ValueError('native p shell must have three AO components')
        directions[rows] = np.eye(3)
    return tuple(map(tuple, directions))
