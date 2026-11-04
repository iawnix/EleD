"""Electronic-structure helpers for RDKit molecules.

The package keeps the RDKit YAeHMOP interface separate from the xTB backend.
"""

from .xtb import (
    XTBError,
    XTBResult,
    XTBRunError,
    XTBNativeNotFoundError,
    calc_xtb_descriptors,
    run_xtb,
    run_xtb_native,
)
from .descriptors import calc_eht_descriptors, calc_orbital_descriptors
from .mapping import annotate_mol
from .batch import XTBCache, run_xtb_batch, xtb_cache_key
from .options import XTBOptions

from .electronic_descriptors import (
    calc_xtb_atom_descriptors, calc_xtb_bond_descriptors, calc_xtb_fragment_descriptors,
)
from .populations import calc_orbital_populations, calc_frontier_descriptors
from .errors import DescriptorUnavailableError

__all__ = [
    "XTBError",
    "XTBResult",
    "XTBRunError",
    "XTBNativeNotFoundError",
    "calc_xtb_descriptors",
    "run_xtb",
    "run_xtb_native",
    "calc_eht_descriptors",
    "calc_orbital_descriptors",
    "annotate_mol",
    "XTBCache",
    "run_xtb_batch",
    "xtb_cache_key",
    "XTBOptions",
    "calc_xtb_atom_descriptors",
    "calc_xtb_bond_descriptors",
    "calc_xtb_fragment_descriptors",
    "calc_orbital_populations",
    "calc_frontier_descriptors",
    "DescriptorUnavailableError",
]

from .reports import (
    DESCRIPTOR_VERSION, descriptor_columns, assess_descriptor_quality,
    calc_descriptor_report, run_descriptor_batch, export_descriptor_reports,
)
from .reactivity import calc_reactivity_descriptors, run_reactivity, ChargeStateCalculationError

__all__ += [
    'DESCRIPTOR_VERSION', 'descriptor_columns', 'assess_descriptor_quality',
    'calc_descriptor_report', 'run_descriptor_batch', 'export_descriptor_reports',
    'calc_reactivity_descriptors', 'run_reactivity', 'ChargeStateCalculationError',
]
