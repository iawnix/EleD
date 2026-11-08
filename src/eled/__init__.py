"""An RDKit extension for electronic-structure descriptors backed by xTB.

The public package is centered on :class:`XTBCalculator` and the project-built
extended ``libxtb``.  The optional ``eled.descriptors`` module contains a
small adapter for RDKit's pre-existing YAeHMOP/EHT interface; it is not part of
the xTB descriptor protocol.
"""

from .xtb import (
    XTBError,
    XTBResult,
    XTBRunError,
    XTBNativeNotFoundError,
    calc_xtb_descriptors,
)
from .mapping import annotate_mol
from .batch import XTBCache, run_xtb_batch, xtb_cache_key
from .options import XTBOptions
from .calculator import XTBCalculator

from .electronic_descriptors import (
    calc_xtb_atom_descriptors, calc_xtb_bond_descriptors, calc_xtb_fragment_descriptors,
)
from .populations import calc_orbital_populations, calc_frontier_descriptors
from .density_descriptors import (
    calc_density_matrix, calc_shell_populations, calc_mayer_bond_descriptors,
    calc_frontier_bond_descriptors, calc_pi_descriptors, calc_density_descriptors,
)
from .errors import DescriptorUnavailableError

__all__ = [
    "XTBError",
    "XTBResult",
    "XTBRunError",
    "XTBNativeNotFoundError",
    "calc_xtb_descriptors",
    "annotate_mol",
    "XTBCache",
    "run_xtb_batch",
    "xtb_cache_key",
    "XTBOptions",
    "XTBCalculator",
    "calc_xtb_atom_descriptors",
    "calc_xtb_bond_descriptors",
    "calc_xtb_fragment_descriptors",
    "calc_orbital_populations",
    "calc_frontier_descriptors",
    "calc_density_matrix",
    "calc_shell_populations",
    "calc_mayer_bond_descriptors",
    "calc_frontier_bond_descriptors",
    "calc_pi_descriptors",
    "calc_density_descriptors",
    "DescriptorUnavailableError",
]

from .reports import (
    DESCRIPTOR_VERSION, descriptor_columns, assess_descriptor_quality,
    calc_descriptor_report, run_descriptor_batch, export_descriptor_reports,
    descriptor_sensitivity,
)
from .reactivity import calc_reactivity_descriptors, run_reactivity, ChargeStateCalculationError

__all__ += [
    'DESCRIPTOR_VERSION', 'descriptor_columns', 'assess_descriptor_quality',
    'calc_descriptor_report', 'run_descriptor_batch', 'export_descriptor_reports',
    'descriptor_sensitivity',
    'calc_reactivity_descriptors', 'run_reactivity', 'ChargeStateCalculationError',
]

from .persistence import XTBDiskCache, ReactivityCheckpoint
from .reactivity_batch import ReactivityJob, run_reactivity_batch
from .errors import XTBCacheError, XTBEnvironmentError
from .errors import XTBTimeoutError, XTBWorkerError

__all__ += ['XTBDiskCache', 'ReactivityCheckpoint', 'ReactivityJob',
            'run_reactivity_batch', 'XTBCacheError', 'XTBEnvironmentError',
            'XTBTimeoutError', 'XTBWorkerError']

from .conformers import conformer_descriptor_columns, run_conformer_analysis

__all__ += ['conformer_descriptor_columns', 'run_conformer_analysis']

from .geometry import ConformerGenerationResult, prepare_conformers

__all__ += ['ConformerGenerationResult', 'prepare_conformers']

from .optimization import XTBOptimizationResult, optimize_xtb

__all__ += ['XTBOptimizationResult', 'optimize_xtb']

from .eas import (normalize_eas_records, cluster_keys, assign_eas_splits,
                  assign_eas_folds, prepare_eas_dataset, load_frozen_protocol,
                  frozen_protocol_sha256)
from .descriptor_catalog import descriptor_catalog, descriptor_status

__all__ += ['normalize_eas_records', 'cluster_keys', 'assign_eas_splits', 'assign_eas_folds',
            'prepare_eas_dataset', 'load_frozen_protocol', 'frozen_protocol_sha256',
            'descriptor_catalog', 'descriptor_status']
