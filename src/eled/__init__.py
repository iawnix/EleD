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
    CalcXTBDescriptors,
)
from .mapping import AnnotateMol
from .batch import XTBCache, RunXTBBatch, XTBCacheKey
from .options import XTBOptions
from .calculator import XTBCalculator

from .electronic_descriptors import (
    CalcXTBAtomDescriptors, CalcXTBBondDescriptors, CalcXTBFragmentDescriptors,
)
from .populations import CalcOrbitalPopulations, CalcFrontierDescriptors
from .density_descriptors import (
    CalcDensityMatrix, CalcShellPopulations, CalcMayerBondDescriptors,
    CalcFrontierBondDescriptors, CalcPiDescriptors, CalcDensityDescriptors,
)
from .errors import DescriptorUnavailableError

__all__ = [
    "XTBError",
    "XTBResult",
    "XTBRunError",
    "XTBNativeNotFoundError",
    "CalcXTBDescriptors",
    "AnnotateMol",
    "XTBCache",
    "RunXTBBatch",
    "XTBCacheKey",
    "XTBOptions",
    "XTBCalculator",
    "CalcXTBAtomDescriptors",
    "CalcXTBBondDescriptors",
    "CalcXTBFragmentDescriptors",
    "CalcOrbitalPopulations",
    "CalcFrontierDescriptors",
    "CalcDensityMatrix",
    "CalcShellPopulations",
    "CalcMayerBondDescriptors",
    "CalcFrontierBondDescriptors",
    "CalcPiDescriptors",
    "CalcDensityDescriptors",
    "DescriptorUnavailableError",
]

from .reports import (
    DESCRIPTOR_VERSION, DescriptorColumns, AssessDescriptorQuality,
    CalcDescriptorReport, RunDescriptorBatch, ExportDescriptorReports,
    DescriptorSensitivity,
)
from .reactivity import CalcReactivityDescriptors, RunReactivity, ChargeStateCalculationError

__all__ += [
    'DESCRIPTOR_VERSION', 'DescriptorColumns', 'AssessDescriptorQuality',
    'CalcDescriptorReport', 'RunDescriptorBatch', 'ExportDescriptorReports',
    'DescriptorSensitivity',
    'CalcReactivityDescriptors', 'RunReactivity', 'ChargeStateCalculationError',
]

from .persistence import XTBDiskCache, ReactivityCheckpoint
from .reactivity_batch import ReactivityJob, RunReactivityBatch
from .errors import XTBCacheError, XTBEnvironmentError
from .errors import XTBTimeoutError, XTBWorkerError

__all__ += ['XTBDiskCache', 'ReactivityCheckpoint', 'ReactivityJob',
            'RunReactivityBatch', 'XTBCacheError', 'XTBEnvironmentError',
            'XTBTimeoutError', 'XTBWorkerError']

from .conformers import ConformerDescriptorColumns, RunConformerAnalysis

__all__ += ['ConformerDescriptorColumns', 'RunConformerAnalysis']

from .geometry import ConformerGenerationResult, PrepareConformers

__all__ += ['ConformerGenerationResult', 'PrepareConformers']

from .optimization import XTBOptimizationResult, OptimizeXTB

__all__ += ['XTBOptimizationResult', 'OptimizeXTB']

from .eas import (NormalizeEASRecords, ClusterKeys, AssignEASSplits,
                  AssignEASFolds, PrepareEASDataset, LoadFrozenProtocol,
                  FrozenProtocolSHA256)
from .descriptor_catalog import DescriptorCatalog, DescriptorStatus

from .visualization import (
    BuildFragmentVisualData, DrawDescriptorMap, PlotFragmentSummary,
    PlotFragmentInteractions, RenderFragmentAnalysis,
)

__all__ += ['NormalizeEASRecords', 'ClusterKeys', 'AssignEASSplits', 'AssignEASFolds',
            'PrepareEASDataset', 'LoadFrozenProtocol', 'FrozenProtocolSHA256',
            'DescriptorCatalog', 'DescriptorStatus']

__all__ += [
    'BuildFragmentVisualData', 'DrawDescriptorMap', 'PlotFragmentSummary',
    'PlotFragmentInteractions', 'RenderFragmentAnalysis',
]
