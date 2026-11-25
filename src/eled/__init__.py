"""An RDKit extension for electronic-structure descriptors backed by xTB.

The public package is centered on :class:`XTBCalculator`, :class:`XTBResults`
and the project-built extended ``libxtb``, following RDKit's ``PascalCase``
function/method and ``camelCase`` parameter/property conventions. Dataset
splitting utilities for the frozen EAS study live in :mod:`eled.eas`.
"""

from .batch import RunXTBBatch, XTBCache, XTBCacheKey
from .calculator import XTBCalculator
from .conformers import ConformerDescriptorColumns, RunConformerAnalysis
from .density_descriptors import (
    CalcDensityDescriptors,
    CalcDensityMatrix,
    CalcFrontierBondDescriptors,
    CalcMCI,
    CalcMCIDescriptors,
    CalcMayerBondDescriptors,
    CalcPiDescriptors,
    CalcShellPopulations,
)
from .descriptor_catalog import DescriptorCatalog, DescriptorStatus
from .descriptors import CalcEHTDescriptors, CalcOrbitalDescriptors
from .electronic_descriptors import (
    CalcCM5Charges,
    CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors,
    CalcXTBFragmentDescriptors,
)
from .errors import (
    DescriptorUnavailableError,
    XTBCacheError,
    XTBEnvironmentError,
    XTBError,
    XTBNativeNotFoundError,
    XTBRunError,
    XTBTimeoutError,
    XTBWorkerError,
)
from .esp import CalcESP, ESPChargeFit, FitESPCharges, FitRESPCharges, GenerateESPGrid
from .geometry import ConformerGenerationResult, PrepareConformers
from .local_geometry import (
    AnalyzeLocalGeometry,
    CalcBuriedVolume,
    CalcHOMA,
    CalcPlaneDescriptors,
    CalcRingPuckering,
    CalcSterimol,
    RunLocalGeometryBatch,
)
from .mapping import AnnotateMol
from .optimization import OptimizeXTB, XTBOptimizationResult
from .options import XTBOptions
from .persistence import ReactivityCheckpoint, XTBDiskCache
from .populations import CalcFrontierDescriptors, CalcOrbitalPopulations
from .rdkit_descriptors import CalcMolXTBDescriptors, RegisterRDKitProperties
from .reactivity import CalcReactivityDescriptors, ChargeStateCalculationError, RunReactivity
from .reactivity_batch import ReactivityJob, RunReactivityBatch
from .reports import (
    DESCRIPTOR_VERSION,
    AssessDescriptorQuality,
    CalcDescriptorReport,
    DescriptorColumns,
    DescriptorSensitivity,
    ExportDescriptorReports,
    RunDescriptorBatch,
)
from .selections import NormalizeSelections, SelectionSet
from .visualization import (
    BuildFragmentVisualData,
    DrawDescriptorMap,
    PlotFragmentInteractions,
    PlotFragmentSummary,
    RenderFragmentAnalysis,
)
from .xtb import CalcXTBDescriptors, XTBResult, XTBResults

__all__ = [
    "XTBError",
    "XTBResults",
    "XTBResult",
    "XTBRunError",
    "XTBNativeNotFoundError",
    "XTBCacheError",
    "XTBEnvironmentError",
    "XTBTimeoutError",
    "XTBWorkerError",
    "DescriptorUnavailableError",
    "ChargeStateCalculationError",
    "XTBOptions",
    "XTBCalculator",
    "XTBCache",
    "XTBDiskCache",
    "XTBCacheKey",
    "RunXTBBatch",
    "CalcMolXTBDescriptors",
    "RegisterRDKitProperties",
    "CalcXTBDescriptors",
    "CalcXTBAtomDescriptors",
    "CalcXTBBondDescriptors",
    "CalcXTBFragmentDescriptors",
    "CalcCM5Charges",
    "CalcEHTDescriptors",
    "CalcOrbitalDescriptors",
    "AnnotateMol",
    "CalcOrbitalPopulations",
    "CalcFrontierDescriptors",
    "CalcDensityMatrix",
    "CalcShellPopulations",
    "CalcMayerBondDescriptors",
    "CalcFrontierBondDescriptors",
    "CalcPiDescriptors",
    "CalcDensityDescriptors",
    "CalcMCI",
    "CalcMCIDescriptors",
    "DESCRIPTOR_VERSION",
    "DescriptorColumns",
    "AssessDescriptorQuality",
    "CalcDescriptorReport",
    "RunDescriptorBatch",
    "ExportDescriptorReports",
    "DescriptorSensitivity",
    "DescriptorCatalog",
    "DescriptorStatus",
    "CalcReactivityDescriptors",
    "RunReactivity",
    "ReactivityJob",
    "ReactivityCheckpoint",
    "RunReactivityBatch",
    "ConformerDescriptorColumns",
    "RunConformerAnalysis",
    "ConformerGenerationResult",
    "PrepareConformers",
    "XTBOptimizationResult",
    "OptimizeXTB",
    "BuildFragmentVisualData",
    "DrawDescriptorMap",
    "PlotFragmentSummary",
    "PlotFragmentInteractions",
    "RenderFragmentAnalysis",
    "CalcESP",
    "GenerateESPGrid",
    "FitESPCharges",
    "FitRESPCharges",
    "ESPChargeFit",
    "CalcSterimol",
    "CalcBuriedVolume",
    "CalcPlaneDescriptors",
    "CalcHOMA",
    "CalcRingPuckering",
    "AnalyzeLocalGeometry",
    "RunLocalGeometryBatch",
    "SelectionSet",
    "NormalizeSelections",
]
