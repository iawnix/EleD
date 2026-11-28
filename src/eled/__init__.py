"""An RDKit extension for electronic-structure descriptors backed by xTB.

The public package is centered on :class:`XTBCalculator`, :class:`XTBResults`
and the project-built extended ``libxtb``, following RDKit's ``PascalCase``
function/method and ``camelCase`` parameter/property conventions. Datasets,
model training and application studies are maintained outside this package.
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
    CalcPiDescriptors,
    CalcShellPopulations,
)
from .descriptor_catalog import DescriptorCatalog, DescriptorStatus
from .dipro import CalcElectronicCoupling, CalcElectronicCouplingFromResults
from .electronic_descriptors import (
    CalcCM5Charges,
    CalcXTBDescriptors,
    CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors,
    CalcXTBFragmentDescriptors,
)
from .electric_field import CalcBondElectricField, CalcBondElectricFieldFromESP
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
from .esp_descriptors import (
    CalcESPSurfaceDescriptors,
    ESPSurfaceGrid,
    GenerateESPSurfaceGrid,
)
from .fod import CalcFODFromResult, RunFOD
from .fragment_eda import CalcFragmentInteraction
from .geometry import ConformerGenerationResult, PrepareConformers
from .local_geometry import (
    AnalyzeLocalGeometry,
    CalcBuriedVolume,
    CalcConeAngle,
    CalcHOMA,
    CalcPlaneDescriptors,
    CalcPyramidalization,
    CalcRingPuckering,
    CalcSterimol,
    RunLocalGeometryBatch,
)
from .local_modes import CalcLocalBondForceConstants, CalcLocalBondForceConstantsFromHessian
from .localized_orbitals import CalcLocalizedOrbitals
from .mapping import AnnotateMol
from .optimization import OptimizeXTB, XTBOptimizationResult
from .options import XTBOptions
from .orbital_descriptors import (
    CalcAtomicOrbitalEnergies,
    CalcBondOverlapDescriptors,
    CalcLinearResponseMatrix,
    CalcOrbitalWeightedFukui,
)
from .persistence import ReactivityCheckpoint, XTBDiskCache
from .populations import CalcFrontierDescriptors, CalcOrbitalPopulations
from .properties import (
    CalcEnergyDecomposition,
    CalcMultipoleDescriptors,
    CalcNuclearESP,
    CalcPolarizabilityDescriptors,
)
from .rdkit_descriptors import CalcMolXTBDescriptors, RegisterRDKitProperties
from .reactivity import CalcReactivityDescriptors, ChargeStateCalculationError, RunReactivity
from .reactivity_batch import ReactivityJob, RunReactivityBatch
from .realspace import (
    CalcBondDensityTopology,
    CalcDensityIsosurface,
    CalcNCIDescriptors,
    CalcRealSpaceProperties,
    CalcRealSpaceSurfaceDescriptors,
    CalcRealSpaceTopology,
)
from .result import XTBResults
from .reports import (
    DESCRIPTOR_VERSION,
    AssessDescriptorQuality,
    CalcDescriptorReport,
    DescriptorColumns,
    ExportDescriptorReports,
    RunDescriptorBatch,
)
from .selections import NormalizeSelections, SelectionSet
from .solvation import CalcSolvationDescriptors
from .surface_extrema import CalcDirectedSurfaceDescriptors
from .thermochemistry import CalcVibrationalThermochemistry
from .visualization import (
    BuildFragmentVisualData,
    DrawDescriptorMap,
    PlotFragmentInteractions,
    PlotFragmentSummary,
    RenderFragmentAnalysis,
)

__all__ = [
    "XTBError",
    "XTBResults",
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
    "AnnotateMol",
    "CalcOrbitalPopulations",
    "CalcFrontierDescriptors",
    "CalcDensityMatrix",
    "CalcShellPopulations",
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
    "ESPSurfaceGrid",
    "GenerateESPSurfaceGrid",
    "CalcESPSurfaceDescriptors",
    "CalcOrbitalWeightedFukui",
    "CalcAtomicOrbitalEnergies",
    "CalcBondOverlapDescriptors",
    "CalcLinearResponseMatrix",
    "CalcLocalizedOrbitals",
    "CalcFODFromResult",
    "RunFOD",
    "CalcBondElectricField",
    "CalcBondElectricFieldFromESP",
    "CalcRealSpaceProperties",
    "CalcRealSpaceSurfaceDescriptors",
    "CalcDensityIsosurface",
    "CalcRealSpaceTopology",
    "CalcBondDensityTopology",
    "CalcNCIDescriptors",
    "CalcDirectedSurfaceDescriptors",
    "CalcFragmentInteraction",
    "CalcElectronicCoupling",
    "CalcElectronicCouplingFromResults",
    "CalcLocalBondForceConstants",
    "CalcLocalBondForceConstantsFromHessian",
    "CalcEnergyDecomposition",
    "CalcMultipoleDescriptors",
    "CalcNuclearESP",
    "CalcPolarizabilityDescriptors",
    "CalcSolvationDescriptors",
    "CalcVibrationalThermochemistry",
    "CalcSterimol",
    "CalcBuriedVolume",
    "CalcPlaneDescriptors",
    "CalcHOMA",
    "CalcRingPuckering",
    "CalcPyramidalization",
    "CalcConeAngle",
    "AnalyzeLocalGeometry",
    "RunLocalGeometryBatch",
    "SelectionSet",
    "NormalizeSelections",
]
