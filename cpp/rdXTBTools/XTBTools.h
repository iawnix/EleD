#pragma once

#include <GraphMol/ROMol.h>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace RDKit::XTBTools {

struct XTBOptions {
  int confId = -1;
  std::string method = "GFN2-xTB";
  std::optional<int> charge;
  std::optional<int> numUnpairedElectrons;
  bool keepOrbitalCoefficients = false;
  bool keepOverlap = false;
  std::vector<double> espPointsAngstrom; // row-major [point][xyz]; empty disables sampling
  std::vector<double> realSpacePointsAngstrom; // row-major [point][xyz]; empty disables sampling
  double realSpaceOccupationThreshold = 1e-7;
  double realSpaceLEAWindowEV = 0.0;
  std::string libraryPath;
  std::string parameterDirectory;
  int maxIterations = 250;
  double accuracy = 1.0;
  double electronicTemperature = 300.0;
  std::string solvent;             // empty disables solvation (gas phase)
  std::string solvationModel = "alpb"; // "alpb" or "gbsa" when solvent is nonempty
  int solvationGrid = 230;         // Lebedev angular grid size when solvent is nonempty
};

// Atomic quantities follow ROMol atom order. Energies and matrix meanings are
// deliberately distinct from EHTTools::EHTResults. Check runMol's return value
// (or success) before using numerical fields. Arrays own all their storage.
struct XTBResults {
  bool success = false;
  std::string errorMessage;
  std::string failureKind;  // calculation or environment; empty on success
  std::string method;
  int conformerId = -1;
  unsigned int numAtoms = 0;
  int numOrbitals = 0;
  int numElectrons = 0;       // modeled electrons, sum of orbital occupancies
  int numTotalElectrons = 0;  // sum of nuclear charges minus molecular charge
  int numUnpairedElectrons = 0;
  int charge = 0;
  // Exact semantics of atomicCharges: GFN1/2 are SCC Mulliken charges;
  // GFN0 is the PEEQ/EEQ charge model. They are not CM5 charges.
  std::string chargeModel;
  std::vector<int> atomicNumbers;
  std::vector<double> coordinatesAngstrom; // row-major [atom][xyz]
  int apiVersion = 0;
  int maxIterations = 250;
  double accuracy = 1.0;
  double electronicTemperature = 300.0;
  double effectiveElectronicTemperature = 300.0;
  bool gfn0StateFixVerified = false;
  int scfIterations = -1;  // unavailable for GFN0 or older libraries
  double scfEnergyChangeHartree = -1., scfMixerRMS = -1.;
  double scfEnergyThresholdHartree = -1., scfMixerThreshold = -1.;
  double totalEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  std::string libraryPath;
  std::string librarySHA256;
  std::string parameterFile;
  std::string parameterSHA256;
  std::string sourceCommit;
  std::vector<double> espPointsAngstrom; // row-major [point][xyz]
  std::vector<double> espValuesAU; // Hartree/e, point order; empty when not requested
  std::vector<double> gradientHartreePerBohr; // row-major [atom][xyz]
  std::vector<double> orbitalEnergiesEV;
  std::vector<double> orbitalOccupancies;
  std::vector<double> atomicCharges;         // elementary charge
  std::vector<double> cm5Charges;            // elementary charge, GFN1/2 only; empty for GFN0
  std::vector<double> dipoleAtomicUnits;     // e * Bohr, length 3
  std::vector<double> wibergBondOrders;      // row-major [atom][atom]
  bool basisExportVerified = false;
  std::string inputProvenanceJSON;  // Python wrapper snapshot; empty for standalone C++.
  std::string initialization = "independent";
  std::vector<double> overlapMatrix;
  std::vector<int> aoAtomMap, aoShellMap, shellAtomMap, shellAngularMomenta;
  std::vector<double> orbitalCoefficients;   // row-major [AO][MO], optional

  // Implicit solvation provenance (empty solvent means gas phase).
  std::string solvent;
  std::string solvationModel;
  int solvationGrid = 0;

  // Extended properties ABI v1 (all atomic units; NaN/empty when unavailable).
  bool propertiesExportVerified = false;
  int energyDecompositionFlags = 0;
  double electronicEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double repulsionEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double electrostaticEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double anisotropicElectrostaticEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double anisotropicXCEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double dispersionEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double halogenBondEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double solvationEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double solvationBornEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double solvationSASAEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double solvationHBEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double solvationShiftEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  double fermiLevelAlphaHartree = std::numeric_limits<double>::quiet_NaN();
  double fermiLevelBetaHartree = std::numeric_limits<double>::quiet_NaN();
  double fermiLevelAlphaEV = std::numeric_limits<double>::quiet_NaN();
  double fermiLevelBetaEV = std::numeric_limits<double>::quiet_NaN();

  bool polarizabilityAvailable = false;
  double molecularPolarizabilityAU = std::numeric_limits<double>::quiet_NaN(); // Bohr^3
  double molecularC6AU = std::numeric_limits<double>::quiet_NaN();             // Hartree*Bohr^6
  double molecularC8AU = std::numeric_limits<double>::quiet_NaN();             // Hartree*Bohr^8
  std::vector<double> atomicPolarizabilitiesAU; // [atom], Bohr^3
  std::vector<double> atomicC6AU;               // [atom], Hartree*Bohr^6

  bool multipolesAvailable = false;
  std::vector<double> atomicDipolesAU;          // row-major [atom][3], e*Bohr
  std::vector<double> atomicQuadrupolesAU;      // row-major [atom][6] (xx,xy,yy,xz,yz,zz), e*Bohr^2
  std::vector<double> molecularQuadrupoleAU;    // [6] (xx,xy,yy,xz,yz,zz), e*Bohr^2

  bool nuclearESPAvailable = false;
  std::vector<double> nuclearESPAU;             // [atom], Hartree/e

  // Real-space wavefunction ABI v1 (exact Cartesian GTO evaluation).
  bool realSpaceExportVerified = false;
  std::vector<double> realSpacePointsAngstrom;    // row-major [point][xyz], Angstrom
  std::vector<double> realSpaceDensityAU;         // [point], e/Bohr^3
  std::vector<double> realSpaceALIEEV;            // [point], eV
  std::vector<double> realSpaceVirtualDensityAU;  // [point], 1/Bohr^3
  std::vector<double> realSpaceLEAEV;             // [point], eV
};

// Gas-phase single point. Invalid molecular input/options throw
// std::invalid_argument; unavailable/incompatible libraries throw
// std::runtime_error. xTB calculation failures return false and an errorMessage.
// Calls through this C++ library are serialized. The molecule is never changed.
bool runMol(const ROMol &mol, XTBResults &result,
            const XTBOptions &options = XTBOptions{});

// Fixed-geometry reference, ionized, attached states. Each neighbor starts from
// its own copy of the reference when requested. Stops at the first failed state;
// returned entries include that failure, but no unattempted states. No raw
// engine handles escape this call. Only GFN1/2; adjacent spins are explicit.
std::vector<XTBResults> runChargeStates(const ROMol &mol, const XTBOptions &options,
    int ionizedUnpaired, int attachedUnpaired, bool reuseReference = true);

// Compute xTB CM5 charges in elementary-charge units from atomic numbers (1..118),
// row-major [atom][xyz] Angstrom coordinates, and SCC Mulliken charges.
std::vector<double> calcCM5Charges(const std::vector<int> &atomicNumbers,
                                   const std::vector<double> &coordinatesAngstrom,
                                   const std::vector<double> &mullikenCharges);

}  // namespace RDKit::XTBTools
