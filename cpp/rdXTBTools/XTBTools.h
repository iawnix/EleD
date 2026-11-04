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
  std::string libraryPath;
  std::string parameterDirectory;
  int maxIterations = 250;
  double accuracy = 1.0;
  double electronicTemperature = 300.0;
};

// Atomic quantities follow ROMol atom order. Energies and matrix meanings are
// deliberately distinct from EHTTools::EHTResults. Check runMol's return value
// (or success) before using numerical fields. Arrays own all their storage.
struct XTBResults {
  bool success = false;
  std::string errorMessage;
  std::string method;
  int conformerId = -1;
  unsigned int numAtoms = 0;
  int numOrbitals = 0;
  int numElectrons = 0;       // modeled electrons, sum of orbital occupancies
  int numTotalElectrons = 0;  // sum of nuclear charges minus molecular charge
  int numUnpairedElectrons = 0;
  int charge = 0;
  int apiVersion = 0;
  int maxIterations = 250;
  double accuracy = 1.0;
  double electronicTemperature = 300.0;
  double effectiveElectronicTemperature = 300.0;
  bool gfn0StateFixVerified = false;
  double totalEnergyHartree = std::numeric_limits<double>::quiet_NaN();
  std::string libraryPath;
  std::string librarySHA256;
  std::string parameterFile;
  std::string parameterSHA256;
  std::string sourceCommit;
  std::vector<double> orbitalEnergiesEV;
  std::vector<double> orbitalOccupancies;
  std::vector<double> atomicCharges;         // elementary charge
  std::vector<double> dipoleAtomicUnits;     // e * Bohr, length 3
  std::vector<double> wibergBondOrders;      // row-major [atom][atom]
  bool basisExportVerified = false;
  std::string inputProvenanceJSON;  // Python wrapper snapshot; empty for standalone C++.
  std::vector<double> overlapMatrix;
  std::vector<int> aoAtomMap, aoShellMap, shellAtomMap, shellAngularMomenta;
  std::vector<double> orbitalCoefficients;   // row-major [AO][MO], optional
};

// Gas-phase single point. Invalid molecular input/options throw
// std::invalid_argument; unavailable/incompatible libraries throw
// std::runtime_error. xTB calculation failures return false and an errorMessage.
// Calls through this C++ library are serialized. The molecule is never changed.
bool runMol(const ROMol &mol, XTBResults &result,
            const XTBOptions &options = XTBOptions{});

}  // namespace RDKit::XTBTools
