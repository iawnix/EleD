#include "XTBTools.h"
#include <GraphMol/RWMol.h>
#include <GraphMol/Atom.h>
#include <GraphMol/Conformer.h>
#include <cmath>
#include <future>
#include <iostream>
#include <stdexcept>

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}

int main() {
  try {
    RDKit::RWMol mol;
    for (int z : {8, 1, 1}) {
      auto *atom = new RDKit::Atom(z);
      atom->setNoImplicit(true);
      mol.addAtom(atom, false, true);
    }
    mol.addBond(0, 1, RDKit::Bond::SINGLE);
    mol.addBond(0, 2, RDKit::Bond::SINGLE);
    mol.updatePropertyCache();
    auto *conf = new RDKit::Conformer(3);
    conf->set3D(true);
    conf->setId(7);
    conf->setAtomPos(0, RDGeom::Point3D(0, 0, 0));
    conf->setAtomPos(1, RDGeom::Point3D(.758, .587, 0));
    conf->setAtomPos(2, RDGeom::Point3D(-.758, .587, 0));
    mol.addConformer(conf, false);
    RDKit::XTBTools::XTBOptions options;
    options.confId = 7;
    options.keepOrbitalCoefficients = true;
    RDKit::XTBTools::XTBResults before, after, intermediate;
    require(RDKit::XTBTools::runMol(mol, before, options), "GFN2 calculation failed");
    require(before.numAtoms == 3 && before.numElectrons == 8 && before.numTotalElectrons == 10 && before.conformerId == 7, "Input identity or electron semantics changed");
    require(before.cm5Charges.size() == 3, "GFN2 CM5 charges missing");
    require(std::abs(before.cm5Charges[0] + before.cm5Charges[1] + before.cm5Charges[2]) < 1e-12, "CM5 charge conservation failed");
    require(std::abs(before.cm5Charges[0] - before.atomicCharges[0]) > 1e-4, "CM5 correction was zero on water");
    require(before.gradientHartreePerBohr.size() == 9, "Gradient dimensions incorrect");
    for (double value : before.gradientHartreePerBohr)
      require(std::isfinite(value), "Gradient contains a non-finite value");
    require(before.orbitalCoefficients.size() == static_cast<size_t>(before.numOrbitals * before.numOrbitals), "Coefficient dimensions incorrect");
    options.method = "GFN0-xTB";
    if (before.gfn0StateFixVerified) {
      require(RDKit::XTBTools::runMol(mol, intermediate, options), "GFN0 calculation failed");
      require(intermediate.cm5Charges.empty(), "GFN0 should not expose CM5 charges");
    } else {
      bool rejected = false;
      try { RDKit::XTBTools::runMol(mol, intermediate, options); }
      catch (const std::invalid_argument &) { rejected = true; }
      require(rejected, "Unverified GFN0 was not rejected");
    }
    options.method = "GFN2-xTB";
    require(RDKit::XTBTools::runMol(mol, after, options), "GFN2 after GFN0 failed");
    require(std::abs(after.totalEnergyHartree - before.totalEnergyHartree) < 1e-10, "GFN0 contaminated GFN2");
    require(before.orbitalEnergiesEV == after.orbitalEnergiesEV, "GFN2 orbitals changed");
    options.keepOverlap = true;
    if (before.basisExportVerified) {
      require(RDKit::XTBTools::runMol(mol, intermediate, options), "Basis calculation failed");
      const int n = intermediate.numOrbitals;
      require(intermediate.overlapMatrix.size() == static_cast<size_t>(n*n), "Overlap dimensions incorrect");
      require(intermediate.aoAtomMap.size() == static_cast<size_t>(n), "AO map dimensions incorrect");
      const auto &s = intermediate.overlapMatrix;
      const auto &c = intermediate.orbitalCoefficients;
      for (int a = 0; a < n; ++a) {
        require(intermediate.aoAtomMap[a] >= 0 && intermediate.aoAtomMap[a] < 3, "Invalid atom map");
        require(intermediate.shellAtomMap.at(intermediate.aoShellMap.at(a)) == intermediate.aoAtomMap[a], "Inconsistent shell map");
        for (int b = 0; b < n; ++b) {
          double value = 0.;
          for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j)
            value += c[i*n+a] * s[i*n+j] * c[j*n+b];
          require(std::abs(value - (a == b)) < 1e-8, "C^T S C is not identity");
        }
      }
    } else {
      bool rejected = false;
      try { RDKit::XTBTools::runMol(mol, intermediate, options); }
      catch (const std::invalid_argument &) { rejected = true; }
      require(rejected, "Unverified basis export was not rejected");
    }
    options.espPointsAngstrom = {2., 3., 4., -3., 2., 4.};
    require(RDKit::XTBTools::runMol(mol, intermediate, options), "ESP sampling failed");
    require(intermediate.espValuesAU.size() == 2 && intermediate.espPointsAngstrom == options.espPointsAngstrom,
            "ESP point/value dimensions changed");
    require(std::abs(intermediate.espValuesAU[0]-intermediate.espValuesAU[1]) > 1e-6,
            "ESP did not depend on sample coordinates");
    auto sampled = intermediate;
    options.espPointsAngstrom.clear();
    require(RDKit::XTBTools::runMol(mol, intermediate, options), "Calculation after ESP failed");
    require(intermediate.espValuesAU.empty() && sampled.espValuesAU.size() == 2, "ESP result ownership failed");
    options.keepOverlap = false;
    options.maxIterations = 1;
    require(!RDKit::XTBTools::runMol(mol, intermediate, options), "Expected SCF failure");
    require(!intermediate.errorMessage.empty(), "Failure has no diagnostic");
    options.maxIterations = 250;
    std::vector<std::future<double>> pending;
    for (int i = 0; i < 4; ++i) {
      pending.push_back(std::async(std::launch::async, [&mol, options] {
        RDKit::XTBTools::XTBResults r;
        if (!RDKit::XTBTools::runMol(mol, r, options)) throw std::runtime_error(r.errorMessage);
        return r.totalEnergyHartree;
      }));
    }
    for (auto &task : pending)
      require(std::abs(task.get() - before.totalEnergyHartree) < 1e-10, "Threaded C++ calculation differs");
    options.accuracy = 0.1;
    options.maxIterations = 333;
    options.electronicTemperature = 0.;
    require(RDKit::XTBTools::runMol(mol, intermediate, options), "Nondefault settings failed");
    require(intermediate.accuracy == 0.1 && intermediate.maxIterations == 333 &&
            intermediate.electronicTemperature == 0. && intermediate.effectiveElectronicTemperature == 1e-6,
            "Numerical settings were not retained");
    options.accuracy = 1e-5;
    bool clampingRejected = false;
    try { RDKit::XTBTools::runMol(mol, intermediate, options); }
    catch (const std::invalid_argument &) { clampingRejected = true; }
    require(clampingRejected, "Accuracy requiring engine clamping was accepted");
    options.accuracy = 0.1;
    options.electronicTemperature = 300.;
    auto states = RDKit::XTBTools::runChargeStates(mol, options, 1, 1);
    require(states.size() == 3 && states.back().success, "Reference restart triplet failed");
    require(states[0].initialization == "independent" && states[1].initialization == "reference" &&
            states[2].initialization == "reference", "Wrong initialization provenance");
    require(states[1].numElectrons == states[0].numElectrons - 1 &&
            states[2].numElectrons == states[0].numElectrons + 1, "Restart retained reference electron count");
    auto repeated = RDKit::XTBTools::runChargeStates(mol, options, 1, 1);
    require(states[0].atomicCharges == repeated[0].atomicCharges &&
            states[1].atomicCharges == repeated[1].atomicCharges &&
            states[2].atomicCharges == repeated[2].atomicCharges, "Restart contaminated later calculations");
    options.maxIterations = 1;
    auto failed = RDKit::XTBTools::runChargeStates(mol, options, 1, 1);
    require(failed.size() == 1 && !failed[0].success, "Failed reference did not stop the sequence");
    options.maxIterations = 250;
    auto recovered = RDKit::XTBTools::runChargeStates(mol, options, 1, 1, false);
    require(recovered.size() == 3 && recovered.back().success &&
            recovered.back().initialization == "independent", "Independent initialization recovery failed");
    require(mol.getNumAtoms() == 3 && mol.getConformer(7).getAtomPos(1).x == .758, "Input molecule changed");
    // Extended properties ABI & implicit solvation checks.
    require(before.propertiesExportVerified, "Extended properties ABI was not verified");
    require(before.energyDecompositionFlags != 0 &&
            std::isfinite(before.electronicEnergyHartree) &&
            std::isfinite(before.repulsionEnergyHartree) &&
            std::isfinite(before.anisotropicElectrostaticEnergyHartree) &&
            std::isfinite(before.anisotropicXCEnergyHartree) &&
            !std::isfinite(before.halogenBondEnergyHartree) &&
            !std::isfinite(before.solvationEnergyHartree) &&
            std::isfinite(before.fermiLevelAlphaHartree) &&
            std::abs((before.electronicEnergyHartree + before.repulsionEnergyHartree) -
                     before.totalEnergyHartree) < 1e-10,
            "GFN2 energy decomposition failed");
    require(before.polarizabilityAvailable && before.molecularPolarizabilityAU > 0.0 &&
            before.atomicPolarizabilitiesAU.size() == 3 && before.atomicC6AU.size() == 3,
            "GFN2 polarizability missing");
    require(std::abs(before.atomicPolarizabilitiesAU[0] + before.atomicPolarizabilitiesAU[1] +
                     before.atomicPolarizabilitiesAU[2] - before.molecularPolarizabilityAU) < 1e-10,
            "GFN2 atomic polarizability sum mismatch");
    require(before.multipolesAvailable && before.atomicDipolesAU.size() == 9 &&
            before.atomicQuadrupolesAU.size() == 18 && before.molecularQuadrupoleAU.size() == 6,
            "GFN2 multipoles missing");
    require(std::abs(before.molecularQuadrupoleAU[0] + before.molecularQuadrupoleAU[2] +
                     before.molecularQuadrupoleAU[5]) < 1e-12,
            "Molecular quadrupole is not traceless");
    require(before.nuclearESPAvailable && before.nuclearESPAU.size() == 3 &&
            before.nuclearESPAU[0] < before.nuclearESPAU[1] &&
            std::abs(before.nuclearESPAU[1] - before.nuclearESPAU[2]) < 1e-10,
            "GFN2 nuclear ESP failed symmetry/ordering check");
    // GFN1: multipoles available, polarizability unavailable, halogenBond available, AES NaN.
    RDKit::XTBTools::XTBOptions gfn1Opts;
    gfn1Opts.confId = 7;
    gfn1Opts.method = "GFN1-xTB";
    RDKit::XTBTools::XTBResults gfn1Res;
    require(RDKit::XTBTools::runMol(mol, gfn1Res, gfn1Opts), "GFN1 calculation failed");
    require(!gfn1Res.polarizabilityAvailable && gfn1Res.atomicPolarizabilitiesAU.empty() &&
            !std::isfinite(gfn1Res.molecularPolarizabilityAU), "GFN1 should not expose D4 polarizability");
    require(gfn1Res.multipolesAvailable && gfn1Res.nuclearESPAvailable &&
            std::isfinite(gfn1Res.halogenBondEnergyHartree) &&
            !std::isfinite(gfn1Res.anisotropicElectrostaticEnergyHartree),
            "GFN1 property availability mismatch");
    // GFN0: polarizability available, multipoles unavailable, nuclear ESP available.
    RDKit::XTBTools::XTBOptions gfn0Opts;
    gfn0Opts.confId = 7;
    gfn0Opts.method = "GFN0-xTB";
    RDKit::XTBTools::XTBResults gfn0Res;
    require(RDKit::XTBTools::runMol(mol, gfn0Res, gfn0Opts), "GFN0 calculation failed");
    require(gfn0Res.polarizabilityAvailable && !gfn0Res.multipolesAvailable &&
            gfn0Res.nuclearESPAvailable, "GFN0 property availability mismatch");
    // Implicit solvation (ALPB water) + invalid solvation option guards.
    RDKit::XTBTools::XTBOptions solvOpts;
    solvOpts.confId = 7;
    solvOpts.solvent = "water";
    solvOpts.solvationModel = "alpb";
    solvOpts.solvationGrid = 230;
    RDKit::XTBTools::XTBResults solvRes;
    require(RDKit::XTBTools::runMol(mol, solvRes, solvOpts), "ALPB water single point failed");
    require(solvRes.solvent == "water" && solvRes.solvationModel == "alpb" && solvRes.solvationGrid == 230,
            "Solvation provenance mismatch");
    require(std::isfinite(solvRes.solvationEnergyHartree) && solvRes.solvationEnergyHartree < 0.0 &&
            std::isfinite(solvRes.solvationBornEnergyHartree) &&
            std::abs(solvRes.solvationEnergyHartree -
                     (solvRes.solvationBornEnergyHartree + solvRes.solvationSASAEnergyHartree +
                      solvRes.solvationHBEnergyHartree + solvRes.solvationShiftEnergyHartree)) < 1e-10,
            "ALPB solvation energy decomposition mismatch");
    bool badSolvRejected = false;
    try {
      auto badOpts = solvOpts;
      badOpts.solvationModel = "invalid_model";
      RDKit::XTBTools::runMol(mol, intermediate, badOpts);
    } catch (const std::invalid_argument &) { badSolvRejected = true; }
    require(badSolvRejected, "Invalid solvationModel was not rejected");
    std::cout << "C++ GFN2, GFN0 guard/state, properties ABI, solvation, SCF failure recovery, and serialization passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
