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
    require(before.orbitalCoefficients.size() == static_cast<size_t>(before.numOrbitals * before.numOrbitals), "Coefficient dimensions incorrect");
    options.method = "GFN0-xTB";
    if (before.gfn0StateFixVerified) {
      require(RDKit::XTBTools::runMol(mol, intermediate, options), "GFN0 calculation failed");
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
    std::cout << "C++ GFN2, GFN0 guard/state, SCF failure recovery, and serialization passed\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
