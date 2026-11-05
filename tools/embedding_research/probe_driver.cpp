// Included after the instrumented upstream Embedder.cpp. Research only.
// All hooks observe const data; no coordinates, constraints, or RNG are changed.
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>

namespace RDKit::DGeomHelpers::probe {
unsigned int attempt = 0;
std::map<std::string, unsigned int> snapshots;

void number(double value) {
  if (std::isfinite(value)) std::cout << value;
  else std::cout << "null";
}

template <typename Container>
void array(const Container &values) {
  std::cout << '[';
  bool first = true;
  for (const auto &v : values) {
    if (!first) std::cout << ',';
    first = false;
    number(v);
  }
  std::cout << ']';
}

void chiralSets(const DistGeom::VECT_CHIRALSET &sets) {
  std::cout << '[';
  bool first = true;
  for (const auto &c : sets) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << "{\"atoms\":";
    array(std::array{c->d_idx0, c->d_idx1, c->d_idx2, c->d_idx3, c->d_idx4});
    std::cout << ",\"lower\":" << c->getLowerVolumeBound()
              << ",\"upper\":" << c->getUpperVolumeBound()
              << ",\"flags\":" << c->d_structureFlags << '}';
  }
  std::cout << ']';
}

void setup(const detail::EmbedArgs &e) {
  std::cout << "{\"kind\":\"constraints\",\"bounds\":[";
  bool first = true;
  for (unsigned int i = 0; i < e.mmat->numRows(); ++i) {
    for (unsigned int j = i + 1; j < e.mmat->numRows(); ++j) {
      if (!first) std::cout << ',';
      first = false;
      std::cout << '[' << i << ',' << j << ',' << e.mmat->getLowerBound(i, j)
                << ',' << e.mmat->getUpperBound(i, j) << ']';
    }
  }
  std::cout << "],\"chiral_sets\":";
  chiralSets(*e.chiralCenters);
  std::cout << ",\"tetrahedral_sets\":";
  chiralSets(*e.tetrahedralCarbons);
  std::cout << ",\"torsion_atoms\":[";
  first = true;
  for (const auto &atoms : e.etkdgDetails->expTorsionAtoms) {
    if (!first) std::cout << ',';
    first = false;
    array(atoms);
  }
  std::cout << "]}" << std::endl;
}

void snapshot(const RDGeom::PointPtrVect &p, const detail::EmbedArgs &e,
              const char *stage) {
  // Initial eigendecomposition failures do not have usable coordinates.
  if (std::string(stage) == "INITIAL_COORDS") return;
  if (++snapshots[stage] > 3) return;
  std::cout << "{\"kind\":\"snapshot\",\"attempt\":" << attempt
            << ",\"stage\":\"" << stage << "\",\"positions\":[";
  for (std::size_t i = 0; i < p.size(); ++i) {
    if (i) std::cout << ',';
    std::cout << '[';
    for (unsigned int d = 0; d < p[i]->dimension(); ++d) {
      if (d) std::cout << ',';
      number((*p[i])[d]);
    }
    std::cout << ']';
  }
  std::cout << "],\"chiral_checks\":[";
  bool first = true;
  for (const auto &c : *e.chiralCenters) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << "{\"volume\":";
    number(DistGeom::calcChiralVolume(c->d_idx1, c->d_idx2, c->d_idx3, c->d_idx4, p));
    std::cout << ",\"center_inside\":" << (_centerInVolume(c, p) ? "true" : "false") << '}';
  }
  std::cout << "],\"tetrahedral_checks\":[";
  first = true;
  for (const auto &c : *e.tetrahedralCarbons) {
    if (!first) std::cout << ',';
    first = false;
    std::cout << "{\"volume_ok\":" << (_volumeTest(c, p) ? "true" : "false")
              << ",\"center_inside\":"
              << (_centerInVolume(c, p, TETRAHEDRAL_CENTERINVOLUME_TOL) ? "true" : "false") << '}';
  }
  std::cout << "]}" << std::endl;
}
}  // namespace RDKit::DGeomHelpers::probe

int main(int argc, char **argv) {
  if (argc != 7) {
    std::cerr << "usage: embedding_probe SMILES PRESET MODE SEED MAX_ITERATIONS UFF_OR_MMFF\n";
    return 2;
  }
  try {
    using namespace RDKit;
    using namespace RDKit::DGeomHelpers;
    std::cout << std::setprecision(17);
    std::unique_ptr<RWMol> mol(SmilesToMol(argv[1]));
    if (!mol) throw std::invalid_argument("invalid SMILES");
    MolOps::addHs(*mol);
    std::string preset(argv[2]), mode(argv[3]);
    EmbedParameters p;
    if (preset == "ETKDGv3") p = ETKDGv3;
    else if (preset == "ETKDGv4") p = ETKDGv4;
    else if (preset == "srETKDGv3") p = srETKDGv3;
    else throw std::invalid_argument("unknown preset");
    if (mode == "internal") p.initialEmbeddingMode = InitialEmbeddingMode::INTERNAL_COORDINATE_EMBEDDING;
    else if (mode == "random") p.useRandomCoords = true;
    else if (mode != "dg") throw std::invalid_argument("unknown mode");
    const std::string forceField(argv[6]);
    if (forceField == "MMFF") p.embedForceField = EmbedFF::MMFF;
    else if (forceField != "UFF") throw std::invalid_argument("unknown embedding force field");
    p.randomSeed = std::stoi(argv[4]);
    p.maxIterations = std::stoi(argv[5]);
    if (p.randomSeed < 0 || p.maxIterations < 1) throw std::invalid_argument("invalid budget or seed");
    p.numThreads = 1;
    p.trackFailures = true;
    p.enforceChirality = true;
    p.callback = [](unsigned int n) { probe::attempt = n; };
    int cid = -1;
    bool optimizerException = false;
    try {
      cid = EmbedMolecule(*mol, p);
    } catch (const std::exception &e) {
      if (std::string(e.what()).find("bad direction in linearSearch") == std::string::npos) throw;
      optimizerException = true;
      std::cerr << e.what() << std::endl;
    }
    std::cout << "{\"kind\":\"result\",\"optimizer_exception\":"
              << (optimizerException ? "true" : "false")
              << ",\"conformer_id\":" << cid << ",\"failure_counts\":";
    probe::array(p.failures);
    std::cout << ",\"atomic_numbers\":[";
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
      if (i) std::cout << ',';
      std::cout << mol->getAtomWithIdx(i)->getAtomicNum();
    }
    std::cout << "]}" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    return 1;
  }
}
