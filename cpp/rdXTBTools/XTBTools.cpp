#include "XTBTools.h"

#include <GraphMol/Atom.h>
#include <GraphMol/Conformer.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <openssl/evp.h>
#include <xtb.h>
#include <dlfcn.h>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <memory>
#include <mutex>
#include <numeric>
#include <sstream>
#include <stdexcept>

#ifndef RDXTB_DEFAULT_LIBRARY
#define RDXTB_DEFAULT_LIBRARY ""
#endif

namespace RDKit::XTBTools {
namespace {
namespace fs = std::filesystem;
constexpr double angstromToBohr = 1.8897261254578281;
constexpr double hartreeToEV = 27.211386245988;
constexpr const char *upstreamGFN0Commit = "33fe381adde775c4992610c0f90d8914b7572bb5";
constexpr const char *upstreamD4SHA256 = "5c4022078296e46dfe6a21a36505ec07402bb1eee88448d735c588676dd6e61c";
constexpr const char *d4PatchSHA256 =
    "2f503f401c410b1456824efe906b89d5f932e8a8f9c63f817669a259db04c859";
constexpr const char *basisPatchSHA256 = "30f93f19b52a688727cb8af30cd0b856fe47977a96bef3c52436be707bcb75a2";
std::mutex calculationMutex;

std::string sha256(const fs::path &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input) throw std::runtime_error("Cannot read " + path.string());
  std::unique_ptr<EVP_MD_CTX, decltype(&EVP_MD_CTX_free)> ctx(EVP_MD_CTX_new(), EVP_MD_CTX_free);
  if (!ctx || EVP_DigestInit_ex(ctx.get(), EVP_sha256(), nullptr) != 1)
    throw std::runtime_error("Cannot initialize SHA-256");
  std::array<char, 65536> block{};
  while (input.read(block.data(), block.size()) || input.gcount()) {
    if (EVP_DigestUpdate(ctx.get(), block.data(), input.gcount()) != 1)
      throw std::runtime_error("Cannot hash " + path.string());
  }
  if (input.bad()) throw std::runtime_error("Cannot read " + path.string());
  std::array<unsigned char, EVP_MAX_MD_SIZE> bytes{};
  unsigned int size = 0;
  if (EVP_DigestFinal_ex(ctx.get(), bytes.data(), &size) != 1)
    throw std::runtime_error("Cannot finish SHA-256");
  std::ostringstream output;
  for (unsigned int i = 0; i < size; ++i)
    output << std::hex << std::setfill('0') << std::setw(2) << static_cast<int>(bytes[i]);
  return output.str();
}

struct EngineError : std::runtime_error { using std::runtime_error::runtime_error; };

struct CloseLibrary {
  void operator()(void *handle) const { dlclose(handle); }
};

struct API {
  std::unique_ptr<void, CloseLibrary> handle;
  std::string path, hash;
  int version = 0;
  int (*getBasisAPIVersion)() = nullptr;
  void (*setKeepOverlap)(xtb_TEnvironment, xtb_TResults, bool) = nullptr;
  void (*getBasisSize)(xtb_TEnvironment, xtb_TResults, int *) = nullptr;
  void (*getOverlap)(xtb_TEnvironment, xtb_TResults, int, double *) = nullptr;
  void (*getBasisMap)(xtb_TEnvironment, xtb_TResults, int, int, int *, int *, int *, int *) = nullptr;
#define XTB_FUNCTION(name) decltype(&xtb_##name) name = nullptr
  XTB_FUNCTION(getAPIVersion);
  XTB_FUNCTION(newEnvironment); XTB_FUNCTION(delEnvironment);
  XTB_FUNCTION(checkEnvironment); XTB_FUNCTION(getError); XTB_FUNCTION(setVerbosity);
  XTB_FUNCTION(newMolecule); XTB_FUNCTION(delMolecule);
  XTB_FUNCTION(newCalculator); XTB_FUNCTION(delCalculator);
  XTB_FUNCTION(loadGFN0xTB); XTB_FUNCTION(loadGFN1xTB); XTB_FUNCTION(loadGFN2xTB);
  XTB_FUNCTION(setMaxIter); XTB_FUNCTION(setAccuracy); XTB_FUNCTION(setElectronicTemp);
  XTB_FUNCTION(newResults); XTB_FUNCTION(delResults); XTB_FUNCTION(singlepoint);
  XTB_FUNCTION(getEnergy); XTB_FUNCTION(getNao); XTB_FUNCTION(getDipole);
  XTB_FUNCTION(getCharges); XTB_FUNCTION(getBondOrders);
  XTB_FUNCTION(getOrbitalEigenvalues); XTB_FUNCTION(getOrbitalOccupations);
  XTB_FUNCTION(getOrbitalCoefficients);
#undef XTB_FUNCTION

  template <typename T> void bind(T &function, const char *name) {
    function = reinterpret_cast<T>(dlsym(handle.get(), name));
    if (!function) throw std::runtime_error(std::string("xTB library lacks ") + name);
  }
  explicit API(const std::string &requested) {
    handle.reset(dlopen(requested.c_str(), RTLD_NOW | RTLD_LOCAL));
    if (!handle) throw std::runtime_error("Cannot load xTB: " + std::string(dlerror()));
    bind(getAPIVersion, "xtb_getAPIVersion");
    version = getAPIVersion();
    if (version / 10000 != 1 && version / 10000 != 2)
      throw std::runtime_error("Unsupported xTB C API version " + std::to_string(version));
    Dl_info info{};
    if (!dladdr(reinterpret_cast<void *>(getAPIVersion), &info) || !info.dli_fname)
      throw std::runtime_error("Cannot identify the loaded xTB library");
    path = fs::canonical(info.dli_fname).string();
    hash = sha256(path);
#define BIND(name) bind(name, "xtb_" #name)
    BIND(newEnvironment); BIND(delEnvironment); BIND(checkEnvironment); BIND(getError);
    BIND(setVerbosity); BIND(newMolecule); BIND(delMolecule);
    BIND(newCalculator); BIND(delCalculator);
    BIND(loadGFN0xTB); BIND(loadGFN1xTB); BIND(loadGFN2xTB);
    BIND(setMaxIter); BIND(setAccuracy); BIND(setElectronicTemp);
    BIND(newResults); BIND(delResults); BIND(singlepoint);
    BIND(getEnergy); BIND(getNao); BIND(getDipole); BIND(getCharges); BIND(getBondOrders);
    BIND(getOrbitalEigenvalues); BIND(getOrbitalOccupations); BIND(getOrbitalCoefficients);
#undef BIND
#define OPTIONAL(name) name = reinterpret_cast<decltype(name)>(dlsym(handle.get(), "xtb_rdkit_" #name))
    OPTIONAL(getBasisAPIVersion); OPTIONAL(setKeepOverlap); OPTIONAL(getBasisSize);
    OPTIONAL(getOverlap); OPTIONAL(getBasisMap);
#undef OPTIONAL
  }
  void check(xtb_TEnvironment env, const char *operation) const {
    if (!checkEnvironment(env)) return;
    std::array<char, 4096> text{};
    const int size = text.size();
    getError(env, text.data(), &size);
    throw EngineError(std::string(operation) + ": " + text.data());
  }
};

// Persistent handles avoid unloading/reloading Fortran runtimes between calls.
// Access is protected by calculationMutex. Replacing a library requires restart.
std::map<std::string, std::shared_ptr<API>> libraries;
API &loadAPI(std::string requested) {
  if (requested.empty()) {
    const char *environment = std::getenv("XTB_LIBRARY");
    requested = environment && *environment ? environment : RDXTB_DEFAULT_LIBRARY;
    if (requested.empty()) requested = "libxtb.so";
  }
  if (fs::exists(requested)) requested = fs::canonical(requested).string();
  auto found = libraries.find(requested);
  if (found == libraries.end()) found = libraries.emplace(requested, std::make_shared<API>(requested)).first;
  return *found->second;
}

bool verifyGFN0(const API &api, std::string &commit) {
  try {
    const auto prefix = fs::path(api.path).parent_path().parent_path();
    boost::property_tree::ptree record;
    boost::property_tree::read_json((prefix / "share/xtb/rdkit-eth/build.json").string(), record);
    const bool upstream = record.get<std::string>("source.commit", "") == upstreamGFN0Commit &&
      record.get<std::string>("source.gfn0_fix.kind", "") == "upstream" &&
      record.get<std::string>("source.gfn0_fix.commit", "") == upstreamGFN0Commit &&
      record.get<std::string>("source.gfn0_fix.file_sha256", "") == upstreamD4SHA256;
    const bool legacy = record.get<std::string>("source.patch_sha256", "") == d4PatchSHA256;
    if (record.get<int>("schema_version", 0) != 1 ||
        record.get<std::string>("capabilities.native_gfn0", "") != "true" ||
        (!upstream && !legacy) ||
        fs::canonical(prefix / "lib/libxtb.so") != fs::path(api.path)) return false;
    // A slash is literal in ptree paths, but dots in libxtb.so are not.
    const auto &files = record.get_child("files");
    const auto entry = files.find("lib/libxtb.so");
    if (entry == files.not_found() || entry->second.get_value<std::string>() != api.hash ||
        sha256(api.path) != api.hash) return false;
    commit = record.get<std::string>("source.commit", "");
    return true;
  } catch (const std::exception &) { return false; }
}

bool verifyBasis(const API &api) {
  if (!api.getBasisAPIVersion || !api.setKeepOverlap || !api.getBasisSize || !api.getOverlap || !api.getBasisMap)
    return false;
  try {
    const auto prefix = fs::path(api.path).parent_path().parent_path();
    boost::property_tree::ptree record;
    boost::property_tree::read_json((prefix / "share/xtb/rdkit-eth/build.json").string(), record);
    const auto &files = record.get_child("files");
    const auto entry = files.find("lib/libxtb.so");
    return record.get<int>("schema_version", 0) == 1 &&
      record.get<int>("capabilities.orbital_basis_api", 0) == 1 &&
      record.get<int>("source.basis_extension.api_version", 0) == 1 &&
      record.get<std::string>("source.basis_extension.patch_sha256", "") == basisPatchSHA256 &&
      fs::canonical(prefix / "lib/libxtb.so") == fs::path(api.path) &&
      entry != files.not_found() && entry->second.get_value<std::string>() == api.hash &&
      sha256(api.path) == api.hash && api.getBasisAPIVersion() == 1;
  } catch (const std::exception &) { return false; }
}

std::string parameterFile(const API &api, const XTBOptions &options, int method) {
  const std::string filename = "param_gfn" + std::to_string(method) + "-xtb.txt";
  std::vector<fs::path> roots;
  if (!options.parameterDirectory.empty()) {
    auto file = fs::path(options.parameterDirectory) / filename;
    if (!fs::is_regular_file(file)) throw std::invalid_argument("Parameter file not found: " + file.string());
    return fs::canonical(file).string();
  }
  if (const char *paths = std::getenv("XTBPATH")) {
    std::istringstream stream(paths);
    for (std::string part; std::getline(stream, part, ':');)
      if (!part.empty()) roots.emplace_back(part);
  }
  roots.push_back(fs::path(api.path).parent_path().parent_path() / "share/xtb");
  for (const auto &root : roots)
    if (fs::is_regular_file(root / filename)) return fs::canonical(root / filename).string();
  if (method == 0) throw std::invalid_argument("GFN0 needs param_gfn0-xtb.txt; set parameterDirectory or XTBPATH");
  return {};  // GFN1/2 can use built-in parameters.
}

struct Handles {
  API &api;
  xtb_TEnvironment env = nullptr;
  xtb_TMolecule mol = nullptr;
  xtb_TCalculator calc = nullptr;
  xtb_TResults res = nullptr;
  explicit Handles(API &library) : api(library) {}
  ~Handles() {
    if (res) api.delResults(&res);
    if (calc) api.delCalculator(&calc);
    if (mol) api.delMolecule(&mol);
    if (env) api.delEnvironment(&env);
  }
};

int methodNumber(std::string name) {
  std::string key;
  for (unsigned char c : name)
    if (c != '-' && c != '_' && !std::isspace(c)) key += std::tolower(c);
  for (int i = 0; i <= 2; ++i) {
    const auto base = "gfn" + std::to_string(i);
    if (key == base || key == base + "xtb") return i;
  }
  throw std::invalid_argument("method must be GFN0-xTB, GFN1-xTB, or GFN2-xTB");
}

void requireFinite(const std::vector<double> &values, const char *name) {
  if (std::any_of(values.begin(), values.end(), [](double v) { return !std::isfinite(v); }))
    throw EngineError(std::string("Nonfinite ") + name);
}
}  // namespace

bool runMol(const ROMol &mol, XTBResults &result, const XTBOptions &options) {
  std::lock_guard<std::mutex> lock(calculationMutex);
  result = XTBResults{};
  const int method = methodNumber(options.method);
  if (!mol.getNumAtoms()) throw std::invalid_argument("Molecule must contain atoms");
  if (mol.getNumAtoms() > static_cast<unsigned int>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Too many atoms for the xTB C API");
  if (!mol.getNumConformers()) throw std::invalid_argument("Molecule must have a 3D conformer");
  if (options.confId < -1) throw std::invalid_argument("confId must be -1 or a conformer ID");
  if (options.maxIterations < 1 || !std::isfinite(options.accuracy) || options.accuracy <= 0 ||
      !std::isfinite(options.electronicTemperature) || options.electronicTemperature < 0)
    throw std::invalid_argument("Invalid iteration limit, accuracy, or electronic temperature");
  const Conformer *conf = nullptr;
  try { conf = &mol.getConformer(options.confId); }
  catch (const std::exception &) { throw std::invalid_argument("Conformer ID does not exist"); }
  if (!conf->is3D()) throw std::invalid_argument("Conformer must be marked as 3D");
  const int n = mol.getNumAtoms();
  std::vector<int> numbers(n);
  std::vector<double> xyz(3 * static_cast<size_t>(n));
  long long inferredCharge = 0, inferredUHF = 0, nuclearCharge = 0;
  for (const auto atom : mol.atoms()) {
    if (atom->getNumImplicitHs() || atom->getNumExplicitHs())
      throw std::invalid_argument("Hydrogen counts must be expanded to atoms; call Chem.AddHs");
    if (!atom->getAtomicNum()) throw std::invalid_argument("Dummy atoms are not supported");
    auto i = atom->getIdx();
    numbers[i] = atom->getAtomicNum();
    inferredCharge += atom->getFormalCharge();
    inferredUHF += atom->getNumRadicalElectrons();
    nuclearCharge += atom->getAtomicNum();
    auto p = conf->getAtomPos(i);
    for (unsigned int k = 0; k < 3; ++k) {
      if (!std::isfinite(p[k])) throw std::invalid_argument("Coordinates must be finite");
      xyz[3 * static_cast<size_t>(i) + k] = p[k] * angstromToBohr;
    }
  }
  const long long charge = options.charge ? *options.charge : inferredCharge;
  const long long unpaired = options.numUnpairedElectrons ? *options.numUnpairedElectrons : inferredUHF;
  if (charge < std::numeric_limits<int>::min() || charge > std::numeric_limits<int>::max() ||
      unpaired < 0 || unpaired > std::numeric_limits<int>::max())
    throw std::invalid_argument("Charge or unpaired electron count is out of range");
  result.charge = static_cast<int>(charge);
  result.numUnpairedElectrons = static_cast<int>(unpaired);
  const long long electrons = nuclearCharge - result.charge;
  if (electrons < 1 || electrons > std::numeric_limits<int>::max() ||
      result.numUnpairedElectrons < 0 || result.numUnpairedElectrons > electrons ||
      (electrons - result.numUnpairedElectrons) % 2)
    throw std::invalid_argument("Charge and unpaired electrons are inconsistent with electron count");
  result.numTotalElectrons = static_cast<int>(electrons);
  result.method = "GFN" + std::to_string(method) + "-xTB";
  result.numAtoms = n;
  result.conformerId = conf->getId();
  auto &api = loadAPI(options.libraryPath);
  result.libraryPath = api.path;
  result.librarySHA256 = api.hash;
  result.apiVersion = api.version;
  result.gfn0StateFixVerified = verifyGFN0(api, result.sourceCommit);
  result.basisExportVerified = verifyBasis(api);
  if (options.keepOverlap && !result.basisExportVerified)
    throw std::invalid_argument("keepOverlap requires a verified native basis extension; rebuild with tools/build_patched_xtb.py");
  if (method == 0 && !result.gfn0StateFixVerified)
    throw std::invalid_argument("Native GFN0 requires a verified patched xTB build; run tools/build_patched_xtb.sh");
  result.parameterFile = parameterFile(api, options, method);
  if (!result.parameterFile.empty()) result.parameterSHA256 = sha256(result.parameterFile);
  Handles h(api);
  try {
    h.env = api.newEnvironment();
    if (!h.env) throw EngineError("Cannot allocate xTB environment");
    api.setVerbosity(h.env, 0);
    const double charge = result.charge;
    h.mol = api.newMolecule(h.env, &n, numbers.data(), xyz.data(), &charge,
                           &result.numUnpairedElectrons, nullptr, nullptr);
    api.check(h.env, "Molecule setup");
    if (!h.mol) throw EngineError("Cannot allocate xTB molecule");
    h.calc = api.newCalculator();
    if (!h.calc) throw EngineError("Cannot allocate xTB calculator");
    auto loader = method == 0 ? api.loadGFN0xTB : method == 1 ? api.loadGFN1xTB : api.loadGFN2xTB;
    loader(h.env, h.mol, h.calc, result.parameterFile.empty() ? nullptr : result.parameterFile.data());
    api.check(h.env, "Parameter loading");
    api.setMaxIter(h.env, h.calc, options.maxIterations);
    api.setAccuracy(h.env, h.calc, options.accuracy);
    api.setElectronicTemp(h.env, h.calc, options.electronicTemperature);
    api.check(h.env, "Calculator settings");
    h.res = api.newResults();
    if (!h.res) throw EngineError("Cannot allocate xTB results");
    if (options.keepOverlap) {
      api.setKeepOverlap(h.env, h.res, true);
      api.check(h.env, "Basis retention setup");
    }
    api.singlepoint(h.env, h.mol, h.calc, h.res);
    api.check(h.env, "Single point");
    api.getEnergy(h.env, h.res, &result.totalEnergyHartree);
    api.getNao(h.env, h.res, &result.numOrbitals);
    api.check(h.env, "Energy and basis size");
    if (!std::isfinite(result.totalEnergyHartree) || result.numOrbitals < 1)
      throw EngineError("Invalid total energy or orbital count");
    const size_t nao = result.numOrbitals;
    result.atomicCharges.resize(n);
    result.dipoleAtomicUnits.resize(3);
    result.wibergBondOrders.resize(static_cast<size_t>(n) * n);
    result.orbitalEnergiesEV.resize(nao);
    result.orbitalOccupancies.resize(nao);
    api.getCharges(h.env, h.res, result.atomicCharges.data());
    api.getDipole(h.env, h.res, result.dipoleAtomicUnits.data());
    api.getBondOrders(h.env, h.res, result.wibergBondOrders.data());
    api.getOrbitalEigenvalues(h.env, h.res, result.orbitalEnergiesEV.data());
    api.getOrbitalOccupations(h.env, h.res, result.orbitalOccupancies.data());
    if (options.keepOrbitalCoefficients || options.keepOverlap) {
      std::vector<double> coefficients(nao * nao);
      api.getOrbitalCoefficients(h.env, h.res, coefficients.data());
      result.orbitalCoefficients.resize(nao * nao);
      for (size_t ao = 0; ao < nao; ++ao)
        for (size_t mo = 0; mo < nao; ++mo)
          result.orbitalCoefficients[ao * nao + mo] = coefficients[mo * nao + ao];
    }
    if (options.keepOverlap) {
      std::array<int, 3> sizes{};
      api.getBasisSize(h.env, h.res, sizes.data());
      api.check(h.env, "Basis size");
      if (sizes[0] != n || sizes[2] != result.numOrbitals || sizes[1] < 1 ||
          sizes[1] > result.numOrbitals || nao * nao > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw EngineError("Native basis dimensions do not match result");
      const int nshell = sizes[1];
      result.overlapMatrix.resize(nao * nao);
      result.aoAtomMap.resize(nao); result.aoShellMap.resize(nao);
      result.shellAtomMap.resize(nshell); result.shellAngularMomenta.resize(nshell);
      api.getOverlap(h.env, h.res, static_cast<int>(nao * nao), result.overlapMatrix.data());
      api.getBasisMap(h.env, h.res, result.numOrbitals, nshell, result.aoAtomMap.data(),
                     result.aoShellMap.data(), result.shellAtomMap.data(), result.shellAngularMomenta.data());
      api.check(h.env, "Basis extraction");
      for (size_t i=0; i<nao; ++i) for (size_t j=0; j<i; ++j)
        std::swap(result.overlapMatrix[i*nao+j], result.overlapMatrix[j*nao+i]);
      requireFinite(result.overlapMatrix, "overlap");
      for (size_t i=0; i<nao; ++i) {
        const int atom=result.aoAtomMap[i], shell=result.aoShellMap[i];
        if (atom<0 || atom>=n || shell<0 || shell>=nshell || result.shellAtomMap[shell]!=atom)
          throw EngineError("Invalid native basis mapping");
      }
    }
    api.check(h.env, "Result extraction");
    for (auto &energy : result.orbitalEnergiesEV) energy *= hartreeToEV;
    // The Fortran WBO array is column-major; expose atom rows and columns.
    for (size_t i = 0; i < static_cast<size_t>(n); ++i)
      for (size_t j = 0; j < i; ++j)
        std::swap(result.wibergBondOrders[i * n + j], result.wibergBondOrders[j * n + i]);
    requireFinite(result.atomicCharges, "atomic charges");
    requireFinite(result.dipoleAtomicUnits, "dipole");
    requireFinite(result.wibergBondOrders, "bond orders");
    requireFinite(result.orbitalEnergiesEV, "orbital energies");
    requireFinite(result.orbitalOccupancies, "occupancies");
    requireFinite(result.orbitalCoefficients, "orbital coefficients");
    const double occupied = std::accumulate(result.orbitalOccupancies.begin(), result.orbitalOccupancies.end(), 0.0);
    if (occupied < 0 || occupied > result.numTotalElectrons + 1e-5 ||
        std::abs(occupied - std::round(occupied)) > 1e-5)
      throw EngineError("Invalid orbital electron count");
    result.numElectrons = static_cast<int>(std::round(occupied));
    result.success = true;
    return true;
  } catch (const EngineError &error) {
    result.errorMessage = error.what();
    return false;
  }
}
}  // namespace RDKit::XTBTools
