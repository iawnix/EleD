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
constexpr const char *basisPatchSHA256 = "3c18bd5e5c418910e9873e915478458026dba6228cbe19eb0da32fde6f106e2c";
constexpr const char *scfPatchSHA256 = "d8db37846002ebc19ed6808356cbd83304cc1608915e69956c1a7afc174d18fd";
constexpr const char *propertiesPatchSHA256 = "56ea1fab2c624e37f851ce3d59d11d73d1a9abd6f9553c5464983c36c07ff817";
constexpr const char *realSpacePatchSHA256 = "2a69bc44fbebd842619ca182cf97e17fbde487c555978d4b5d7158c2c4f6857b";
constexpr std::array<int, 32> validSolvationGrids = {
    6, 14, 26, 38, 50, 74, 86, 110, 146, 170, 194, 230, 266, 302, 350, 434,
    590, 770, 974, 1202, 1454, 1730, 2030, 2354, 2702, 3074, 3470, 3890, 4334, 4802, 5294, 5810};
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
  int (*getESPAPIVersion)() = nullptr;
  void (*setKeepESP)(xtb_TEnvironment, xtb_TResults, bool) = nullptr;
  int (*getESP)(xtb_TEnvironment, xtb_TResults, int, const double *, int, double *) = nullptr;
  int (*getSCFAPIVersion)() = nullptr;
  int (*getSCFDiagnostics)(xtb_TEnvironment, xtb_TResults, int *, double *) = nullptr;
  int (*getPropertiesAPIVersion)() = nullptr;
  int (*getEnergyDecomposition)(xtb_TEnvironment, xtb_TResults, int *, int, double *) = nullptr;
  int (*getPolarizability)(xtb_TEnvironment, xtb_TResults, int, double *, double *, double *) = nullptr;
  int (*getMultipoles)(xtb_TEnvironment, xtb_TResults, int, double *, double *, double *) = nullptr;
  int (*getNuclearESP)(xtb_TEnvironment, xtb_TResults, int, double *) = nullptr;
  int (*getRealSpaceAPIVersion)() = nullptr;
  int (*evalRealSpace)(xtb_TEnvironment, xtb_TResults, int, const double *, double, double,
                       double *, double *, double *, double *) = nullptr;
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
  XTB_FUNCTION(setSolvent); XTB_FUNCTION(releaseSolvent);
  XTB_FUNCTION(setMaxIter); XTB_FUNCTION(setAccuracy); XTB_FUNCTION(setElectronicTemp);
  XTB_FUNCTION(newResults); XTB_FUNCTION(copyResults); XTB_FUNCTION(delResults); XTB_FUNCTION(singlepoint);
  XTB_FUNCTION(getEnergy); XTB_FUNCTION(getNao); XTB_FUNCTION(getDipole);
  XTB_FUNCTION(getGradient);
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
    BIND(setSolvent); BIND(releaseSolvent);
    BIND(setMaxIter); BIND(setAccuracy); BIND(setElectronicTemp);
    BIND(newResults); BIND(copyResults); BIND(delResults); BIND(singlepoint);
    BIND(getEnergy); BIND(getNao); BIND(getDipole); BIND(getGradient);
    BIND(getCharges); BIND(getBondOrders);
    BIND(getOrbitalEigenvalues); BIND(getOrbitalOccupations); BIND(getOrbitalCoefficients);
#undef BIND
#define OPTIONAL(name) name = reinterpret_cast<decltype(name)>(dlsym(handle.get(), "xtb_rdkit_" #name))
    OPTIONAL(getBasisAPIVersion); OPTIONAL(setKeepOverlap); OPTIONAL(getBasisSize);
    OPTIONAL(getOverlap); OPTIONAL(getBasisMap);
    OPTIONAL(getESPAPIVersion); OPTIONAL(setKeepESP); OPTIONAL(getESP);
    OPTIONAL(getSCFAPIVersion); OPTIONAL(getSCFDiagnostics);
    OPTIONAL(getPropertiesAPIVersion); OPTIONAL(getEnergyDecomposition);
    OPTIONAL(getPolarizability); OPTIONAL(getMultipoles); OPTIONAL(getNuclearESP);
    OPTIONAL(getRealSpaceAPIVersion); OPTIONAL(evalRealSpace);
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
    boost::property_tree::read_json((prefix / "share/xtb/eled/build.json").string(), record);
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
    boost::property_tree::read_json((prefix / "share/xtb/eled/build.json").string(), record);
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

bool verifyDiagnostics(const API &api) {
  if (!api.getSCFAPIVersion || !api.getSCFDiagnostics) return false;
  try {
    const auto prefix = fs::path(api.path).parent_path().parent_path();
    boost::property_tree::ptree record;
    boost::property_tree::read_json((prefix / "share/xtb/eled/build.json").string(), record);
    const auto &extension = record.get_child("source.diagnostics_extension");
    return record.get<int>("schema_version", 0) == 1 &&
      record.get<int>("capabilities.scf_diagnostics_api", 0) == 1 &&
      extension.get<int>("api_version", 0) == 1 &&
      extension.get<std::string>("patch_sha256", "") == scfPatchSHA256 &&
      api.getSCFAPIVersion() == 1;
  } catch (const std::exception &) { return false; }
}

bool verifyESP(const API &api) {
  if (!api.getESPAPIVersion || !api.setKeepESP || !api.getESP) return false;
  try {
    const auto prefix = fs::path(api.path).parent_path().parent_path();
    boost::property_tree::ptree record;
    boost::property_tree::read_json((prefix / "share/xtb/eled/build.json").string(), record);
    std::string commit;
    return verifyGFN0(api, commit) && record.get<int>("capabilities.esp_api", 0) == 1 &&
      record.get<int>("source.esp_extension.api_version", 0) == 1 &&
      record.get<std::string>("source.esp_extension.patch_sha256", "") ==
        "aca69861eaf5e25b5437553c1a7f3c6117091c7abc3e6b3ab900109f50b70cd0" &&
      api.getESPAPIVersion() == 1;
  } catch (const std::exception &) { return false; }
}

bool verifyProperties(const API &api) {
  if (!api.getPropertiesAPIVersion || !api.getEnergyDecomposition ||
      !api.getPolarizability || !api.getMultipoles || !api.getNuclearESP)
    return false;
  try {
    const auto prefix = fs::path(api.path).parent_path().parent_path();
    boost::property_tree::ptree record;
    boost::property_tree::read_json((prefix / "share/xtb/eled/build.json").string(), record);
    std::string commit;
    return verifyGFN0(api, commit) && record.get<int>("capabilities.properties_api", 0) == 1 &&
      record.get<int>("source.properties_extension.api_version", 0) == 1 &&
      record.get<std::string>("source.properties_extension.patch_sha256", "") == propertiesPatchSHA256 &&
      api.getPropertiesAPIVersion() == 1;
  } catch (const std::exception &) { return false; }
}

bool verifyRealSpace(const API &api) {
  if (!api.getRealSpaceAPIVersion || !api.evalRealSpace) return false;
  try {
    const auto prefix = fs::path(api.path).parent_path().parent_path();
    boost::property_tree::ptree record;
    boost::property_tree::read_json((prefix / "share/xtb/eled/build.json").string(), record);
    std::string commit;
    if (!verifyGFN0(api, commit)) return false;
    const int cap = record.get<int>("capabilities.realspace_api", 0);
    const int ver = record.get<int>("source.realspace_extension.api_version", 0);
    const std::string sha = record.get<std::string>("source.realspace_extension.patch_sha256", "");
    const int live = api.getRealSpaceAPIVersion();
    if (cap == 2 && ver == 2 && live == 2 && sha == realSpacePatchSHA256) return true;
    return false;
  } catch (const std::exception &) { return false; }
}

std::string trimString(const std::string &value) {
  size_t start = 0;
  while (start < value.size() && std::isspace(static_cast<unsigned char>(value[start]))) ++start;
  size_t end = value.size();
  while (end > start && std::isspace(static_cast<unsigned char>(value[end - 1]))) --end;
  return value.substr(start, end - start);
}

std::string normalizeSolvationModel(const std::string &model) {
  std::string out;
  for (unsigned char c : trimString(model)) out += static_cast<char>(std::tolower(c));
  if (out.empty()) out = "alpb";
  if (out != "alpb" && out != "gbsa")
    throw std::invalid_argument("solvationModel must be 'alpb' or 'gbsa'");
  return out;
}

int validateSolvationGrid(int grid) {
  if (grid == 0) return 230;
  if (std::find(validSolvationGrids.begin(), validSolvationGrids.end(), grid) == validSolvationGrids.end())
    throw std::invalid_argument("solvationGrid must be a valid xTB Lebedev grid size (e.g. 110, 230, 302, 434)");
  return grid;
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

static bool runMolImpl(const ROMol &mol, XTBResults &result, const XTBOptions &options,
                       xtb_TResults restart = nullptr, xtb_TResults *saveRestart = nullptr) {
  result = XTBResults{};
  const int method = methodNumber(options.method);
  if (options.espPointsAngstrom.size() % 3 || options.espPointsAngstrom.size()/3 >
      static_cast<size_t>(std::numeric_limits<int>::max()) ||
      !std::all_of(options.espPointsAngstrom.begin(), options.espPointsAngstrom.end(),
                   [](double v) { return std::isfinite(v); }))
    throw std::invalid_argument("ESP points must be finite [npoints,3] Angstrom coordinates");
  if (options.realSpacePointsAngstrom.size() % 3 || options.realSpacePointsAngstrom.size()/3 >
      static_cast<size_t>(std::numeric_limits<int>::max()) ||
      !std::all_of(options.realSpacePointsAngstrom.begin(), options.realSpacePointsAngstrom.end(),
                   [](double v) { return std::isfinite(v); }))
    throw std::invalid_argument("Real-space points must be finite [npoints,3] Angstrom coordinates");
  if (!std::isfinite(options.realSpaceOccupationThreshold) || options.realSpaceOccupationThreshold < 0.0 ||
      !std::isfinite(options.realSpaceLEAWindowEV) || options.realSpaceLEAWindowEV < 0.0)
    throw std::invalid_argument("Real-space occupation threshold and LEA window must be finite and nonnegative");
  if (!mol.getNumAtoms()) throw std::invalid_argument("Molecule must contain atoms");
  if (mol.getNumAtoms() > static_cast<unsigned int>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("Too many atoms for the xTB C API");
  if (!mol.getNumConformers()) throw std::invalid_argument("Molecule must have a 3D conformer");
  if (options.confId < -1) throw std::invalid_argument("confId must be -1 or a conformer ID");
  if (options.maxIterations < 1 || !std::isfinite(options.accuracy) || (options.accuracy < 1e-4 || options.accuracy > 1e3) ||
      !std::isfinite(options.electronicTemperature) || options.electronicTemperature < 0)
    throw std::invalid_argument("Invalid iteration limit, accuracy, or electronic temperature");
  const std::string solvModel = normalizeSolvationModel(options.solvationModel);
  const int solvGrid = validateSolvationGrid(options.solvationGrid);
  const std::string solventTrimmed = trimString(options.solvent);
  if (!options.solvent.empty() && solventTrimmed.empty())
    throw std::invalid_argument("solvent must be a nonempty solvent name or empty for gas phase");
  if (!solventTrimmed.empty()) {
    result.solvent = solventTrimmed;
    result.solvationModel = solvModel;
    result.solvationGrid = solvGrid;
  }
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
  result.maxIterations = options.maxIterations;
  result.accuracy = options.accuracy;
  result.electronicTemperature = options.electronicTemperature;
  result.effectiveElectronicTemperature = std::max(options.electronicTemperature, 1e-6);
  result.initialization = restart ? "reference" : "independent";
  result.method = "GFN" + std::to_string(method) + "-xTB";
  result.chargeModel = method == 0 ? "gfn0_eeq" : "xtb_scc_mulliken";
  result.numAtoms = n;
  result.conformerId = conf->getId();
  result.atomicNumbers = numbers;
  result.coordinatesAngstrom.resize(xyz.size());
  std::transform(xyz.begin(), xyz.end(), result.coordinatesAngstrom.begin(),
                 [](double value) { return value / angstromToBohr; });
  auto &api = loadAPI(options.libraryPath);
  result.libraryPath = api.path;
  result.librarySHA256 = api.hash;
  result.apiVersion = api.version;
  result.gfn0StateFixVerified = verifyGFN0(api, result.sourceCommit);
  result.basisExportVerified = verifyBasis(api);
  result.propertiesExportVerified = verifyProperties(api);
  result.realSpaceExportVerified = verifyRealSpace(api);
  if (method == 0 && !result.gfn0StateFixVerified)
    throw std::invalid_argument("Native GFN0 requires a verified patched xTB build; run tools/build_patched_xtb.sh");
  if (!result.basisExportVerified || !verifyDiagnostics(api))
    throw std::invalid_argument(
      "rdXTBTools requires the reviewed extended libxtb build with basis-export and SCC-diagnostics APIs; "
      "run tools/install_native.py or tools/build_patched_xtb.py");
  if (!options.espPointsAngstrom.empty() && !verifyESP(api))
    throw std::invalid_argument("ESP requires the verified ESP extension; rebuild libxtb");
  if (!options.realSpacePointsAngstrom.empty() && !result.realSpaceExportVerified)
    throw std::invalid_argument("Real-space evaluation requires the verified realspace extension; rebuild libxtb");
  result.parameterFile = parameterFile(api, options, method);
  if (!result.parameterFile.empty()) result.parameterSHA256 = sha256(result.parameterFile);
  Handles h(api);
  bool inSinglePoint = false;
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
    if (!result.solvent.empty()) {
      std::string solvArg = result.solvent;
      std::string modelArg = result.solvationModel;
      int gridArg = result.solvationGrid;
      api.setSolvent(h.env, h.calc, solvArg.data(), nullptr, nullptr, &gridArg, modelArg.data());
      api.check(h.env, "Solvation setup");
    }
    api.setMaxIter(h.env, h.calc, options.maxIterations);
    api.setAccuracy(h.env, h.calc, options.accuracy);
    api.setElectronicTemp(h.env, h.calc, options.electronicTemperature);
    api.check(h.env, "Calculator settings");
    h.res = restart ? api.copyResults(restart) : api.newResults();
    if (!h.res) throw EngineError(restart ? "Cannot copy reference wavefunction" : "Cannot allocate xTB results");
    if (options.keepOverlap) {
      api.setKeepOverlap(h.env, h.res, true);
      api.check(h.env, "Basis retention setup");
    }
    if (!options.espPointsAngstrom.empty()) {
      api.setKeepESP(h.env, h.res, true);
      api.check(h.env, "ESP retention setup");
    }
    inSinglePoint = true;
    api.singlepoint(h.env, h.mol, h.calc, h.res);
    api.check(h.env, "Single point");
    inSinglePoint = false;
    if (method != 0 && api.getSCFAPIVersion && api.getSCFDiagnostics && result.gfn0StateFixVerified) {
      boost::property_tree::ptree record;
      boost::property_tree::read_json((fs::path(api.path).parent_path().parent_path() /
                                      "share/xtb/eled/build.json").string(), record);
      if (record.get<int>("capabilities.scf_diagnostics_api", 0) == 1 &&
          record.get<std::string>("source.diagnostics_extension.patch_sha256", "") == scfPatchSHA256) {
        if (api.getSCFAPIVersion() != 1) throw EngineError("Unsupported SCC diagnostics API");
        std::array<double, 4> data{};
        if (api.getSCFDiagnostics(h.env, h.res, &result.scfIterations, data.data())) {
          if (result.scfIterations < 1 || !std::all_of(data.begin(), data.end(),
                [](double v){ return std::isfinite(v) && v >= 0.; }) || data[2] <= 0. || data[3] <= 0.)
            throw EngineError("Invalid SCC diagnostics");
          result.scfEnergyChangeHartree = data[0]; result.scfMixerRMS = data[1];
          result.scfEnergyThresholdHartree = data[2]; result.scfMixerThreshold = data[3];
        }
      }
    }
    api.getEnergy(h.env, h.res, &result.totalEnergyHartree);
    api.getNao(h.env, h.res, &result.numOrbitals);
    api.check(h.env, "Energy and basis size");
    if (!std::isfinite(result.totalEnergyHartree) || result.numOrbitals < 1)
      throw EngineError("Invalid total energy or orbital count");
    const size_t nao = result.numOrbitals;
    result.atomicCharges.resize(n);
    result.gradientHartreePerBohr.resize(static_cast<size_t>(3) * n);
    result.dipoleAtomicUnits.resize(3);
    result.wibergBondOrders.resize(static_cast<size_t>(n) * n);
    result.orbitalEnergiesEV.resize(nao);
    result.orbitalOccupancies.resize(nao);
    api.getCharges(h.env, h.res, result.atomicCharges.data());
    api.getGradient(h.env, h.res, result.gradientHartreePerBohr.data());
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
    if (!options.espPointsAngstrom.empty()) {
      const int count = static_cast<int>(options.espPointsAngstrom.size()/3);
      auto points = options.espPointsAngstrom;
      for (auto &v : points) v *= angstromToBohr;
      std::vector<double> values(count);
      const int ok = api.getESP(h.env, h.res, count, points.data(), count, values.data());
      api.check(h.env, "ESP evaluation");
      if (!ok) throw EngineError("ESP evaluation failed without a diagnostic");
      requireFinite(values, "ESP");
      result.espPointsAngstrom = options.espPointsAngstrom;
      result.espValuesAU = std::move(values);
    }
    if (!options.realSpacePointsAngstrom.empty()) {
      const int count = static_cast<int>(options.realSpacePointsAngstrom.size() / 3);
      auto points = options.realSpacePointsAngstrom;
      for (auto &v : points) v *= angstromToBohr;
      std::vector<double> density(count), alie(count), virDensity(count), lea(count);
      const int ok = api.evalRealSpace(
          h.env, h.res, count, points.data(),
          options.realSpaceOccupationThreshold, options.realSpaceLEAWindowEV,
          density.data(), alie.data(), virDensity.data(), lea.data());
      api.check(h.env, "Real-space evaluation");
      if (!ok) throw EngineError("Real-space evaluation failed without a diagnostic");
      requireFinite(density, "real-space density");
      requireFinite(alie, "real-space ALIE");
      requireFinite(virDensity, "real-space virtual density");
      requireFinite(lea, "real-space LEA");
      result.realSpacePointsAngstrom = options.realSpacePointsAngstrom;
      result.realSpaceDensityAU = std::move(density);
      result.realSpaceALIEEV = std::move(alie);
      result.realSpaceVirtualDensityAU = std::move(virDensity);
      result.realSpaceLEAEV = std::move(lea);
    }
    if (result.propertiesExportVerified) {
      int eFlags = 0;
      std::array<double, 15> eVals{};
      if (api.getEnergyDecomposition(h.env, h.res, &eFlags, static_cast<int>(eVals.size()), eVals.data())) {
        result.energyDecompositionFlags = eFlags;
        if (eFlags & 2) result.electronicEnergyHartree = eVals[1];
        if (eFlags & 4) result.repulsionEnergyHartree = eVals[2];
        if (eFlags & 8) result.electrostaticEnergyHartree = eVals[3];
        if (eFlags & 16) result.anisotropicElectrostaticEnergyHartree = eVals[4];
        if (eFlags & 32) result.anisotropicXCEnergyHartree = eVals[5];
        if (eFlags & 64) result.dispersionEnergyHartree = eVals[6];
        if (eFlags & 128) result.halogenBondEnergyHartree = eVals[7];
        if (eFlags & 256) result.solvationEnergyHartree = eVals[8];
        if (eFlags & 512) {
          result.fermiLevelAlphaHartree = eVals[9];
          result.fermiLevelBetaHartree = eVals[10];
          result.fermiLevelAlphaEV = eVals[9] * hartreeToEV;
          result.fermiLevelBetaEV = eVals[10] * hartreeToEV;
        }
        if (eFlags & 1024) {
          result.solvationBornEnergyHartree = eVals[11];
          result.solvationSASAEnergyHartree = eVals[12];
          result.solvationHBEnergyHartree = eVals[13];
          result.solvationShiftEnergyHartree = eVals[14];
        }
      }
      std::array<double, 3> molPol{};
      std::vector<double> alphaAtom(n), c6Atom(n);
      if (api.getPolarizability(h.env, h.res, n, molPol.data(), alphaAtom.data(), c6Atom.data())) {
        requireFinite(alphaAtom, "atomic polarizabilities");
        requireFinite(c6Atom, "atomic C6");
        if (!std::isfinite(molPol[0]) || !std::isfinite(molPol[1]) || !std::isfinite(molPol[2]))
          throw EngineError("Nonfinite molecular polarizability");
        result.polarizabilityAvailable = true;
        result.molecularPolarizabilityAU = molPol[0];
        result.molecularC6AU = molPol[1];
        result.molecularC8AU = molPol[2];
        result.atomicPolarizabilitiesAU = std::move(alphaAtom);
        result.atomicC6AU = std::move(c6Atom);
      }
      std::vector<double> dipAtom(static_cast<size_t>(3) * n);
      std::vector<double> qpAtom(static_cast<size_t>(6) * n);
      std::vector<double> molQp(6);
      if (api.getMultipoles(h.env, h.res, n, dipAtom.data(), qpAtom.data(), molQp.data())) {
        requireFinite(dipAtom, "atomic dipoles");
        requireFinite(qpAtom, "atomic quadrupoles");
        requireFinite(molQp, "molecular quadrupole");
        result.multipolesAvailable = true;
        result.atomicDipolesAU = std::move(dipAtom);
        result.atomicQuadrupolesAU = std::move(qpAtom);
        result.molecularQuadrupoleAU = std::move(molQp);
      }
      std::vector<double> nuclearESP(n);
      const int epnOk = api.getNuclearESP(h.env, h.res, n, nuclearESP.data());
      api.check(h.env, "Nuclear ESP evaluation");
      if (!epnOk) throw EngineError("Nuclear ESP evaluation failed without a diagnostic");
      requireFinite(nuclearESP, "nuclear ESP");
      result.nuclearESPAvailable = true;
      result.nuclearESPAU = std::move(nuclearESP);
    }
    api.check(h.env, "Result extraction");
    requireFinite(result.gradientHartreePerBohr, "gradient");
    for (auto &energy : result.orbitalEnergiesEV) energy *= hartreeToEV;
    // The Fortran WBO array is column-major; expose atom rows and columns.
    for (size_t i = 0; i < static_cast<size_t>(n); ++i)
      for (size_t j = 0; j < i; ++j)
        std::swap(result.wibergBondOrders[i * n + j], result.wibergBondOrders[j * n + i]);
    requireFinite(result.atomicCharges, "atomic charges");
    if (result.chargeModel == "xtb_scc_mulliken") {
      result.cm5Charges = calcCM5Charges(result.atomicNumbers, result.coordinatesAngstrom,
                                         result.atomicCharges);
    }
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
    if (saveRestart) {
      *saveRestart = h.res;
      h.res = nullptr;
    }
    return true;
  } catch (const EngineError &error) {
    result.errorMessage = error.what();
    result.failureKind = inSinglePoint ? "calculation" : "environment";
    return false;
  }
}

bool runMol(const ROMol &mol, XTBResults &result, const XTBOptions &options) {
  std::lock_guard<std::mutex> lock(calculationMutex);
  return runMolImpl(mol, result, options);
}

std::vector<XTBResults> runChargeStates(const ROMol &mol, const XTBOptions &options,
    int ionizedUnpaired, int attachedUnpaired, bool reuseReference) {
  std::lock_guard<std::mutex> lock(calculationMutex);
  if (methodNumber(options.method) == 0)
    throw std::invalid_argument("Charge-state SCF requires GFN1/2");
  auto &api = loadAPI(options.libraryPath);
  Handles reference(api);
  std::vector<XTBResults> results;
  for (int index = 0; index < 3; ++index) {
    auto selected = options;
    if (index) {
      const long long q = static_cast<long long>(results.front().charge) + (index == 1 ? 1 : -1);
      if (q < std::numeric_limits<int>::min() || q > std::numeric_limits<int>::max())
        throw std::invalid_argument("Adjacent charge is out of range");
      selected.charge = static_cast<int>(q);
      selected.numUnpairedElectrons = index == 1 ? ionizedUnpaired : attachedUnpaired;
    }
    results.emplace_back();
    if (!runMolImpl(mol, results.back(), selected,
        index && reuseReference ? reference.res : nullptr,
        !index && reuseReference ? &reference.res : nullptr)) break;
  }
  return results;
}

std::vector<double> calcCM5Charges(const std::vector<int> &atomicNumbers,
                                   const std::vector<double> &coordinatesAngstrom,
                                   const std::vector<double> &mullikenCharges) {
  static constexpr std::array<double, 118> kRadiiAngstrom = {
      0.32, 0.37, 1.30, 0.99, 0.84, 0.75, 0.71, 0.64, 0.60, 0.62,
      1.60, 1.40, 1.24, 1.14, 1.09, 1.04, 1.00, 1.01, 2.00, 1.74,
      1.59, 1.48, 1.44, 1.30, 1.29, 1.24, 1.18, 1.17, 1.22, 1.20,
      1.23, 1.20, 1.20, 1.18, 1.17, 1.16, 2.15, 1.90, 1.76, 1.64,
      1.56, 1.46, 1.38, 1.36, 1.34, 1.30, 1.36, 1.40, 1.42, 1.40,
      1.40, 1.37, 1.36, 1.36, 2.38, 2.06, 1.94, 1.84, 1.90, 1.88,
      1.86, 1.85, 1.83, 1.82, 1.81, 1.80, 1.79, 1.77, 1.77, 1.78,
      1.74, 1.64, 1.58, 1.50, 1.41, 1.36, 1.32, 1.30, 1.30, 1.32,
      1.44, 1.45, 1.50, 1.42, 1.48, 1.46, 2.42, 2.11, 2.01, 1.90,
      1.84, 1.83, 1.80, 1.80, 1.73, 1.68, 1.68, 1.68, 1.65, 1.67,
      1.73, 1.76, 1.61, 1.57, 1.49, 1.43, 1.41, 1.34, 1.29, 1.28,
      1.21, 1.22, 1.36, 1.43, 1.62, 1.75, 1.65, 1.57,
  };
  static constexpr std::array<double, 118> kA0 = {
      0.0056,  -0.1543, 0.0000,  0.0333,  -0.1030, -0.0446, -0.1072,
      -0.0802, -0.0629, -0.1088, 0.0184,  0.0000,  -0.0726, -0.0790,
      -0.0756, -0.0565, -0.0444, -0.0767, 0.0130,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  -0.0512, -0.0557, -0.0533, -0.0399, -0.0313,
      -0.0541, 0.0092,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  -0.0361,
      -0.0393, -0.0376, -0.0281, -0.0220, -0.0381, 0.0065,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  -0.0255, -0.0277, -0.0265, -0.0198, -0.0155,
      -0.0269, 0.0046,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,  0.0000,
      -0.0179, -0.0195, -0.0187, -0.0140, -0.0110, -0.0189,
  };
  const size_t n = atomicNumbers.size();
  if (coordinatesAngstrom.size() != 3 * n || mullikenCharges.size() != n)
    throw std::invalid_argument("CM5 inputs must have matching atom counts");
  for (int z : atomicNumbers)
    if (z < 1 || z > 118) throw std::invalid_argument("CM5 supports atomic numbers 1..118");
  for (double v : coordinatesAngstrom)
    if (!std::isfinite(v)) throw std::invalid_argument("CM5 coordinates must be finite 3-vectors");
  for (double v : mullikenCharges)
    if (!std::isfinite(v)) throw std::invalid_argument("Mulliken charges must be finite");
  std::vector<double> out(mullikenCharges);
  for (size_t i = 0; i < n; ++i) {
    const int zi = atomicNumbers[i];
    for (size_t j = 0; j < i; ++j) {
      const int zj = atomicNumbers[j];
      const double dx = coordinatesAngstrom[3 * i] - coordinatesAngstrom[3 * j];
      const double dy = coordinatesAngstrom[3 * i + 1] - coordinatesAngstrom[3 * j + 1];
      const double dz = coordinatesAngstrom[3 * i + 2] - coordinatesAngstrom[3 * j + 2];
      const double distance = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (distance == 0.0) throw std::invalid_argument("CM5 requires distinct atom coordinates");
      if (zi == zj) continue;
      const int zMin = std::min(zi, zj);
      const int zMax = std::max(zi, zj);
      double base = 0.0;
      bool paired = true;
      if (zMin == 1 && zMax == 6) base = 0.0502;
      else if (zMin == 1 && zMax == 7) base = 0.1747;
      else if (zMin == 1 && zMax == 8) base = 0.1671;
      else if (zMin == 6 && zMax == 7) base = 0.0556;
      else if (zMin == 6 && zMax == 8) base = 0.0234;
      else if (zMin == 7 && zMax == 8) base = -0.0346;
      else paired = false;
      const double pij = paired ? (zi < zj ? base : -base) : (kA0[zi - 1] - kA0[zj - 1]);
      const double weight = std::exp(-2.4740 * (distance - kRadiiAngstrom[zi - 1] - kRadiiAngstrom[zj - 1]));
      out[i] += weight * pij;
      out[j] -= weight * pij;
    }
  }
  return out;
}
}  // namespace RDKit::XTBTools
