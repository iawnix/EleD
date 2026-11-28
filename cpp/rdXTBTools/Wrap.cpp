#define PY_ARRAY_UNIQUE_SYMBOL rdxtb_array_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#define BOOST_PYTHON_MAX_ARITY 24
#include <boost/python.hpp>
#include <numpy/arrayobject.h>
#include <RDGeneral/versions.h>
#include "XTBTools.h"

#include <cstring>
#include <limits>
#include <memory>

namespace py = boost::python;
using RDKit::XTBTools::XTBOptions;
using RDKit::XTBTools::XTBResults;

namespace {
void requireSuccess(const XTBResults &result) {
  if (!result.success) throw std::runtime_error("No numerical results: " + result.errorMessage);
}

py::object array(const std::vector<double> &values, std::vector<npy_intp> shape) {
  PyObject *object = PyArray_SimpleNew(shape.size(), shape.data(), NPY_DOUBLE);
  if (!object) py::throw_error_already_set();
  if (!values.empty()) {
    std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject *>(object)), values.data(),
                values.size() * sizeof(double));
  }
  return py::object(py::handle<>(object));
}

py::object charges(const XTBResults &r) {
  requireSuccess(r); return array(r.atomicCharges, {r.numAtoms});
}
py::object espPoints(const XTBResults &r) {
  requireSuccess(r);
  if (r.espValuesAU.empty()) throw std::invalid_argument("ESP was not sampled; pass espPointsAngstrom");
  return array(r.espPointsAngstrom, {static_cast<npy_intp>(r.espValuesAU.size()), 3});
}
py::object espValues(const XTBResults &r) {
  requireSuccess(r);
  if (r.espValuesAU.empty()) throw std::invalid_argument("ESP was not sampled; pass espPointsAngstrom");
  return array(r.espValuesAU, {static_cast<npy_intp>(r.espValuesAU.size())});
}
std::vector<double> readESPPoints(const py::object &value) {
  if (value.is_none()) return {};
  PyObject *rawObj = PyArray_FromAny(value.ptr(), nullptr, 0, 0, NPY_ARRAY_ENSUREARRAY, nullptr);
  if (!rawObj) {
    PyErr_Clear();
    throw std::invalid_argument("espPointsAngstrom must be a real [npoints, 3] array");
  }
  py::handle<> rawHandle(rawObj);
  auto *raw = reinterpret_cast<PyArrayObject *>(rawObj);
  const int ndim = PyArray_NDIM(raw);
  const npy_intp *dims = PyArray_DIMS(raw);
  if (PyArray_SIZE(raw) == 0 &&
      ((ndim == 1 && dims[0] == 0) || (ndim == 2 && dims[0] == 0 && dims[1] == 3))) {
    return {};
  }
  const char kind = PyArray_DESCR(raw)->kind;
  if ((kind != 'i' && kind != 'u' && kind != 'f') || ndim != 2 || dims[1] != 3)
    throw std::invalid_argument("espPointsAngstrom must be a real [npoints, 3] array");
  if (dims[0] > static_cast<npy_intp>(std::numeric_limits<int>::max()))
    throw std::invalid_argument("espPointsAngstrom must be finite and fit the C API point count");
  PyObject *contigObj = PyArray_FROM_OTF(rawObj, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY);
  if (!contigObj) py::throw_error_already_set();
  py::handle<> contigHandle(contigObj);
  const auto *data = static_cast<const double *>(
      PyArray_DATA(reinterpret_cast<PyArrayObject *>(contigObj)));
  const size_t total = static_cast<size_t>(dims[0]) * 3;
  std::vector<double> output(data, data + total);
  for (double v : output) {
    if (!std::isfinite(v))
      throw std::invalid_argument("espPointsAngstrom must be finite and fit the C API point count");
  }
  return output;
}
bool hasESP(const XTBResults &r) { return r.success && !r.espValuesAU.empty(); }
bool hasRealSpace(const XTBResults &r) { return r.success && !r.realSpaceDensityAU.empty(); }
py::object realSpacePoints(const XTBResults &r) {
  requireSuccess(r);
  if (r.realSpaceDensityAU.empty())
    throw std::invalid_argument("Real-space properties were not sampled; pass realSpacePointsAngstrom");
  return array(r.realSpacePointsAngstrom, {static_cast<npy_intp>(r.realSpaceDensityAU.size()), 3});
}
py::object realSpaceDensity(const XTBResults &r) {
  requireSuccess(r);
  if (r.realSpaceDensityAU.empty())
    throw std::invalid_argument("Real-space properties were not sampled; pass realSpacePointsAngstrom");
  return array(r.realSpaceDensityAU, {static_cast<npy_intp>(r.realSpaceDensityAU.size())});
}
py::object realSpaceALIE(const XTBResults &r) {
  requireSuccess(r);
  if (r.realSpaceALIEEV.empty())
    throw std::invalid_argument("Real-space properties were not sampled; pass realSpacePointsAngstrom");
  return array(r.realSpaceALIEEV, {static_cast<npy_intp>(r.realSpaceALIEEV.size())});
}
py::object realSpaceVirtualDensity(const XTBResults &r) {
  requireSuccess(r);
  if (r.realSpaceVirtualDensityAU.empty())
    throw std::invalid_argument("Real-space properties were not sampled; pass realSpacePointsAngstrom");
  return array(r.realSpaceVirtualDensityAU, {static_cast<npy_intp>(r.realSpaceVirtualDensityAU.size())});
}
py::object realSpaceLEA(const XTBResults &r) {
  requireSuccess(r);
  if (r.realSpaceLEAEV.empty())
    throw std::invalid_argument("Real-space properties were not sampled; pass realSpacePointsAngstrom");
  return array(r.realSpaceLEAEV, {static_cast<npy_intp>(r.realSpaceLEAEV.size())});
}
py::object cm5Charges(const XTBResults &r) {
  requireSuccess(r);
  if (r.chargeModel != "xtb_scc_mulliken" || r.cm5Charges.size() != r.numAtoms)
    throw std::invalid_argument("CM5 charges require GFN1/2 SCC Mulliken charges; GFN0 exposes PEEQ/EEQ");
  return array(r.cm5Charges, {r.numAtoms});
}
py::object energies(const XTBResults &r) {
  requireSuccess(r); return array(r.orbitalEnergiesEV, {r.numOrbitals});
}
py::object occupancies(const XTBResults &r) {
  requireSuccess(r); return array(r.orbitalOccupancies, {r.numOrbitals});
}
py::object dipole(const XTBResults &r) {
  requireSuccess(r); return array(r.dipoleAtomicUnits, {3});
}
py::object gradient(const XTBResults &r) {
  requireSuccess(r); return array(r.gradientHartreePerBohr, {r.numAtoms, 3});
}
py::object bondOrders(const XTBResults &r) {
  requireSuccess(r); return array(r.wibergBondOrders, {r.numAtoms, r.numAtoms});
}
py::object coefficients(const XTBResults &r) {
  requireSuccess(r);
  if (r.orbitalCoefficients.empty())
    throw std::invalid_argument("Orbital coefficients were not retained; set keepOrbitalCoefficients=True");
  return array(r.orbitalCoefficients, {r.numOrbitals, r.numOrbitals});
}
py::object overlap(const XTBResults &r) {
  requireSuccess(r);
  if (r.overlapMatrix.empty()) throw std::invalid_argument("Overlap was not retained; set keepOverlap=True");
  return array(r.overlapMatrix, {r.numOrbitals, r.numOrbitals});
}
py::object integerArray(const std::vector<int> &values) {
  npy_intp size = values.size();
  PyObject *object = PyArray_SimpleNew(1, &size, NPY_INT);
  if (!object) py::throw_error_already_set();
  std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject *>(object)), values.data(), values.size()*sizeof(int));
  return py::object(py::handle<>(object));
}
#define MAP_GETTER(name, field) \
py::object name(const XTBResults &r) { \
  requireSuccess(r); \
  if (r.field.empty()) throw std::invalid_argument("Basis maps were not retained; set keepOverlap=True"); \
  return integerArray(r.field); \
}
MAP_GETTER(aoAtoms, aoAtomMap)
MAP_GETTER(aoShells, aoShellMap)
MAP_GETTER(shellAtoms, shellAtomMap)
MAP_GETTER(shellL, shellAngularMomenta)
#undef MAP_GETTER
bool hasOverlap(const XTBResults &r) { return r.success && !r.overlapMatrix.empty(); }
bool hasEnergyDecomp(const XTBResults &r) { return r.success && r.energyDecompositionFlags != 0; }
bool hasPolarizability(const XTBResults &r) { return r.success && r.polarizabilityAvailable; }
bool hasMultipoles(const XTBResults &r) { return r.success && r.multipolesAvailable; }
bool hasNuclearESP(const XTBResults &r) { return r.success && r.nuclearESPAvailable; }

py::object optionalDouble(const XTBResults &r, double value) {
  requireSuccess(r);
  return std::isfinite(value) ? py::object(value) : py::object();
}
py::object getSolvent(const XTBResults &r) {
  return r.solvent.empty() ? py::object() : py::object(r.solvent);
}
py::object getSolvationModel(const XTBResults &r) {
  return r.solvationModel.empty() ? py::object() : py::object(r.solvationModel);
}
py::object getSolvationGrid(const XTBResults &r) {
  return r.solvent.empty() ? py::object() : py::object(r.solvationGrid);
}

#define OPT_DOUBLE_GETTER(name, field) \
py::object name(const XTBResults &r) { return optionalDouble(r, r.field); }
OPT_DOUBLE_GETTER(eElec, electronicEnergyHartree)
OPT_DOUBLE_GETTER(eRep, repulsionEnergyHartree)
OPT_DOUBLE_GETTER(eES, electrostaticEnergyHartree)
OPT_DOUBLE_GETTER(eAES, anisotropicElectrostaticEnergyHartree)
OPT_DOUBLE_GETTER(eAXC, anisotropicXCEnergyHartree)
OPT_DOUBLE_GETTER(eDisp, dispersionEnergyHartree)
OPT_DOUBLE_GETTER(eXB, halogenBondEnergyHartree)
OPT_DOUBLE_GETTER(gSolv, solvationEnergyHartree)
OPT_DOUBLE_GETTER(gBorn, solvationBornEnergyHartree)
OPT_DOUBLE_GETTER(gSASA, solvationSASAEnergyHartree)
OPT_DOUBLE_GETTER(gHB, solvationHBEnergyHartree)
OPT_DOUBLE_GETTER(gShift, solvationShiftEnergyHartree)
OPT_DOUBLE_GETTER(efaHa, fermiLevelAlphaHartree)
OPT_DOUBLE_GETTER(efbHa, fermiLevelBetaHartree)
OPT_DOUBLE_GETTER(efaEV, fermiLevelAlphaEV)
OPT_DOUBLE_GETTER(efbEV, fermiLevelBetaEV)
OPT_DOUBLE_GETTER(molPol, molecularPolarizabilityAU)
OPT_DOUBLE_GETTER(molC6, molecularC6AU)
OPT_DOUBLE_GETTER(molC8, molecularC8AU)
#undef OPT_DOUBLE_GETTER

py::object atomicPolarizabilities(const XTBResults &r) {
  requireSuccess(r);
  if (!r.polarizabilityAvailable || r.atomicPolarizabilitiesAU.size() != r.numAtoms)
    throw std::invalid_argument("D4 polarizabilities are unavailable for this calculation");
  return array(r.atomicPolarizabilitiesAU, {r.numAtoms});
}
py::object atomicC6(const XTBResults &r) {
  requireSuccess(r);
  if (!r.polarizabilityAvailable || r.atomicC6AU.size() != r.numAtoms)
    throw std::invalid_argument("D4 atomic C6 coefficients are unavailable for this calculation");
  return array(r.atomicC6AU, {r.numAtoms});
}
py::object atomicDipoles(const XTBResults &r) {
  requireSuccess(r);
  if (!r.multipolesAvailable || r.atomicDipolesAU.size() != static_cast<size_t>(3) * r.numAtoms)
    throw std::invalid_argument("Atomic dipoles are unavailable for this calculation");
  return array(r.atomicDipolesAU, {r.numAtoms, 3});
}
py::object atomicQuadrupoles(const XTBResults &r) {
  requireSuccess(r);
  if (!r.multipolesAvailable || r.atomicQuadrupolesAU.size() != static_cast<size_t>(6) * r.numAtoms)
    throw std::invalid_argument("Atomic quadrupoles are unavailable for this calculation");
  return array(r.atomicQuadrupolesAU, {r.numAtoms, 6});
}
py::object molecularQuadrupole(const XTBResults &r) {
  requireSuccess(r);
  if (!r.multipolesAvailable || r.molecularQuadrupoleAU.size() != 6)
    throw std::invalid_argument("Molecular quadrupole is unavailable for this calculation");
  return array(r.molecularQuadrupoleAU, {6});
}
py::object nuclearESP(const XTBResults &r) {
  requireSuccess(r);
  if (!r.nuclearESPAvailable || r.nuclearESPAU.size() != r.numAtoms)
    throw std::invalid_argument("Nuclear ESP is unavailable for this calculation");
  return array(r.nuclearESPAU, {r.numAtoms});
}

std::string parseOptionalSolvent(const py::object &value) {
  if (value.is_none()) return {};
  if (!PyUnicode_Check(value.ptr()))
    throw std::invalid_argument("solvent must be a string or None");
  std::string text = py::extract<std::string>(value);
  if (text.find_first_not_of(" \t\r\n") == std::string::npos)
    throw std::invalid_argument("solvent must be a nonempty solvent name or None");
  return text;
}

double totalEnergy(const XTBResults &r) {
  requireSuccess(r); return r.totalEnergyHartree;
}

std::optional<int> optionalInteger(const py::object &value, const char *name) {
  if (value.is_none()) return std::nullopt;
  if (PyBool_Check(value.ptr()) || !PyIndex_Check(value.ptr()))
    throw std::invalid_argument(std::string(name) + " must be an integer or None");
  py::handle<> indexed(PyNumber_Index(value.ptr()));
  if (!indexed) py::throw_error_already_set();
  const long number = PyLong_AsLong(indexed.get());
  if (PyErr_Occurred()) py::throw_error_already_set();
  if (number < std::numeric_limits<int>::min() || number > std::numeric_limits<int>::max())
    throw std::invalid_argument(std::string(name) + " is out of range");
  return static_cast<int>(number);
}

int parseMaxIterations(const py::object &value) {
  if (PyBool_Check(value.ptr()) || !PyIndex_Check(value.ptr()))
    throw std::invalid_argument("maxIterations (iteration limit) must be an integer in [1, 2147483647]");
  py::handle<> indexed(PyNumber_Index(value.ptr()));
  if (!indexed) py::throw_error_already_set();
  const long number = PyLong_AsLong(indexed.get());
  if (PyErr_Occurred()) {
    PyErr_Clear();
    throw std::invalid_argument("maxIterations (iteration limit) must be an integer in [1, 2147483647]");
  }
  if (number < 1 || number > std::numeric_limits<int>::max())
    throw std::invalid_argument("maxIterations (iteration limit) must be an integer in [1, 2147483647]");
  return static_cast<int>(number);
}

double parseReal(const py::object &value, const char *message) {
  if (PyBool_Check(value.ptr()) || PyComplex_Check(value.ptr()) ||
      PyUnicode_Check(value.ptr()) || PyBytes_Check(value.ptr()) || PyByteArray_Check(value.ptr()))
    throw std::invalid_argument(message);
  PyObject *flt = PyNumber_Float(value.ptr());
  if (!flt) {
    PyErr_Clear();
    throw std::invalid_argument(message);
  }
  py::handle<> handle(flt);
  const double result = PyFloat_AsDouble(flt);
  if (PyErr_Occurred() || !std::isfinite(result)) {
    PyErr_Clear();
    throw std::invalid_argument(message);
  }
  return result;
}

void applySettings(XTBOptions &options, const py::object &maxIterations,
                   const py::object &accuracy, const py::object &electronicTemperature) {
  options.maxIterations = parseMaxIterations(maxIterations);
  options.accuracy = parseReal(accuracy, "accuracy must be a finite real number");
  if (options.accuracy < 1e-4 || options.accuracy > 1e3)
    throw std::invalid_argument("accuracy must be in [0.0001, 1000]; avoid silent engine clamping");
  options.electronicTemperature = parseReal(electronicTemperature,
                                            "electronic temperature must be a finite real number");
  if (options.electronicTemperature < 0.0)
    throw std::invalid_argument("electronic temperature must be nonnegative (kelvin)");
}

py::object optionalModule(const char *name) {
  PyObject *mod = PyImport_ImportModule(name);
  if (mod) return py::object(py::handle<>(mod));
  if (PyErr_ExceptionMatches(PyExc_ModuleNotFoundError)) {
    PyErr_Clear();
    return py::object();
  }
  py::throw_error_already_set();
  return py::object();
}

void recordInputProvenance(const RDKit::ROMol &snapshot, XTBResults &result) {
  py::object input = optionalModule("eled.input");
  if (input.is_none()) return;
  py::object snapshotObject(py::ptr(&snapshot));
  py::object xyz = input.attr("_xyz_from_mol")(snapshotObject, result.conformerId);
  py::object metadata = input.attr("_input_metadata")(snapshotObject, result.conformerId,
                                                      result.charge, result.numUnpairedElectrons, xyz);
  result.inputProvenanceJSON = py::extract<std::string>(py::import("json").attr("dumps")(metadata));
}

struct ReleaseGIL {
  PyThreadState *state = PyEval_SaveThread();
  ~ReleaseGIL() { PyEval_RestoreThread(state); }
};

py::tuple run(const RDKit::ROMol &mol, int confId, const std::string &method,
              const py::object &charge, const py::object &unpaired, bool keepCoefficients,
              const std::string &library, const std::string &parameterDirectory,
              const py::object &maxIterations, const py::object &accuracy,
              const py::object &electronicTemperature, bool keepOverlap, const py::object &espPointsAngstrom,
              const py::object &solvent, const std::string &solvationModel, const py::object &solvationGrid,
              const py::object &realSpacePointsAngstrom, const py::object &realSpaceOccupationThreshold,
              const py::object &realSpaceLEAWindowEV) {
  XTBOptions options;
  options.confId = confId;
  options.method = method;
  options.charge = optionalInteger(charge, "charge");
  options.numUnpairedElectrons = optionalInteger(unpaired, "numUnpairedElectrons");
  options.keepOrbitalCoefficients = keepCoefficients;
  options.keepOverlap = keepOverlap;
  options.espPointsAngstrom = readESPPoints(espPointsAngstrom);
  options.realSpacePointsAngstrom = readESPPoints(realSpacePointsAngstrom);
  options.realSpaceOccupationThreshold = parseReal(
      realSpaceOccupationThreshold, "realSpaceOccupationThreshold must be a finite real number");
  options.realSpaceLEAWindowEV = parseReal(
      realSpaceLEAWindowEV, "realSpaceLEAWindowEV must be a finite real number");
  if (options.realSpaceOccupationThreshold < 0.0 || options.realSpaceLEAWindowEV < 0.0)
    throw std::invalid_argument("realSpaceOccupationThreshold and realSpaceLEAWindowEV must be nonnegative");
  options.solvent = parseOptionalSolvent(solvent);
  options.solvationModel = solvationModel;
  auto gridOpt = optionalInteger(solvationGrid, "solvationGrid");
  options.solvationGrid = gridOpt ? *gridOpt : 230;
  // Use the active Python runtime when eled is installed; an explicit path
  // remains authoritative, and standalone C++/Python fallback is preserved.
  if (library.empty()) {
    py::object runtime = optionalModule("eled.runtime");
    options.libraryPath = runtime.is_none()
        ? std::string()
        : std::string(py::extract<std::string>(runtime.attr("_find_xtb_library")()));
  } else {
    options.libraryPath = library;
  }
  options.parameterDirectory = parameterDirectory;
  applySettings(options, maxIterations, accuracy, electronicTemperature);
  // Snapshot the molecule under the GIL before allowing Python threads to run.
  RDKit::ROMol snapshot(mol);
  auto result = std::make_shared<XTBResults>();
  // Coordinate with the existing ctypes adapter when present as well as the core C++ mutex.
  py::object nativeMod = optionalModule("eled.backends.native");
  py::object lock = nativeMod.is_none() ? py::object() : nativeMod.attr("_native_lock");
  if (!lock.is_none()) lock.attr("acquire")();
  bool ok = false;
  try {
    ReleaseGIL released;
    ok = RDKit::XTBTools::runMol(snapshot, *result, options);
  } catch (...) {
    if (!lock.is_none()) lock.attr("release")();
    throw;
  }
  if (!lock.is_none()) lock.attr("release")();
  if (ok) recordInputProvenance(snapshot, *result);
  return py::make_tuple(ok, result);
}

void invalidArgument(const std::invalid_argument &error) { PyErr_SetString(PyExc_ValueError, error.what()); }
void runtimeError(const std::runtime_error &error) { PyErr_SetString(PyExc_RuntimeError, error.what()); }

// Internal bridge for the validated high-level reactivity protocol. Engine
// handles remain scoped inside the C++ call; Python only receives owned data.
py::list runStates(const RDKit::ROMol &mol, const py::object &config,
                   const py::object &ionized, const py::object &attached, bool reuseReference) {
  XTBOptions options;
  options.confId = py::extract<int>(config.attr("confId"));
  options.method = py::extract<std::string>(config.attr("method"));
  options.charge = optionalInteger(config.attr("charge"), "charge");
  options.numUnpairedElectrons = optionalInteger(config.attr("numUnpairedElectrons"), "numUnpairedElectrons");
  py::object libObj = config.attr("xtbLibrary");
  options.libraryPath = libObj.is_none() ? std::string() : std::string(py::extract<std::string>(libObj));
  options.keepOverlap = py::extract<bool>(config.attr("keepOverlap"));
  options.espPointsAngstrom = readESPPoints(config.attr("espPointsAngstrom"));
  if (PyObject_HasAttrString(config.ptr(), "realSpacePointsAngstrom")) {
    options.realSpacePointsAngstrom = readESPPoints(config.attr("realSpacePointsAngstrom"));
    options.realSpaceOccupationThreshold = parseReal(
        py::object(config.attr("realSpaceOccupationThreshold")),
        "realSpaceOccupationThreshold must be a finite real number");
    options.realSpaceLEAWindowEV = parseReal(
        py::object(config.attr("realSpaceLEAWindowEV")),
        "realSpaceLEAWindowEV must be a finite real number");
  }
  if (PyObject_HasAttrString(config.ptr(), "solvent")) {
    options.solvent = parseOptionalSolvent(config.attr("solvent"));
    options.solvationModel = py::extract<std::string>(config.attr("solvationModel"));
    auto gridOpt = optionalInteger(config.attr("solvationGrid"), "solvationGrid");
    options.solvationGrid = gridOpt ? *gridOpt : 230;
  }
  applySettings(options, py::object(config.attr("maxIterations")),
                py::object(config.attr("accuracy")),
                py::object(config.attr("electronicTemperature")));
  auto ionSpin = optionalInteger(ionized, "ionizedUnpaired");
  auto attachSpin = optionalInteger(attached, "attachedUnpaired");
  if (!ionSpin || !attachSpin) throw std::invalid_argument("Adjacent spins must be explicit");
  RDKit::ROMol snapshot(mol);
  std::vector<XTBResults> results;
  py::object nativeMod = optionalModule("eled.backends.native");
  py::object lock = nativeMod.is_none() ? py::object() : nativeMod.attr("_native_lock");
  if (!lock.is_none()) lock.attr("acquire")();
  try {
    ReleaseGIL released;
    results = RDKit::XTBTools::runChargeStates(snapshot, options, *ionSpin, *attachSpin, reuseReference);
  } catch (...) {
    if (!lock.is_none()) lock.attr("release")();
    throw;
  }
  if (!lock.is_none()) lock.attr("release")();
  py::list output;
  for (auto &value : results) {
    if (value.success) recordInputProvenance(snapshot, value);
    output.append(std::make_shared<XTBResults>(std::move(value)));
  }
  return output;
}
}  // namespace

BOOST_PYTHON_MODULE(rdXTBTools) {
  py::import("rdkit.Chem.rdchem");  // Register the ROMol Python converters.
  if (_import_array() < 0) py::throw_error_already_set();
  const std::string installed = py::extract<std::string>(py::import("rdkit.rdBase").attr("rdkitVersion"));
  if (installed != RDXTB_RDKIT_VERSION) {
    PyErr_SetString(PyExc_ImportError, "rdXTBTools was built for a different RDKit version; rebuild the extension");
    py::throw_error_already_set();
  }
  py::register_exception_translator<std::invalid_argument>(invalidArgument);
  py::register_exception_translator<std::runtime_error>(runtimeError);
  py::scope().attr("__doc__") =
      "RDKit ROMol to libxtb C++ interface. Gas-phase GFN0/1/2 single points. "
      "Total energies are Hartree; orbital energies are eV; dipoles are e*Bohr. "
      "The reviewed extended libxtb is required. Independent of YAeHMOP.";
  py::scope().attr("rdkitVersion") = RDXTB_RDKIT_VERSION;
  py::scope().attr("totalEnergyUnit") = "Hartree";
  py::scope().attr("orbitalEnergyUnit") = "eV";
  py::scope().attr("dipoleUnit") = "e*Bohr";
  py::scope().attr("coefficientAxes") = py::make_tuple("AO", "MO");

  py::class_<XTBResults, std::shared_ptr<XTBResults>, boost::noncopyable>("XTBResults", py::no_init)
      .def_readonly("success", &XTBResults::success)
      .def_readonly("errorMessage", &XTBResults::errorMessage)
      .def_readonly("failureKind", &XTBResults::failureKind)
      .def_readonly("method", &XTBResults::method)
      .def_readonly("conformerId", &XTBResults::conformerId)
      .def_readonly("numAtoms", &XTBResults::numAtoms)
      .def_readonly("numOrbitals", &XTBResults::numOrbitals)
      .def_readonly("numElectrons", &XTBResults::numElectrons)
      .def_readonly("numTotalElectrons", &XTBResults::numTotalElectrons)
      .def_readonly("numUnpairedElectrons", &XTBResults::numUnpairedElectrons)
      .def_readonly("charge", &XTBResults::charge)
      .def_readonly("chargeModel", &XTBResults::chargeModel)
      .def_readonly("apiVersion", &XTBResults::apiVersion)
      .def_readonly("maxIterations", &XTBResults::maxIterations)
      .def_readonly("accuracy", &XTBResults::accuracy)
      .def_readonly("electronicTemperature", &XTBResults::electronicTemperature)
      .def_readonly("effectiveElectronicTemperature", &XTBResults::effectiveElectronicTemperature)
      .def_readonly("gfn0StateFixVerified", &XTBResults::gfn0StateFixVerified)
      .def_readonly("scfIterations", &XTBResults::scfIterations)
      .def_readonly("scfEnergyChangeHartree", &XTBResults::scfEnergyChangeHartree)
      .def_readonly("scfMixerRMS", &XTBResults::scfMixerRMS)
      .def_readonly("scfEnergyThresholdHartree", &XTBResults::scfEnergyThresholdHartree)
      .def_readonly("scfMixerThreshold", &XTBResults::scfMixerThreshold)
      .def_readonly("basisExportVerified", &XTBResults::basisExportVerified)
      .def_readonly("inputProvenanceJSON", &XTBResults::inputProvenanceJSON)
      .def_readonly("initialization", &XTBResults::initialization)
      .add_property("overlapRetained", hasOverlap)
      .def("GetOverlapMatrix", overlap)
      .def("GetAOAtomMap", aoAtoms)
      .def("GetAOShellMap", aoShells)
      .def("GetShellAtomMap", shellAtoms)
      .def("GetShellAngularMomenta", shellL)
      .def_readonly("libraryPath", &XTBResults::libraryPath)
      .def_readonly("librarySHA256", &XTBResults::librarySHA256)
      .def_readonly("parameterFile", &XTBResults::parameterFile)
      .def_readonly("parameterSHA256", &XTBResults::parameterSHA256)
      .def_readonly("sourceCommit", &XTBResults::sourceCommit)
      .add_property("totalEnergyHartree", totalEnergy, "xTB total energy in Hartree; not EHT totalEnergy in eV.")
      .add_property("espSampled", hasESP)
      .def("GetESPPoints", espPoints, "Copy of ESP points in Angstrom, [point,xyz].")
      .def("GetESPValues", espValues, "Copy of ESP in Hartree/e, [point].")
      .def_readonly("realSpaceExportVerified", &XTBResults::realSpaceExportVerified)
      .add_property("realSpaceSampled", hasRealSpace)
      .def("GetRealSpacePoints", realSpacePoints, "Copy of real-space points in Angstrom, [point,xyz].")
      .def("GetRealSpaceDensity", realSpaceDensity, "Copy of valence electron density in e/Bohr^3, [point].")
      .def("GetRealSpaceALIE", realSpaceALIE, "Copy of Average Local Ionization Energy in eV, [point].")
      .def("GetRealSpaceVirtualDensity", realSpaceVirtualDensity, "Copy of virtual orbital density in 1/Bohr^3, [point].")
      .def("GetRealSpaceLEA", realSpaceLEA, "Copy of Local Electron Affinity in eV, [point].")
      .def("GetAtomicCharges", charges,
           "Copy of xTB atomic charges in elementary-charge units, [atom]. "
           "GFN1/2: converged SCC Mulliken charges; GFN0: PEEQ/EEQ charges; never CM5.")
      .def("GetCM5Charges", cm5Charges,
           "Return xTB CM5 charges in elementary-charge units, [atom]. "
           "Only defined for GFN1/2 SCC Mulliken charges; GFN0 has no CM5 result.")
      .def("GetOrbitalEnergies", energies, "Copy of molecular orbital energies in eV, [MO].")
      .def("GetOrbitalOccupancies", occupancies, "Copy of orbital occupancies in electrons, [MO].")
      .def("GetDipole", dipole, "Copy of dipole in e*Bohr, [x,y,z].")
      .def("GetGradient", gradient, "Copy of nuclear gradient in Hartree/Bohr, [atom,xyz].")
      .def("GetWibergBondOrders", bondOrders, "Copy of dimensionless Wiberg bond orders, [atom,atom].")
      .def("GetOrbitalCoefficients", coefficients, "Copy of optional MO coefficients, [AO,MO].")
      .add_property("solvent", getSolvent)
      .add_property("solvationModel", getSolvationModel)
      .add_property("solvationGrid", getSolvationGrid)
      .def_readonly("propertiesExportVerified", &XTBResults::propertiesExportVerified)
      .add_property("energyDecompositionAvailable", hasEnergyDecomp)
      .add_property("polarizabilityAvailable", hasPolarizability)
      .add_property("multipolesAvailable", hasMultipoles)
      .add_property("nuclearESPAvailable", hasNuclearESP)
      .add_property("electronicEnergyHartree", eElec)
      .add_property("repulsionEnergyHartree", eRep)
      .add_property("electrostaticEnergyHartree", eES)
      .add_property("anisotropicElectrostaticEnergyHartree", eAES)
      .add_property("anisotropicXCEnergyHartree", eAXC)
      .add_property("dispersionEnergyHartree", eDisp)
      .add_property("halogenBondEnergyHartree", eXB)
      .add_property("solvationEnergyHartree", gSolv)
      .add_property("solvationBornEnergyHartree", gBorn)
      .add_property("solvationSASAEnergyHartree", gSASA)
      .add_property("solvationHBEnergyHartree", gHB)
      .add_property("solvationShiftEnergyHartree", gShift)
      .add_property("fermiLevelAlphaHartree", efaHa)
      .add_property("fermiLevelBetaHartree", efbHa)
      .add_property("fermiLevelAlphaEV", efaEV)
      .add_property("fermiLevelBetaEV", efbEV)
      .def_readonly("energyDecompositionFlags", &XTBResults::energyDecompositionFlags)
      .add_property("molecularPolarizabilityAU", molPol)
      .add_property("molecularC6AU", molC6)
      .add_property("molecularC8AU", molC8)
      .def("GetAtomicPolarizabilities", atomicPolarizabilities,
           "Copy of D4 atomic static dipole polarizabilities in Bohr^3, [atom].")
      .def("GetAtomicC6", atomicC6,
           "Copy of D4 homoatomic dispersion coefficients C6_AA in Hartree*Bohr^6, [atom].")
      .def("GetAtomicDipoles", atomicDipoles,
           "Copy of CAMM atomic dipoles in e*Bohr, [atom,3].")
      .def("GetAtomicQuadrupoles", atomicQuadrupoles,
           "Copy of traceless CAMM atomic quadrupoles in e*Bohr^2 (xx,xy,yy,xz,yz,zz), [atom,6].")
      .def("GetMolecularQuadrupole", molecularQuadrupole,
           "Copy of traceless molecular quadrupole in e*Bohr^2 (xx,xy,yy,xz,yz,zz), [6].")
      .def("GetNuclearESP", nuclearESP,
           "Copy of exact electrostatic potential at nuclei (EPN, V_N) in Hartree/e, [atom].");

  py::def("_RunChargeStates", runStates);
  py::def("RunMol", run,
          (py::arg("mol"), py::arg("confId") = -1, py::arg("method") = "GFN2-xTB",
           py::arg("charge") = py::object(), py::arg("numUnpairedElectrons") = py::object(),
           py::arg("keepOrbitalCoefficients") = false, py::arg("xtbLibrary") = "",
           py::arg("parameterDirectory") = "", py::arg("maxIterations") = 250,
           py::arg("accuracy") = 1.0, py::arg("electronicTemperature") = 300.0, py::arg("keepOverlap") = false,
           py::arg("espPointsAngstrom") = py::object(), py::arg("solvent") = py::object(),
           py::arg("solvationModel") = "alpb", py::arg("solvationGrid") = py::object(),
           py::arg("realSpacePointsAngstrom") = py::object(),
           py::arg("realSpaceOccupationThreshold") = 1e-7,
           py::arg("realSpaceLEAWindowEV") = 0.0),
          "Run a gas-phase xTB single point on an explicit-H 3D RDKit molecule. "
          "Returns (success, XTBResults). xTB failures return False and errorMessage; "
          "numerical getters then raise RuntimeError. Invalid input raises ValueError; "
          "library/ABI errors raise RuntimeError. Coordinates use Angstrom. "
          "Charge and unpaired electrons default to formal charges and radicals. "
          "Calls are serialized and do not modify mol. The GIL is released during C++ calculation.");
}
