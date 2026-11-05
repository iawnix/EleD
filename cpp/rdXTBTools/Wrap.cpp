#define PY_ARRAY_UNIQUE_SYMBOL rdxtb_array_API
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
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
  std::memcpy(PyArray_DATA(reinterpret_cast<PyArrayObject *>(object)), values.data(),
              values.size() * sizeof(double));
  return py::object(py::handle<>(object));
}

py::object charges(const XTBResults &r) {
  requireSuccess(r); return array(r.atomicCharges, {r.numAtoms});
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

struct ReleaseGIL {
  PyThreadState *state = PyEval_SaveThread();
  ~ReleaseGIL() { PyEval_RestoreThread(state); }
};

py::tuple run(const RDKit::ROMol &mol, int confId, const std::string &method,
              const py::object &charge, const py::object &unpaired, bool keepCoefficients,
              const std::string &library, const std::string &parameterDirectory,
              const py::object &maxIterations, const py::object &accuracy,
              const py::object &electronicTemperature, bool keepOverlap) {
  XTBOptions options;
  options.confId = confId;
  options.method = method;
  options.charge = optionalInteger(charge, "charge");
  options.numUnpairedElectrons = optionalInteger(unpaired, "numUnpairedElectrons");
  options.keepOrbitalCoefficients = keepCoefficients;
  options.keepOverlap = keepOverlap;
  // Use the active Python runtime, including relocated distributions. An explicit
  // path remains authoritative; the standalone C++ default is unchanged.
  options.libraryPath = library.empty() ? py::extract<std::string>(
      py::import("rdkit_eth.runtime").attr("_find_xtb_library")()) : library;
  options.parameterDirectory = parameterDirectory;
  py::object settings = py::import("rdkit_eth.options").attr("calculation_settings")(
      maxIterations, accuracy, electronicTemperature);
  options.maxIterations = py::extract<int>(settings["max_iterations"]);
  options.accuracy = py::extract<double>(settings["accuracy"]);
  options.electronicTemperature = py::extract<double>(settings["electronic_temperature"]);
  // Snapshot the molecule under the GIL before allowing Python threads to run.
  RDKit::ROMol snapshot(mol);
  auto result = std::make_shared<XTBResults>();
  // Coordinate with the existing ctypes adapter as well as the core C++ mutex.
  py::object lock = py::import("rdkit_eth.backends.native").attr("_native_lock");
  lock.attr("acquire")();
  bool ok = false;
  try {
    ReleaseGIL released;
    ok = RDKit::XTBTools::runMol(snapshot, *result, options);
  } catch (...) {
    lock.attr("release")();
    throw;
  }
  lock.attr("release")();
  if (ok) {
    // Hash the immutable snapshot actually passed to C++, not the caller's
    // molecule (which another Python thread could change while GIL is released).
    py::object input = py::import("rdkit_eth.input");
    py::object snapshotObject(py::ptr(&snapshot));
    py::object xyz = input.attr("_xyz_from_mol")(snapshotObject, result->conformerId);
    py::object metadata = input.attr("_input_metadata")(snapshotObject, result->conformerId,
                                  result->charge, result->numUnpairedElectrons, xyz);
    result->inputProvenanceJSON = py::extract<std::string>(py::import("json").attr("dumps")(metadata));
  }
  return py::make_tuple(ok, result);
}

void invalidArgument(const std::invalid_argument &error) { PyErr_SetString(PyExc_ValueError, error.what()); }
void runtimeError(const std::runtime_error &error) { PyErr_SetString(PyExc_RuntimeError, error.what()); }

// Internal bridge for the validated high-level reactivity protocol. Engine
// handles remain scoped inside the C++ call; Python only receives owned data.
py::list runStates(const RDKit::ROMol &mol, const py::object &config,
                   const py::object &ionized, const py::object &attached, bool reuseReference) {
  XTBOptions options;
  options.confId = py::extract<int>(config.attr("conf_id"));
  options.method = py::extract<std::string>(config.attr("method"));
  options.charge = optionalInteger(config.attr("charge"), "charge");
  options.numUnpairedElectrons = optionalInteger(config.attr("uhf"), "uhf");
  options.libraryPath = py::extract<std::string>(config.attr("xtb_library"));
  options.keepOverlap = py::extract<bool>(config.attr("keep_overlap"));
  py::object settings = py::import("rdkit_eth.options").attr("calculation_settings")(
      py::object(config.attr("max_iterations")), py::object(config.attr("accuracy")),
      py::object(config.attr("electronic_temperature")));
  options.maxIterations = py::extract<int>(settings["max_iterations"]);
  options.accuracy = py::extract<double>(settings["accuracy"]);
  options.electronicTemperature = py::extract<double>(settings["electronic_temperature"]);
  auto ionSpin = optionalInteger(ionized, "ionized_uhf");
  auto attachSpin = optionalInteger(attached, "attached_uhf");
  if (!ionSpin || !attachSpin) throw std::invalid_argument("Adjacent spins must be explicit");
  RDKit::ROMol snapshot(mol);
  std::vector<XTBResults> results;
  py::object lock = py::import("rdkit_eth.backends.native").attr("_native_lock");
  lock.attr("acquire")();
  try {
    ReleaseGIL released;
    results = RDKit::XTBTools::runChargeStates(snapshot, options, *ionSpin, *attachSpin, reuseReference);
  } catch (...) {
    lock.attr("release")();
    throw;
  }
  lock.attr("release")();
  py::list output;
  for (auto &value : results) {
    if (value.success) {
      py::object input = py::import("rdkit_eth.input");
      py::object snapshotObject(py::ptr(&snapshot));
      py::object xyz = input.attr("_xyz_from_mol")(snapshotObject, value.conformerId);
      py::object metadata = input.attr("_input_metadata")(snapshotObject, value.conformerId,
                                       value.charge, value.numUnpairedElectrons, xyz);
      value.inputProvenanceJSON = py::extract<std::string>(py::import("json").attr("dumps")(metadata));
    }
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
      "GFN0 requires a verified patched library. Independent of YAeHMOP.";
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
      .def("GetAtomicCharges", charges, "Copy of atomic partial charges in elementary-charge units, [atom].")
      .def("GetOrbitalEnergies", energies, "Copy of molecular orbital energies in eV, [MO].")
      .def("GetOrbitalOccupancies", occupancies, "Copy of orbital occupancies in electrons, [MO].")
      .def("GetDipole", dipole, "Copy of dipole in e*Bohr, [x,y,z].")
      .def("GetWibergBondOrders", bondOrders, "Copy of dimensionless Wiberg bond orders, [atom,atom].")
      .def("GetOrbitalCoefficients", coefficients, "Copy of optional MO coefficients, [AO,MO].");

  py::def("_RunChargeStates", runStates);
  py::def("RunMol", run,
          (py::arg("mol"), py::arg("confId") = -1, py::arg("method") = "GFN2-xTB",
           py::arg("charge") = py::object(), py::arg("numUnpairedElectrons") = py::object(),
           py::arg("keepOrbitalCoefficients") = false, py::arg("xtbLibrary") = "",
           py::arg("parameterDirectory") = "", py::arg("maxIterations") = 250,
           py::arg("accuracy") = 1.0, py::arg("electronicTemperature") = 300.0, py::arg("keepOverlap") = false),
          "Run a gas-phase xTB single point on an explicit-H 3D RDKit molecule. "
          "Returns (success, XTBResults). xTB failures return False and errorMessage; "
          "numerical getters then raise RuntimeError. Invalid input raises ValueError; "
          "library/ABI errors raise RuntimeError. Coordinates use Angstrom. "
          "Charge and unpaired electrons default to formal charges and radicals. "
          "Calls are serialized and do not modify mol. The GIL is released during C++ calculation.");
}
