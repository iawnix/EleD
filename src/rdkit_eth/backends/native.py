"""ctypes binding for the xTB C API, loaded only when requested."""

from __future__ import annotations
import ctypes
import os
import time
from functools import wraps
from threading import RLock
from typing import Any, Sequence

from ..errors import XTBNativeNotFoundError, XTBRunError, XTBEnvironmentError, DescriptorUnavailableError
from ..options import calculation_settings
from ..result import XTBResult, _frontier_energies
from ..input import _xyz_from_mol, _method_flag, _effective_charge_uhf, _input_metadata
from ..runtime import _find_xtb_library, _engine_identity, _parameter_file

HARTREE_TO_EV = 27.211386245988
ANGSTROM_TO_BOHR = 1.8897261254578281

class _NativeXTB:
    """Small ctypes binding for the stable xTB C API.

    The binding is loaded lazily.  Installing rdkit-eth therefore does not
    require a shared xTB library until a native calculation is requested.
    """

    _handle = ctypes.c_void_p

    def __init__(self, library: str):
        self.library = library
        self.lib = ctypes.CDLL(library)
        self.get_api_version = self._bind("xtb_getAPIVersion", ctypes.c_int, [])
        api_version = self.get_api_version()
        # API 2 changes xtb_setSolvent, which this binding does not use. The
        # gas-phase molecule/calculator/result functions bound here retain
        # their API 1 signatures (verified against upstream a9040512).
        if api_version // 10000 not in (1, 2):
            raise XTBNativeNotFoundError(f'unsupported xTB C API version {api_version}')
        self._bind_api()

    def _bind(self, name: str, restype: Any, argtypes: Sequence[Any]) -> Any:
        try:
            function = getattr(self.lib, name)
        except AttributeError as exc:
            raise XTBNativeNotFoundError(f"xTB library is missing required symbol {name}") from exc
        function.restype = restype
        function.argtypes = list(argtypes)
        return function

    def _bind_api(self) -> None:
        h = self._handle
        self.new_environment = self._bind("xtb_newEnvironment", h, [])
        self.del_environment = self._bind("xtb_delEnvironment", None, [ctypes.POINTER(h)])
        self.check_environment = self._bind("xtb_checkEnvironment", ctypes.c_int, [h])
        self.get_error = self._bind(
            "xtb_getError", None, [h, ctypes.c_char_p, ctypes.POINTER(ctypes.c_int)]
        )
        self.set_verbosity = self._bind("xtb_setVerbosity", None, [h, ctypes.c_int])
        self.new_molecule = self._bind(
            "xtb_newMolecule",
            h,
            [
                h,
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_bool),
            ],
        )
        self.del_molecule = self._bind("xtb_delMolecule", None, [ctypes.POINTER(h)])
        self.new_calculator = self._bind("xtb_newCalculator", h, [])
        self.del_calculator = self._bind("xtb_delCalculator", None, [ctypes.POINTER(h)])
        for method in ("GFN0xTB", "GFN1xTB", "GFN2xTB", "GFNFF"):
            setattr(
                self,
                "load_" + method.lower(),
                self._bind("xtb_load" + method, None, [h, h, h, ctypes.c_char_p]),
            )
        self.set_max_iterations = self._bind('xtb_setMaxIter', None, [h, h, ctypes.c_int])
        self.set_accuracy = self._bind('xtb_setAccuracy', None, [h, h, ctypes.c_double])
        self.set_temperature = self._bind('xtb_setElectronicTemp', None, [h, h, ctypes.c_double])
        self.new_results = self._bind("xtb_newResults", h, [])
        self.copy_results = self._bind("xtb_copyResults", h, [h])
        self.del_results = self._bind("xtb_delResults", None, [ctypes.POINTER(h)])
        self.singlepoint = self._bind("xtb_singlepoint", None, [h, h, h, h])
        self.get_energy = self._bind("xtb_getEnergy", None, [h, h, ctypes.POINTER(ctypes.c_double)])
        self.get_dipole = self._bind(
            "xtb_getDipole", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )
        self.get_charges = self._bind(
            "xtb_getCharges", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )
        self.get_bond_orders = self._bind(
            "xtb_getBondOrders", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )
        self.get_nao = self._bind("xtb_getNao", None, [h, h, ctypes.POINTER(ctypes.c_int)])
        self.get_orbital_eigenvalues = self._bind(
            "xtb_getOrbitalEigenvalues", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )
        self.get_orbital_occupations = self._bind(
            "xtb_getOrbitalOccupations", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )
        self.get_orbital_coefficients = self._bind(
            "xtb_getOrbitalCoefficients", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )

    def bind_basis_api(self):
        h, i, d = self._handle, ctypes.c_int, ctypes.c_double
        version = self._bind('xtb_rdkit_getBasisAPIVersion', i, [])()
        if version != 1:
            raise DescriptorUnavailableError(f'unsupported native basis API {version}')
        self.set_keep_overlap = self._bind('xtb_rdkit_setKeepOverlap', None, [h, h, ctypes.c_bool])
        self.get_basis_size = self._bind('xtb_rdkit_getBasisSize', None, [h, h, ctypes.POINTER(i)])
        self.get_overlap = self._bind('xtb_rdkit_getOverlap', None, [h, h, i, ctypes.POINTER(d)])
        self.get_basis_map = self._bind('xtb_rdkit_getBasisMap', None,
            [h, h, i, i, ctypes.POINTER(i), ctypes.POINTER(i), ctypes.POINTER(i), ctypes.POINTER(i)])

    def basis_data(self, environment, results, natoms, nao):
        sizes = (ctypes.c_int * 3)()
        self.get_basis_size(environment, results, sizes)
        self.check(environment, 'basis size')
        nshell = sizes[1]
        if sizes[0] != natoms or sizes[2] != nao or not 0 < nshell <= nao or nao * nao > 2147483647:
            raise XTBEnvironmentError('native basis dimensions do not match the result')
        overlap = (ctypes.c_double * (nao * nao))()
        aa, ash = (ctypes.c_int * nao)(), (ctypes.c_int * nao)()
        sa, sl = (ctypes.c_int * nshell)(), (ctypes.c_int * nshell)()
        self.get_overlap(environment, results, nao * nao, overlap)
        self.get_basis_map(environment, results, nao, nshell, aa, ash, sa, sl)
        self.check(environment, 'basis extraction')
        return dict(overlap_matrix=tuple(tuple(overlap[j * nao + i] for j in range(nao)) for i in range(nao)),
                    ao_atom_indices=tuple(aa), ao_shell_indices=tuple(ash),
                    shell_atom_indices=tuple(sa), shell_angular_momenta=tuple(sl))

    def error(self, environment: Any) -> str:
        buffer = ctypes.create_string_buffer(4096)
        size = ctypes.c_int(len(buffer))
        self.get_error(environment, buffer, ctypes.byref(size))
        return buffer.value.decode(errors="replace").strip()

    def check(self, environment: Any, operation: str) -> None:
        status = int(self.check_environment(environment))
        if status:
            detail = self.error(environment)
            error_type = XTBRunError if operation == 'single point' else XTBEnvironmentError
            raise error_type(f"native xTB {operation} failed" + (f": {detail}" if detail else ""))

    @staticmethod
    def _delete(function: Any, handle: Any) -> None:
        if handle:
            if not isinstance(handle, ctypes.c_void_p):
                handle = ctypes.c_void_p(handle)
            function(ctypes.byref(handle))


_native_lock = RLock()


def _serialized(function):
    @wraps(function)
    def call(*args, **kwargs):
        with _native_lock:
            return function(*args, **kwargs)
    return call


@_serialized
def _run_xtb_native(
    mol: Any,
    *,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    uhf: int | None = None,
    xtb_library: str | os.PathLike[str] | None = None,
    keep_overlap: bool = False,
    max_iterations: int = 250,
    accuracy: float = 1.0,
    electronic_temperature: float = 300.0,
    _restart=None,
    _save_restart=None,
) -> XTBResult:
    """Run a single point through a shared ``libxtb`` C API library.

    This adapter supports gas-phase single points only; geometry optimization
    and solvent setup are outside its current API.
    """

    from rdkit import Chem
    settings = calculation_settings(max_iterations, accuracy, electronic_temperature)
    mol = Chem.Mol(mol)  # Own the exact geometry/topology while ctypes releases the GIL.
    if not isinstance(keep_overlap, bool):
        raise ValueError('keep_overlap must be a bool')
    xyz = _xyz_from_mol(mol, conf_id)
    method_name, _ = _method_flag(method)
    if method_name == 'GFN-FF':
        raise ValueError('GFN-FF does not provide molecular orbital descriptors')
    charge, uhf = _effective_charge_uhf(mol, charge, uhf)
    library = _NativeXTB(_find_xtb_library(xtb_library))
    engine = _engine_identity('native', library.library)
    if keep_overlap:
        if not engine['orbital_basis']:
            raise DescriptorUnavailableError('keep_overlap requires a verified native basis extension; rebuild with tools/build_patched_xtb.py')
        library.bind_basis_api()
    # Check again here because callers can use run_xtb_native directly. Do
    # this before creating handles or loading any method's parameters.
    if method_name == 'GFN0-xTB' and not engine['native_gfn0']:
        raise ValueError('native GFN0 requires a verified patched xTB build; '
                         'run tools/build_patched_xtb.sh')
    parameter = _parameter_file(method_name, library.library)
    started = time.perf_counter()
    natoms = mol.GetNumAtoms()
    conf = mol.GetConformer(conf_id)
    atom_count = ctypes.c_int(natoms)
    atomic_numbers = (ctypes.c_int * natoms)(*(int(atom.GetAtomicNum()) for atom in mol.GetAtoms()))
    coordinates = (ctypes.c_double * (3 * natoms))()
    for atom in mol.GetAtoms():
        position = conf.GetAtomPosition(atom.GetIdx())
        offset = 3 * atom.GetIdx()
        coordinates[offset : offset + 3] = (
            position.x * ANGSTROM_TO_BOHR,
            position.y * ANGSTROM_TO_BOHR,
            position.z * ANGSTROM_TO_BOHR,
        )
    charge_value = ctypes.c_double(charge)
    uhf_value = ctypes.c_int(uhf)
    lattice = (ctypes.c_double * 9)(*([0.0] * 9))
    periodic = (ctypes.c_bool * 3)(False, False, False)
    environment = library.new_environment()
    molecule = None
    calculator = None
    results = None
    try:
        if not environment:
            raise XTBEnvironmentError("native xTB could not create an environment")
        library.set_verbosity(environment, 0)
        molecule = library.new_molecule(
            environment,
            ctypes.byref(atom_count),
            atomic_numbers,
            coordinates,
            ctypes.byref(charge_value),
            ctypes.byref(uhf_value),
            lattice,
            periodic,
        )
        library.check(environment, "molecule setup")
        if not molecule:
            raise XTBEnvironmentError('native xTB could not create a molecule')
        calculator = library.new_calculator()
        if not calculator:
            raise XTBEnvironmentError("native xTB could not create a calculator")
        loader_name = {
            "GFN0-xTB": "load_gfn0xtb",
            "GFN1-xTB": "load_gfn1xtb",
            "GFN2-xTB": "load_gfn2xtb",
        }[method_name]
        parameter_path = os.fsencode(parameter['path']) if parameter else None
        getattr(library, loader_name)(environment, molecule, calculator, parameter_path)
        library.check(environment, "parameter loading")
        library.set_max_iterations(environment, calculator, settings['max_iterations'])
        library.set_accuracy(environment, calculator, settings['accuracy'])
        library.set_temperature(environment, calculator, settings['electronic_temperature'])
        library.check(environment, 'calculator settings')
        results = library.copy_results(_restart) if _restart else library.new_results()
        if not results:
            if _restart:
                raise XTBEnvironmentError('native xTB could not copy the reference wavefunction')
            raise XTBEnvironmentError("native xTB could not create a results object")
        if keep_overlap:
            library.set_keep_overlap(environment, results, True)
            library.check(environment, "basis retention setup")
        library.singlepoint(environment, molecule, calculator, results)
        library.check(environment, "single point")

        energy = ctypes.c_double()
        dipole = (ctypes.c_double * 3)()
        charges = (ctypes.c_double * natoms)()
        bond_orders = (ctypes.c_double * (natoms * natoms))()
        nao = ctypes.c_int()
        library.get_energy(environment, results, ctypes.byref(energy))
        library.get_dipole(environment, results, dipole)
        library.get_charges(environment, results, charges)
        library.get_bond_orders(environment, results, bond_orders)
        library.get_nao(environment, results, ctypes.byref(nao))
        library.check(environment, "result extraction")
        if nao.value <= 0:
            raise XTBEnvironmentError("native xTB returned no orbital basis functions")
        orbital_energies = (ctypes.c_double * nao.value)()
        occupations = (ctypes.c_double * nao.value)()
        coefficients = (ctypes.c_double * (nao.value * nao.value))()
        library.get_orbital_eigenvalues(environment, results, orbital_energies)
        library.get_orbital_occupations(environment, results, occupations)
        library.get_orbital_coefficients(environment, results, coefficients)
        library.check(environment, "orbital extraction")
        energies_ev = tuple(value * HARTREE_TO_EV for value in orbital_energies)
        occupation_values = tuple(float(value) for value in occupations)
        homo, lumo = _frontier_energies(energies_ev, occupation_values)
        basis = library.basis_data(environment, results, natoms, nao.value) if keep_overlap else {}
        if _save_restart is not None:
            saved = library.copy_results(results)
            if not saved:
                raise XTBEnvironmentError('native xTB could not copy the reference wavefunction')
            _save_restart.append((library, saved))
        return XTBResult(
            **basis,
            method=method_name,
            version=engine["version"],
            total_energy_eh=float(energy.value),
            homo_lumo_gap_ev=(lumo - homo if homo is not None and lumo is not None else None),
            number_of_molecular_orbitals=int(nao.value),
            number_of_electrons=int(round(sum(occupation_values))),
            number_of_unpaired_electrons=int(uhf),
            orbital_energies_ev=energies_ev,
            orbital_occupations=occupation_values,
            partial_charges=tuple(float(value) for value in charges),
            dipole_au=tuple(float(value) for value in dipole),
            wiberg_bond_orders=tuple(
                tuple(float(bond_orders[row * natoms + column]) for column in range(natoms))
                for row in range(natoms)
            ),
            orbital_coefficients=tuple(
                tuple(float(coefficients[column * nao.value + row]) for column in range(nao.value))
                for row in range(nao.value)
            ),
            backend="native",
            raw_json={"native_api_version": int(library.get_api_version())},
            metadata={
                **_input_metadata(mol, conf_id, charge, uhf, xyz),
                'calculation_settings': settings,
                'calculation_status': 'success',
                'initialization': 'reference' if _restart else 'independent',
                'optimize': False,
                'solvent': None,
                'engine': engine,
                'parameter_file': parameter,
                'native_api_version': int(library.get_api_version()),
                'elapsed_seconds': time.perf_counter() - started,
                'coefficient_axes': ['AO', 'MO'],
                'overlap_source': 'same libxtb calculation' if keep_overlap else None,
            },
        )
    finally:
        library._delete(library.del_results, results)
        library._delete(library.del_calculator, calculator)
        library._delete(library.del_molecule, molecule)
        library._delete(library.del_environment, environment)


def run_xtb_native(mol, *, conf_id=-1, method='GFN2-xTB', charge=None, uhf=None,
                   xtb_library=None, keep_overlap=False, max_iterations=250,
                   accuracy=1.0, electronic_temperature=300.0):
    """Run an independently initialized gas-phase libxtb single point."""
    return _run_xtb_native(mol, conf_id=conf_id, method=method, charge=charge, uhf=uhf,
        xtb_library=xtb_library, keep_overlap=keep_overlap, max_iterations=max_iterations,
        accuracy=accuracy, electronic_temperature=electronic_temperature)


def _native_charge_states(mol, configurations):
    """Scoped reference restart; callers must close this iterator on early exit.

    The fixed-geometry/settings configurations are validated by run_reactivity.
    Each neighbor receives its own deep copy of the untouched reference.
    """
    saved = []
    with _native_lock:
        try:
            for index, config in enumerate(configurations):
                options = config.to_dict()
                options.pop('backend')
                yield _run_xtb_native(mol, **options,
                    _restart=saved[0][1] if index else None,
                    _save_restart=saved if not index else None)
        finally:
            for library, handle in saved:
                library._delete(library.del_results, handle)
