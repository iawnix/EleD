"""ctypes binding for the xTB C API, loaded only when requested."""

from __future__ import annotations

import ctypes
from functools import wraps
import os
from threading import RLock
import time
from typing import Any, Sequence

from ..charges import CHARGE_MODEL_CM5, _calc_cm5_charges, _charge_model_for_method
from ..errors import DescriptorUnavailableError, XTBEnvironmentError, XTBNativeNotFoundError, XTBRunError
from ..input import _effective_charge_uhf, _input_metadata, _method_flag, _xyz_from_mol
from ..options import _validate_calculation_settings, _validate_solvation_settings
from ..result import XTBResults, _frontier_energies
from ..runtime import (
    _engine_identity,
    _find_xtb_library,
    _native_properties_build,
    _native_realspace_build,
    _parameter_file,
    _require_extended_libxtb,
)

HARTREE_TO_EV = 27.211386245988
ANGSTROM_TO_BOHR = 1.8897261254578281


class _NativeXTB:
    """Small ctypes binding for the stable xTB C API."""

    _handle = ctypes.c_void_p

    def __init__(self, library: str):
        self.library = library
        self.lib = ctypes.CDLL(library)
        self.get_api_version = self._bind("xtb_getAPIVersion", ctypes.c_int, [])
        api_version = self.get_api_version()
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
        self.set_solvent = self._bind(
            'xtb_setSolvent',
            None,
            [
                h,
                h,
                ctypes.c_char_p,
                ctypes.POINTER(ctypes.c_int),
                ctypes.POINTER(ctypes.c_double),
                ctypes.POINTER(ctypes.c_int),
                ctypes.c_char_p,
            ],
        )
        self.release_solvent = self._bind('xtb_releaseSolvent', None, [h, h])
        self.new_results = self._bind("xtb_newResults", h, [])
        self.copy_results = self._bind("xtb_copyResults", h, [h])
        self.del_results = self._bind("xtb_delResults", None, [ctypes.POINTER(h)])
        self.singlepoint = self._bind("xtb_singlepoint", None, [h, h, h, h])
        self.get_energy = self._bind("xtb_getEnergy", None, [h, h, ctypes.POINTER(ctypes.c_double)])
        self.get_gradient = self._bind(
            "xtb_getGradient", None, [h, h, ctypes.POINTER(ctypes.c_double)]
        )
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

    def bind_esp_api(self):
        from ..runtime import _native_esp_build
        if not _native_esp_build(self.library):
            raise DescriptorUnavailableError('ESP requires the verified ESP extension; rebuild libxtb')
        h, i, d = self._handle, ctypes.c_int, ctypes.c_double
        if self._bind('xtb_rdkit_getESPAPIVersion', i, [])() != 1:
            raise DescriptorUnavailableError('unsupported native ESP API')
        self.set_keep_esp = self._bind('xtb_rdkit_setKeepESP', None, [h, h, ctypes.c_bool])
        self.get_esp = self._bind('xtb_rdkit_getESP', i,
                                  [h, h, i, ctypes.POINTER(d), i, ctypes.POINTER(d)])

    def esp_data(self, environment, results, points):
        count = len(points)
        xyz = (ctypes.c_double * (3 * count))(*(v * ANGSTROM_TO_BOHR for row in points for v in row))
        values = (ctypes.c_double * count)()
        ok = self.get_esp(environment, results, count, xyz, count, values)
        self.check(environment, 'ESP evaluation')
        if not ok:
            raise XTBEnvironmentError('native ESP evaluation failed without a diagnostic')
        return dict(espPointsAngstrom=points, espValuesAU=tuple(values))

    def bind_realspace_api(self):
        if not _native_realspace_build(self.library):
            raise DescriptorUnavailableError(
                'real-space evaluation requires the verified realspace extension; rebuild libxtb'
            )
        h, i, d = self._handle, ctypes.c_int, ctypes.c_double
        if self._bind('xtb_rdkit_getRealSpaceAPIVersion', i, [])() != 1:
            raise DescriptorUnavailableError('unsupported native real-space API')
        self.eval_realspace = self._bind(
            'xtb_rdkit_evalRealSpace',
            i,
            [
                h,
                h,
                i,
                ctypes.POINTER(d),
                d,
                d,
                ctypes.POINTER(d),
                ctypes.POINTER(d),
                ctypes.POINTER(d),
                ctypes.POINTER(d),
            ],
        )

    def realspace_data(self, environment, results, points, occ_threshold=1e-7, lea_window_ev=0.0):
        count = len(points)
        xyz = (ctypes.c_double * (3 * count))(*(v * ANGSTROM_TO_BOHR for row in points for v in row))
        density = (ctypes.c_double * count)()
        alie = (ctypes.c_double * count)()
        vir_density = (ctypes.c_double * count)()
        lea = (ctypes.c_double * count)()
        ok = self.eval_realspace(
            environment,
            results,
            count,
            xyz,
            float(occ_threshold),
            float(lea_window_ev),
            density,
            alie,
            vir_density,
            lea,
        )
        self.check(environment, 'real-space evaluation')
        if not ok:
            raise XTBEnvironmentError('native real-space evaluation failed without a diagnostic')
        return dict(
            realSpacePointsAngstrom=points,
            realSpaceDensityAU=tuple(float(v) for v in density),
            realSpaceALIEEV=tuple(float(v) for v in alie),
            realSpaceVirtualDensityAU=tuple(float(v) for v in vir_density),
            realSpaceLEAEV=tuple(float(v) for v in lea),
        )

    def bind_properties_api(self):
        if not _native_properties_build(self.library):
            raise DescriptorUnavailableError(
                'extended properties require the verified properties extension; rebuild libxtb'
            )
        h, i, d = self._handle, ctypes.c_int, ctypes.c_double
        if self._bind('xtb_rdkit_getPropertiesAPIVersion', i, [])() != 1:
            raise DescriptorUnavailableError('unsupported native properties API')
        self.get_energy_decomposition = self._bind(
            'xtb_rdkit_getEnergyDecomposition', i, [h, h, ctypes.POINTER(i), i, ctypes.POINTER(d)]
        )
        self.get_polarizability = self._bind(
            'xtb_rdkit_getPolarizability',
            i,
            [h, h, i, ctypes.POINTER(d), ctypes.POINTER(d), ctypes.POINTER(d)],
        )
        self.get_multipoles = self._bind(
            'xtb_rdkit_getMultipoles',
            i,
            [h, h, i, ctypes.POINTER(d), ctypes.POINTER(d), ctypes.POINTER(d)],
        )
        self.get_nuclear_esp = self._bind(
            'xtb_rdkit_getNuclearESP', i, [h, h, i, ctypes.POINTER(d)]
        )

    def properties_data(self, environment, results, natoms):
        out: dict[str, Any] = {}
        flags = ctypes.c_int(0)
        evals = (ctypes.c_double * 15)()
        if self.get_energy_decomposition(environment, results, ctypes.byref(flags), 15, evals):
            f = int(flags.value)
            out['energyDecompositionFlags'] = f
            out['electronicEnergyHartree'] = float(evals[1]) if (f & 2) else None
            out['repulsionEnergyHartree'] = float(evals[2]) if (f & 4) else None
            out['electrostaticEnergyHartree'] = float(evals[3]) if (f & 8) else None
            out['anisotropicElectrostaticEnergyHartree'] = float(evals[4]) if (f & 16) else None
            out['anisotropicXCEnergyHartree'] = float(evals[5]) if (f & 32) else None
            out['dispersionEnergyHartree'] = float(evals[6]) if (f & 64) else None
            out['halogenBondEnergyHartree'] = float(evals[7]) if (f & 128) else None
            out['solvationEnergyHartree'] = float(evals[8]) if (f & 256) else None
            if f & 512:
                out['fermiLevelAlphaHartree'] = float(evals[9])
                out['fermiLevelBetaHartree'] = float(evals[10])
                out['fermiLevelAlphaEV'] = float(evals[9] * HARTREE_TO_EV)
                out['fermiLevelBetaEV'] = float(evals[10] * HARTREE_TO_EV)
            if f & 1024:
                out['solvationBornEnergyHartree'] = float(evals[11])
                out['solvationSASAEnergyHartree'] = float(evals[12])
                out['solvationHBEnergyHartree'] = float(evals[13])
                out['solvationShiftEnergyHartree'] = float(evals[14])
        mol_pol = (ctypes.c_double * 3)()
        alpha_atom = (ctypes.c_double * natoms)()
        c6_atom = (ctypes.c_double * natoms)()
        if self.get_polarizability(environment, results, natoms, mol_pol, alpha_atom, c6_atom):
            out['molecularPolarizabilityAU'] = float(mol_pol[0])
            out['molecularC6AU'] = float(mol_pol[1])
            out['molecularC8AU'] = float(mol_pol[2])
            out['atomicPolarizabilitiesAU'] = tuple(float(v) for v in alpha_atom)
            out['atomicC6AU'] = tuple(float(v) for v in c6_atom)
        dip_atom = (ctypes.c_double * (3 * natoms))()
        qp_atom = (ctypes.c_double * (6 * natoms))()
        mol_qp = (ctypes.c_double * 6)()
        if self.get_multipoles(environment, results, natoms, dip_atom, qp_atom, mol_qp):
            out['atomicDipolesAU'] = tuple(
                tuple(float(dip_atom[3 * i + k]) for k in range(3)) for i in range(natoms)
            )
            out['atomicQuadrupolesAU'] = tuple(
                tuple(float(qp_atom[6 * i + k]) for k in range(6)) for i in range(natoms)
            )
            out['molecularQuadrupoleAU'] = tuple(float(v) for v in mol_qp)
        nuclear_esp = (ctypes.c_double * natoms)()
        epn_ok = self.get_nuclear_esp(environment, results, natoms, nuclear_esp)
        self.check(environment, 'nuclear ESP evaluation')
        if not epn_ok:
            raise XTBEnvironmentError('native nuclear ESP evaluation failed without a diagnostic')
        out['nuclearESPAU'] = tuple(float(v) for v in nuclear_esp)
        return out

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
        from ..basis import _native_p_directions
        return dict(
            aoPDirections=_native_p_directions(tuple(ash), tuple(sl)),
            overlapMatrix=tuple(tuple(overlap[j * nao + i] for j in range(nao)) for i in range(nao)),
            aoAtomMap=tuple(aa),
            aoShellMap=tuple(ash),
            shellAtomMap=tuple(sa),
            shellAngularMomenta=tuple(sl),
        )

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
    confId: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    numUnpairedElectrons: int | None = None,
    xtbLibrary: str | os.PathLike[str] | None = None,
    keepOverlap: bool = False,
    espPointsAngstrom=(),
    realSpacePointsAngstrom=(),
    realSpaceOccupationThreshold: float = 1e-7,
    realSpaceLEAWindowEV: float = 0.0,
    maxIterations: int = 250,
    accuracy: float = 1.0,
    electronicTemperature: float = 300.0,
    solvent: str | None = None,
    solvationModel: str = "alpb",
    solvationGrid: int | None = None,
    _restart=None,
    _save_restart=None,
    _return_gradient=False,
    _keep_esp_for_hook: bool = False,
    _post_singlepoint_hook=None,
) -> XTBResults:
    """Run a single point through a shared ``libxtb`` C API library."""
    import math
    from numbers import Real
    from rdkit import Chem
    from ..esp_input import _normalize_esp_points
    espPointsAngstrom = _normalize_esp_points(espPointsAngstrom)
    realSpacePointsAngstrom = _normalize_esp_points(realSpacePointsAngstrom)
    for label, val in (
        ('realSpaceOccupationThreshold', realSpaceOccupationThreshold),
        ('realSpaceLEAWindowEV', realSpaceLEAWindowEV),
    ):
        if isinstance(val, bool) or not isinstance(val, Real) or not math.isfinite(val) or val < 0:
            raise ValueError(f'{label} must be a finite nonnegative real number')
    realSpaceOccupationThreshold = float(realSpaceOccupationThreshold)
    realSpaceLEAWindowEV = float(realSpaceLEAWindowEV)
    settings = _validate_calculation_settings(maxIterations, accuracy, electronicTemperature)
    solvent_val, solv_model, solv_grid = _validate_solvation_settings(
        solvent, solvationModel, solvationGrid
    )
    eff_solv_grid = (230 if solv_grid is None else solv_grid) if solvent_val is not None else None
    mol = Chem.Mol(mol)  # Own the exact geometry/topology while ctypes releases the GIL.
    if not isinstance(keepOverlap, bool):
        raise ValueError('keepOverlap must be a bool')
    xyz = _xyz_from_mol(mol, confId)
    method_name, _ = _method_flag(method)
    if method_name == 'GFN-FF':
        raise ValueError('GFN-FF does not provide molecular orbital descriptors')
    charge, unpaired = _effective_charge_uhf(mol, charge, numUnpairedElectrons)
    library = _NativeXTB(_find_xtb_library(xtbLibrary))
    _require_extended_libxtb(library.library)
    engine = _engine_identity('native', library.library)
    need_esp = bool(espPointsAngstrom) or bool(_keep_esp_for_hook)
    if need_esp:
        library.bind_esp_api()
    if realSpacePointsAngstrom:
        library.bind_realspace_api()
    if keepOverlap:
        if not engine['orbital_basis']:
            raise DescriptorUnavailableError(
                'keepOverlap requires a verified native basis extension; rebuild with tools/build_patched_xtb.py'
            )
        library.bind_basis_api()
    if method_name == 'GFN0-xTB' and not engine['native_gfn0']:
        raise ValueError('native GFN0 requires a verified patched xTB build; run tools/build_patched_xtb.sh')
    parameter = _parameter_file(method_name, library.library)
    started = time.perf_counter()
    natoms = mol.GetNumAtoms()
    conf = mol.GetConformer(confId)
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
    uhf_value = ctypes.c_int(unpaired)
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
        if solvent_val is not None:
            grid_arg = ctypes.c_int(eff_solv_grid)
            library.set_solvent(
                environment,
                calculator,
                solvent_val.encode('utf-8'),
                None,
                None,
                ctypes.byref(grid_arg),
                solv_model.encode('utf-8'),
            )
            library.check(environment, "solvation setup")
        library.set_max_iterations(environment, calculator, settings['max_iterations'])
        library.set_accuracy(environment, calculator, settings['accuracy'])
        library.set_temperature(environment, calculator, settings['electronic_temperature'])
        library.check(environment, 'calculator settings')
        results = library.copy_results(_restart) if _restart else library.new_results()
        if not results:
            if _restart:
                raise XTBEnvironmentError('native xTB could not copy the reference wavefunction')
            raise XTBEnvironmentError("native xTB could not create a results object")
        if keepOverlap:
            library.set_keep_overlap(environment, results, True)
            library.check(environment, "basis retention setup")
        if need_esp:
            library.set_keep_esp(environment, results, True)
            library.check(environment, "ESP retention setup")
        library.singlepoint(environment, molecule, calculator, results)
        library.check(environment, "single point")

        energy = ctypes.c_double()
        gradient = (ctypes.c_double * (3 * natoms))() if _return_gradient else None
        dipole = (ctypes.c_double * 3)()
        charges = (ctypes.c_double * natoms)()
        bond_orders = (ctypes.c_double * (natoms * natoms))()
        nao = ctypes.c_int()
        library.get_energy(environment, results, ctypes.byref(energy))
        if gradient is not None:
            library.get_gradient(environment, results, gradient)
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
        esp = library.esp_data(environment, results, espPointsAngstrom) if espPointsAngstrom else {}
        realspace = (
            library.realspace_data(
                environment,
                results,
                realSpacePointsAngstrom,
                occ_threshold=realSpaceOccupationThreshold,
                lea_window_ev=realSpaceLEAWindowEV,
            )
            if realSpacePointsAngstrom
            else {}
        )
        if _post_singlepoint_hook is not None:
            hook_out = _post_singlepoint_hook(library, environment, results, natoms, nao.value)
            if isinstance(hook_out, dict):
                if 'esp' in hook_out:
                    esp = hook_out['esp']
                if 'realspace' in hook_out:
                    realspace = hook_out['realspace']
        basis = library.basis_data(environment, results, natoms, nao.value) if keepOverlap else {}
        props = {}
        if _native_properties_build(library.library):
            library.bind_properties_api()
            props = library.properties_data(environment, results, natoms)
        from ..runtime import _native_diagnostics_build
        from ..scf import scf_record
        scf = scf_record(method_name)
        if _native_diagnostics_build(library.library):
            version = library._bind('xtb_rdkit_getSCFAPIVersion', ctypes.c_int, [])
            getter = library._bind(
                'xtb_rdkit_getSCFDiagnostics',
                ctypes.c_int,
                [library._handle, library._handle, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)],
            )
            if version() != 1:
                raise XTBEnvironmentError('unsupported SCC diagnostics API')
            count = ctypes.c_int(-1)
            data = (ctypes.c_double * 4)()
            if getter(environment, results, ctypes.byref(count), data):
                scf = scf_record(method_name, count.value, list(data))
        if _save_restart is not None:
            saved = library.copy_results(results)
            if not saved:
                raise XTBEnvironmentError('native xTB could not copy the reference wavefunction')
            _save_restart.append((library, saved))
        grad_rows = (
            tuple(
                tuple(float(gradient[3 * row + axis]) for axis in range(3))
                for row in range(natoms)
            )
            if gradient is not None
            else ()
        )
        return XTBResults(
            **esp,
            **realspace,
            **basis,
            **props,
            solvent=solvent_val,
            solvationModel=(solv_model if solvent_val is not None else None),
            solvationGrid=eff_solv_grid,
            gradientHartreePerBohr=grad_rows,
            method=method_name,
            version=engine["version"],
            totalEnergyHartree=float(energy.value),
            homoLumoGapEV=(lumo - homo if homo is not None and lumo is not None else None),
            numOrbitals=int(nao.value),
            numElectrons=int(round(sum(occupation_values))),
            numUnpairedElectrons=int(unpaired),
            orbitalEnergiesEV=energies_ev,
            orbitalOccupancies=occupation_values,
            atomicCharges=tuple(float(value) for value in charges),
            chargeModel=_charge_model_for_method(method_name),
            cm5Charges=(
                _calc_cm5_charges(
                    [atom.GetAtomicNum() for atom in mol.GetAtoms()],
                    [
                        [conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y, conf.GetAtomPosition(i).z]
                        for i in range(natoms)
                    ],
                    tuple(float(value) for value in charges),
                )
                if method_name in ('GFN1-xTB', 'GFN2-xTB')
                else ()
            ),
            dipoleAtomicUnits=tuple(float(value) for value in dipole),
            wibergBondOrders=tuple(
                tuple(float(bond_orders[row * natoms + column]) for column in range(natoms))
                for row in range(natoms)
            ),
            orbitalCoefficients=tuple(
                tuple(float(coefficients[column * nao.value + row]) for column in range(nao.value))
                for row in range(nao.value)
            ),
            backend="native",
            rawJSON={"native_api_version": int(library.get_api_version())},
            metadata={
                **_input_metadata(mol, confId, charge, unpaired, xyz),
                'calculation_settings': settings,
                'calculation_status': 'success',
                'scf': scf,
                'initialization': 'reference' if _restart else 'independent',
                'optimize': False,
                'solvent': (
                    None
                    if solvent_val is None
                    else {'solvent': solvent_val, 'model': solv_model, 'grid': eff_solv_grid}
                ),
                **({'properties_api': 1} if props else {}),
                'engine': engine,
                'parameter_file': parameter,
                'native_api_version': int(library.get_api_version()),
                'elapsed_seconds': time.perf_counter() - started,
                'coefficient_axes': ['AO', 'MO'],
                'charge_model': _charge_model_for_method(method_name),
                'cm5_model': CHARGE_MODEL_CM5 if method_name in ('GFN1-xTB', 'GFN2-xTB') else None,
                **(
                    {
                        'esp': {
                            'api_version': 1,
                            'potential_unit': 'hartree/e',
                            'coordinate_unit': 'angstrom',
                            'source': 'native_density',
                        }
                    }
                    if esp
                    else {}
                ),
                **(
                    {
                        'realspace': {
                            'api_version': 1,
                            'density_unit': 'e/bohr^3',
                            'alie_unit': 'eV',
                            'virtual_density_unit': '1/bohr^3',
                            'lea_unit': 'eV',
                            'coordinate_unit': 'angstrom',
                            'occupation_threshold': realSpaceOccupationThreshold,
                            'lea_window_ev': realSpaceLEAWindowEV,
                            'source': 'native_cartesian_gto',
                        }
                    }
                    if realspace
                    else {}
                ),
                'overlap_source': 'same libxtb calculation' if keepOverlap else None,
                **(
                    {'gradient_hartree_per_bohr': tuple(float(value) for value in gradient)}
                    if gradient is not None
                    else {}
                ),
            },
        )
    finally:
        library._delete(library.del_results, results)
        library._delete(library.del_calculator, calculator)
        library._delete(library.del_molecule, molecule)
        library._delete(library.del_environment, environment)


def run_xtb_native(
    mol,
    *,
    confId=-1,
    method='GFN2-xTB',
    charge=None,
    numUnpairedElectrons=None,
    xtbLibrary=None,
    keepOverlap=False,
    espPointsAngstrom=(),
    realSpacePointsAngstrom=(),
    realSpaceOccupationThreshold=1e-7,
    realSpaceLEAWindowEV=0.0,
    maxIterations=250,
    accuracy=1.0,
    electronicTemperature=300.0,
    solvent=None,
    solvationModel='alpb',
    solvationGrid=None,
):
    """Run an independently initialized libxtb single point."""
    return _run_xtb_native(
        mol,
        confId=confId,
        method=method,
        charge=charge,
        numUnpairedElectrons=numUnpairedElectrons,
        xtbLibrary=xtbLibrary,
        keepOverlap=keepOverlap,
        espPointsAngstrom=espPointsAngstrom,
        realSpacePointsAngstrom=realSpacePointsAngstrom,
        realSpaceOccupationThreshold=realSpaceOccupationThreshold,
        realSpaceLEAWindowEV=realSpaceLEAWindowEV,
        maxIterations=maxIterations,
        accuracy=accuracy,
        electronicTemperature=electronicTemperature,
        solvent=solvent,
        solvationModel=solvationModel,
        solvationGrid=solvationGrid,
    )


def _native_charge_states(mol, configurations):
    """Scoped reference restart; callers must close this iterator on early exit."""
    saved = []
    with _native_lock:
        try:
            for index, config in enumerate(configurations):
                options = config.ToDict()
                options.pop('backend')
                yield _run_xtb_native(
                    mol,
                    **options,
                    _restart=saved[0][1] if index else None,
                    _save_restart=saved if not index else None,
                )
        finally:
            for library, handle in saved:
                library._delete(library.del_results, handle)
