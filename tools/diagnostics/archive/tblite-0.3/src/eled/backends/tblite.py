"""Optional in-process tblite backend for RDKit-to-xTB comparisons.

tblite is deliberately imported lazily.  Installing the base package does not
install tblite; this backend is selected only with ``backend='tblite'``.
"""

from __future__ import annotations

import hashlib
import importlib.metadata
import importlib.util
import math
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

from ..errors import TBLiteNotFoundError, XTBRunError
from ..input import _effective_charge_uhf, _input_metadata, _method_flag, _xyz_from_mol
from ..orbital_alignment import validate_tblite_basis
from ..result import XTBResult, _frontier_energies

ANGSTROM_TO_BOHR = 1.8897261254578281
HARTREE_TO_EV = 27.211386245988
KELVIN_TO_HARTREE = 3.166811563e-6


def _load_tblite():
    try:
        from tblite.interface import Calculator
    except (ImportError, OSError) as exc:
        raise TBLiteNotFoundError(
            "tblite backend requires the optional dependency; install "
            "from the source checkout with pip install -e '.[tblite]' "
            "or conda install tblite-python"
        ) from exc
    return Calculator


@lru_cache(maxsize=1)
def _tblite_identity() -> dict[str, Any]:
    """Return a cacheable identity for the optional tblite wheel/library."""
    try:
        version = importlib.metadata.version("tblite")
    except importlib.metadata.PackageNotFoundError as exc:
        raise TBLiteNotFoundError("tblite is not installed") from exc
    spec = importlib.util.find_spec("tblite._libtblite")
    path = Path(spec.origin).resolve() if spec and spec.origin else None
    digest = None
    if path and path.is_file():
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "backend": "tblite",
        "package": "tblite",
        "version": version,
        "path": str(path) if path else None,
        "sha256": digest,
    }


def _one_spin(values: Any, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 1:
        return array
    if array.ndim == 2 and array.shape[0] == 1:
        return array[0]
    raise XTBRunError(
        f"tblite returned spin-polarized {name}; the current XTBResult schema "
        "does not represent separate alpha and beta channels"
    )


def _matrix(values: Any, name: str, natoms: int) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.ndim == 3 and array.shape[-1] == 1:
        array = array[..., 0]
    if array.shape != (natoms, natoms):
        raise XTBRunError(f"tblite returned an invalid {name} shape {array.shape}")
    return array


def run_xtb_tblite(
    mol: Any,
    *,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    uhf: int | None = None,
    save_integrals: bool = False,
    accuracy: float = 1.0,
) -> XTBResult:
    """Run GFN1/2 through tblite's in-process library interface.

    The result uses the same normalized :class:`XTBResult` units as the
    ctypes backend.  tblite's spin-polarized result arrays are rejected
    until the normalized result model grows explicit spin-channel fields.
    Diagnostic callers can retain the overlap in raw_json with save_integrals.
    """
    xyz = _xyz_from_mol(mol, conf_id)
    if not math.isfinite(accuracy) or accuracy <= 0:
        raise ValueError("accuracy must be finite and positive")
    method_name, _ = _method_flag(method)
    if method_name == "GFN0-xTB":
        raise ValueError("tblite does not provide GFN0-xTB; use backend='native'")
    charge, uhf = _effective_charge_uhf(mol, charge, uhf)
    if uhf:
        raise XTBRunError(
            "tblite spin-polarized results are not yet representable by XTBResult; "
            "use backend='native' for radicals"
        )
    try:
        Calculator = _load_tblite()
        coordinates = np.asarray(mol.GetConformer(conf_id).GetPositions(), dtype=float)
        numbers = np.asarray([atom.GetAtomicNum() for atom in mol.GetAtoms()], dtype=np.int32)
        calculator = Calculator(
            method_name,
            numbers,
            coordinates * ANGSTROM_TO_BOHR,
            charge=charge,
            uhf=uhf,
        )
        # Match the native adapter's defaults: xTB accuracy 1 and 300 K filling.
        for attribute, value in (
            ("accuracy", accuracy),
            ("max-iter", 250),
            ("temperature", 300.0 * KELVIN_TO_HARTREE),
            ("verbosity", 0),
        ):
            calculator.set(attribute, value)
        if save_integrals:
            calculator.set("save-integrals", True)
        result = calculator.singlepoint()
        energy = float(np.asarray(result.get("energy")))
        orbital_energies = _one_spin(result.get("orbital-energies"), "orbital energies")
        occupations = _one_spin(result.get("orbital-occupations"), "orbital occupations")
        coefficients = np.asarray(result.get("orbital-coefficients"), dtype=float)
        if coefficients.ndim == 3 and coefficients.shape[0] == 1:
            coefficients = coefficients[0]
        if coefficients.ndim != 2:
            raise XTBRunError(f"tblite returned an invalid orbital coefficient shape {coefficients.shape}")
        charges = np.asarray(result.get("charges"), dtype=float)
        dipole = np.asarray(result.get("dipole"), dtype=float)
        natoms = mol.GetNumAtoms()
        bond_orders = _matrix(result.get("bond-orders"), "bond-order", natoms)
        if charges.shape != (natoms,) or dipole.shape != (3,):
            raise XTBRunError("tblite returned an invalid charge or dipole shape")
        if coefficients.shape != (len(orbital_energies), len(orbital_energies)):
            raise XTBRunError("tblite returned inconsistent orbital dimensions")
        values = [energy, *orbital_energies, *occupations, *charges, *dipole, *bond_orders.ravel(), *coefficients.ravel()]
        if not all(math.isfinite(float(value)) for value in values):
            raise XTBRunError("tblite returned non-finite numerical results")
        if len(orbital_energies) != len(occupations):
            raise XTBRunError("tblite returned inconsistent orbital energies and occupations")
        homo, lumo = _frontier_energies(
            tuple(float(value * HARTREE_TO_EV) for value in orbital_energies),
            tuple(float(value) for value in occupations),
        )
        orbital_map, shell_map, angular_momenta = validate_tblite_basis(
            calculator.get("orbital-map"), calculator.get("shell-map"),
            calculator.get("angular-momenta"), natoms, len(orbital_energies))
        integral_data = {}
        if save_integrals:
            overlap = np.asarray(result.get("overlap-matrix"), dtype=float)
            if (overlap.shape != coefficients.shape or not np.isfinite(overlap).all()
                    or not np.allclose(overlap, overlap.T, atol=1e-12, rtol=0)):
                raise XTBRunError("tblite returned an invalid overlap matrix")
            integral_data["overlap_matrix"] = overlap.tolist()
        identity = dict(_tblite_identity())
        return XTBResult(
            method=method_name,
            version=identity["version"],
            total_energy_eh=energy,
            homo_lumo_gap_ev=(lumo - homo if homo is not None and lumo is not None else None),
            number_of_molecular_orbitals=len(orbital_energies),
            number_of_electrons=int(round(float(occupations.sum()))),
            number_of_unpaired_electrons=uhf,
            orbital_energies_ev=tuple(float(value * HARTREE_TO_EV) for value in orbital_energies),
            orbital_occupations=tuple(float(value) for value in occupations),
            partial_charges=tuple(float(value) for value in charges),
            dipole_au=tuple(float(value) for value in dipole),
            wiberg_bond_orders=tuple(tuple(float(value) for value in row) for row in bond_orders),
            orbital_coefficients=tuple(
                tuple(float(value) for value in row) for row in coefficients
            ),
            overlap_matrix=tuple(map(tuple, integral_data["overlap_matrix"])) if save_integrals else (),
            ao_atom_indices=tuple(int(v) for v in shell_map[orbital_map]),
            backend="tblite",
            raw_json={
                "energy_unit": "Hartree",
                "orbital_energy_unit": "Hartree",
                "orbital_map": orbital_map.tolist(),
                "shell_map": shell_map.tolist(),
                "angular_momenta": angular_momenta.tolist(),
                **integral_data,
            },
            metadata={
                **_input_metadata(mol, conf_id, charge, uhf, xyz),
                "optimize": False,
                "solvent": None,
                "engine": identity,
                "coefficient_axes": ["AO", "MO"],
                "overlap_source": "same tblite calculation" if save_integrals else None,
                "tblite_temperature_hartree": 300.0 * KELVIN_TO_HARTREE,
                "accuracy": accuracy,
            },
        )
    except (TBLiteNotFoundError, XTBRunError):
        raise
    except Exception as exc:
        raise XTBRunError(f"tblite single point failed: {exc}") from exc
