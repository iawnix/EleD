"""Result model shared by xTB backends.

Property and getter conventions match ``rdXTBTools.XTBResults`` and RDKit's
``EHTResults``: scalar attributes use ``camelCase`` and array getters use
``PascalCase`` ``Get*()`` methods returning independent NumPy copies.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field, fields
from typing import Any, Mapping, Sequence

import numpy as np


def _frontier_energies(
    energies: Sequence[float], occupations: Sequence[float], threshold: float = 1.0e-7
) -> tuple[float | None, float | None]:
    """Return HOMO and LUMO from xTB occupations."""
    if not energies or len(energies) != len(occupations):
        return None, None
    occupied = [energy for energy, occupation in zip(energies, occupations) if occupation > threshold]
    virtual = [energy for energy, occupation in zip(energies, occupations) if occupation <= threshold]
    return (occupied[-1] if occupied else None, virtual[0] if virtual else None)


@dataclass(frozen=True)
class XTBResults:
    """Normalized output from one gas-phase xTB single point.

    ``totalEnergyHartree`` uses Hartree; ``orbitalEnergiesEV`` uses eV;
    coordinates use Angstrom; dipoles use ``e*Bohr``. ``atomicCharges`` hold
    engine charges whose exact model is recorded in ``chargeModel``
    (``xtb_scc_mulliken`` for GFN1/2, ``gfn0_eeq`` for GFN0). ``cm5Charges``
    holds derived CM5 charges when available.
    """

    method: str
    version: str | None = None
    totalEnergyHartree: float | None = None
    homoLumoGapEV: float | None = None
    numOrbitals: int | None = None
    numElectrons: int | None = None
    numUnpairedElectrons: int | None = None
    orbitalEnergiesEV: tuple[float, ...] = ()
    orbitalOccupancies: tuple[float, ...] = ()
    atomicCharges: tuple[float, ...] = ()
    chargeModel: str = "unknown"
    cm5Charges: tuple[float, ...] = ()
    dipoleAtomicUnits: tuple[float, ...] = ()
    wibergBondOrders: tuple[tuple[float, ...], ...] = ()
    rawJSON: Mapping[str, Any] = field(default_factory=dict, repr=False)
    orbitalCoefficients: tuple[tuple[float, ...], ...] = ()
    backend: str = "native"
    metadata: Mapping[str, Any] = field(default_factory=dict, repr=False)

    # Same-calculation native basis snapshot (opt-in via keepOverlap=True).
    overlapMatrix: tuple[tuple[float, ...], ...] = ()
    aoAtomMap: tuple[int, ...] = ()
    aoShellMap: tuple[int, ...] = ()
    shellAtomMap: tuple[int, ...] = ()
    shellAngularMomenta: tuple[int, ...] = ()
    aoPDirections: tuple[tuple[float, ...], ...] = ()

    # Optional nuclear gradient [atom][xyz] in Hartree/Bohr.
    gradientHartreePerBohr: tuple[tuple[float, ...], ...] = ()

    # Sampled ESP only; no native pointers. Empty tuples mean not requested.
    espPointsAngstrom: tuple[tuple[float, ...], ...] = ()
    espValuesAU: tuple[float, ...] = ()

    def __post_init__(self) -> None:
        if self.chargeModel == "unknown":
            model = {
                "GFN0-xTB": "gfn0_eeq",
                "GFN1-xTB": "xtb_scc_mulliken",
                "GFN2-xTB": "xtb_scc_mulliken",
            }.get(self.method)
            if model is not None:
                object.__setattr__(self, "chargeModel", model)

    def _require_success(self) -> None:
        if not self.success:
            detail = self.metadata.get("error_message", "calculation did not succeed")
            raise RuntimeError(f"No numerical results: {detail}")

    @property
    def success(self) -> bool:
        return self.metadata.get("calculation_status", "success") == "success"

    @property
    def conformerId(self) -> int:
        return int(self.metadata.get("conformer_id", -1))

    @property
    def numAtoms(self) -> int:
        if "atom_count" in self.metadata and self.metadata["atom_count"] is not None:
            return int(self.metadata["atom_count"])
        return len(self.atomicCharges)

    @property
    def charge(self) -> int | None:
        value = self.metadata.get("charge")
        return None if value is None else int(value)

    @property
    def overlapRetained(self) -> bool:
        return self.success and bool(self.overlapMatrix)

    @property
    def espSampled(self) -> bool:
        return self.success and bool(self.espValuesAU)

    @property
    def homoEV(self) -> float | None:
        homo, _ = _frontier_energies(self.orbitalEnergiesEV, self.orbitalOccupancies)
        return homo

    @property
    def lumoEV(self) -> float | None:
        _, lumo = _frontier_energies(self.orbitalEnergiesEV, self.orbitalOccupancies)
        return lumo

    @property
    def gapEV(self) -> float | None:
        homo, lumo = self.homoEV, self.lumoEV
        if homo is not None and lumo is not None:
            return lumo - homo
        return self.homoLumoGapEV

    @property
    def fermiEnergyEV(self) -> float | None:
        """Use the occupied frontier energy as the discrete Fermi estimate."""
        return self.homoEV

    def GetAtomicCharges(self) -> np.ndarray:
        """Copy of xTB atomic charges in elementary-charge units, [atom]."""
        self._require_success()
        return np.array(self.atomicCharges, dtype=float, copy=True)

    def GetCM5Charges(self) -> np.ndarray:
        """Copy of derived CM5 charges in elementary-charge units, [atom]."""
        self._require_success()
        if self.chargeModel != "xtb_scc_mulliken" or not len(self.cm5Charges):
            raise ValueError("CM5 charges require GFN1/2 SCC Mulliken charges; GFN0 exposes PEEQ/EEQ")
        return np.array(self.cm5Charges, dtype=float, copy=True)

    def GetOrbitalEnergies(self) -> np.ndarray:
        """Copy of molecular orbital energies in eV, [MO]."""
        self._require_success()
        return np.array(self.orbitalEnergiesEV, dtype=float, copy=True)

    def GetOrbitalOccupancies(self) -> np.ndarray:
        """Copy of orbital occupancies in electrons, [MO]."""
        self._require_success()
        return np.array(self.orbitalOccupancies, dtype=float, copy=True)

    def GetDipole(self) -> np.ndarray:
        """Copy of dipole in e*Bohr, [x, y, z]."""
        self._require_success()
        return np.array(self.dipoleAtomicUnits, dtype=float, copy=True)

    def GetGradient(self) -> np.ndarray:
        """Copy of nuclear gradient in Hartree/Bohr, [atom, xyz]."""
        self._require_success()
        if len(self.gradientHartreePerBohr):
            return np.array(self.gradientHartreePerBohr, dtype=float, copy=True)
        stored = self.metadata.get("gradient_hartree_per_bohr")
        if stored is not None and len(stored):
            return np.array(stored, dtype=float, copy=True).reshape(-1, 3)
        raise ValueError("Gradient was not retained")

    def GetWibergBondOrders(self) -> np.ndarray:
        """Copy of dimensionless Wiberg bond orders, [atom, atom]."""
        self._require_success()
        return np.array(self.wibergBondOrders, dtype=float, copy=True)

    def GetOrbitalCoefficients(self) -> np.ndarray:
        """Copy of optional MO coefficients, [AO, MO]."""
        self._require_success()
        if not len(self.orbitalCoefficients):
            raise ValueError("Orbital coefficients were not retained; set keepOverlap=True")
        return np.array(self.orbitalCoefficients, dtype=float, copy=True)

    def GetOverlapMatrix(self) -> np.ndarray:
        """Copy of same-calculation AO overlap matrix, [AO, AO]."""
        self._require_success()
        if not len(self.overlapMatrix):
            raise ValueError("Overlap was not retained; set keepOverlap=True")
        return np.array(self.overlapMatrix, dtype=float, copy=True)

    def GetAOAtomMap(self) -> np.ndarray:
        self._require_success()
        if not len(self.aoAtomMap):
            raise ValueError("Basis maps were not retained; set keepOverlap=True")
        return np.array(self.aoAtomMap, dtype=np.int32, copy=True)

    def GetAOShellMap(self) -> np.ndarray:
        self._require_success()
        if not len(self.aoShellMap):
            raise ValueError("Basis maps were not retained; set keepOverlap=True")
        return np.array(self.aoShellMap, dtype=np.int32, copy=True)

    def GetShellAtomMap(self) -> np.ndarray:
        self._require_success()
        if not len(self.shellAtomMap):
            raise ValueError("Basis maps were not retained; set keepOverlap=True")
        return np.array(self.shellAtomMap, dtype=np.int32, copy=True)

    def GetShellAngularMomenta(self) -> np.ndarray:
        self._require_success()
        if not len(self.shellAngularMomenta):
            raise ValueError("Basis maps were not retained; set keepOverlap=True")
        return np.array(self.shellAngularMomenta, dtype=np.int32, copy=True)

    def GetAOPDirections(self) -> np.ndarray:
        self._require_success()
        if not len(self.aoPDirections):
            raise ValueError("AO p directions were not retained; set keepOverlap=True")
        return np.array(self.aoPDirections, dtype=float, copy=True)

    def GetESPPoints(self) -> np.ndarray:
        """Copy of ESP sampling points in Angstrom, [point, xyz]."""
        self._require_success()
        if not len(self.espValuesAU) or not len(self.espPointsAngstrom):
            raise ValueError("ESP was not sampled; pass espPointsAngstrom")
        return np.array(self.espPointsAngstrom, dtype=float, copy=True)

    def GetESPValues(self) -> np.ndarray:
        """Copy of sampled ESP in Hartree/e, [point]."""
        self._require_success()
        if not len(self.espValuesAU):
            raise ValueError("ESP was not sampled; pass espPointsAngstrom")
        return np.array(self.espValuesAU, dtype=float, copy=True)

    def ToDict(self) -> dict[str, Any]:
        """Return a JSON-serializable dictionary of owned result fields."""
        return asdict(self)

    @classmethod
    def FromDict(cls, data: Mapping[str, Any]) -> XTBResults:
        """Reconstruct ``XTBResults`` from :meth:`ToDict` output."""
        payload = dict(data)

        def to_tuples(value: Any) -> Any:
            return tuple(to_tuples(item) for item in value) if isinstance(value, list) else value

        for item in fields(cls):
            if item.name in payload and isinstance(item.default, tuple):
                payload[item.name] = to_tuples(payload[item.name])
        return cls(**payload)


XTBResult = XTBResults
