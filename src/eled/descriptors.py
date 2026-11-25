"""Optional adapters for RDKit EHT and xTB scalar orbital features.

The main public descriptor protocol is xTB-based and lives in the dedicated
modules exported by :mod:`eled`. ``CalcEHTDescriptors`` and the EHT branches of
``CalcOrbitalDescriptors`` only normalize RDKit's existing YAeHMOP interface
for comparison; they do not modify or replace YAeHMOP.
"""

from __future__ import annotations

from typing import Any, Mapping

from .calculator import XTBCalculator
from .xtb import CalcXTBDescriptors, XTBResults


def _frontier_from_fermi(energies: tuple[float, ...], fermi: float) -> tuple[float | None, float | None]:
    occupied = [energy for energy in energies if energy <= fermi + 1.0e-7]
    virtual = [energy for energy in energies if energy > fermi + 1.0e-7]
    return (occupied[-1] if occupied else None, virtual[0] if virtual else None)


def CalcEHTDescriptors(result: Any) -> dict[str, float | int | str | None]:
    """Extract stable scalar descriptors from an RDKit ``EHTResults`` object."""
    energies = tuple(float(value) for value in result.GetOrbitalEnergies())
    occupied_count = (int(result.numElectrons) + 1) // 2
    if 0 < occupied_count < len(energies):
        homo, lumo = energies[occupied_count - 1], energies[occupied_count]
    else:
        homo, lumo = _frontier_from_fermi(energies, float(result.fermiEnergy))
    return {
        "eht_backend": "yaehmop",
        "eht_homo_ev": homo,
        "eht_lumo_ev": lumo,
        "eht_gap_ev": lumo - homo if homo is not None and lumo is not None else None,
        "eht_fermi_ev": float(result.fermiEnergy),
        "eht_energy_ev": float(result.totalEnergy),
        "eht_n_orb": int(result.numOrbitals),
        "eht_n_e": int(result.numElectrons),
    }


def CalcOrbitalDescriptors(
    mol: Any,
    *,
    backend: str = "xtb",
    confId: int = -1,
    keepMatrices: bool = False,
    xtbBackend: str = "native",
    frontierDepth: int = 2,
    occupationThreshold: float = 1e-7,
    **xtbOptions: Any,
) -> Mapping[str, float | int | str | None]:
    """Calculate method-prefixed scalar orbital descriptors.

    ``backend`` can be ``"eht"``, ``"xtb"``, or ``"both"``. ``xtbOptions`` are
    forwarded to :class:`XTBCalculator`; EHT always uses RDKit's YAeHMOP interface.
    """
    normalized = backend.strip().lower()
    if normalized not in {"eht", "xtb", "both"}:
        raise ValueError("backend must be 'eht', 'xtb', or 'both'")
    descriptors: dict[str, float | int | str | None] = {}
    if normalized in {"eht", "both"}:
        from rdkit.Chem import rdEHTTools

        ok, result = rdEHTTools.RunMol(
            mol,
            confId=confId,
            keepOverlapAndHamiltonianMatrices=keepMatrices,
        )
        if not ok:
            raise RuntimeError("RDKit YAeHMOP EHT calculation failed")
        descriptors.update(CalcEHTDescriptors(result))
    if normalized in {"xtb", "both"}:
        result: XTBResults = XTBCalculator(
            confId=confId,
            backend=xtbBackend,
            **xtbOptions,
        ).Calculate(mol)
        descriptors.update(
            CalcXTBDescriptors(
                result,
                frontierDepth=frontierDepth,
                occupationThreshold=occupationThreshold,
            )
        )
    return descriptors
