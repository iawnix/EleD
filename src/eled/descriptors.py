"""Optional adapters for RDKit EHT and xTB scalar orbital features.

The main public descriptor protocol is xTB-based and lives in the dedicated
modules exported by :mod:`eled`.  ``CalcEHTDescriptors`` and the EHT
branches of ``CalcOrbitalDescriptors`` only normalize RDKit's existing
YAeHMOP interface for comparison; they do not modify or replace YAeHMOP.
"""

from __future__ import annotations

from typing import Any, Mapping

from .xtb import XTBResult, CalcXTBDescriptors
from .calculator import XTBCalculator


def _frontier_from_fermi(energies: tuple[float, ...], fermi: float) -> tuple[float | None, float | None]:
    occupied = [energy for energy in energies if energy <= fermi + 1.0e-7]
    virtual = [energy for energy in energies if energy > fermi + 1.0e-7]
    return (occupied[-1] if occupied else None, virtual[0] if virtual else None)


def CalcEHTDescriptors(result: Any) -> dict[str, float | int | str | None]:
    """Extract stable scalar descriptors from an RDKit ``EHTResults`` object."""

    energies = tuple(float(value) for value in result.GetOrbitalEnergies())
    # YAeHMOP's Fermi value can differ from the returned HOMO by a few ulps.
    # The electron count gives an unambiguous frontier index for its
    # valence-only, doubly occupied ordering.
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
        "eht_total_energy_ev": float(result.totalEnergy),
        "eht_num_orbitals": int(result.numOrbitals),
        "eht_num_electrons": int(result.numElectrons),
    }


def CalcOrbitalDescriptors(
    mol: Any,
    *,
    backend: str = "xtb",
    conf_id: int = -1,
    keep_matrices: bool = False,
    xtb_backend: str = "native",
    frontier_depth: int = 2,
    occupation_threshold: float = 1e-7,
    **xtb_options: Any,
) -> Mapping[str, float | int | str | None]:
    """Calculate method-prefixed scalar orbital descriptors.

    ``backend`` can be ``"eht"``, ``"xtb"``, or ``"both"``. xTB options are
    forwarded to :class:`XTBCalculator`; EHT always uses RDKit's YAeHMOP interface.
    frontier_depth and occupation_threshold configure xTB feature extraction
    only; they do not change the calculation or EHT frontier convention.
    """

    normalized = backend.strip().lower()
    if normalized not in {"eht", "xtb", "both"}:
        raise ValueError("backend must be 'eht', 'xtb', or 'both'")
    descriptors: dict[str, float | int | str | None] = {}
    if normalized in {"eht", "both"}:
        from rdkit.Chem import rdEHTTools

        ok, result = rdEHTTools.RunMol(
            mol,
            confId=conf_id,
            keepOverlapAndHamiltonianMatrices=keep_matrices,
        )
        if not ok:
            raise RuntimeError("RDKit YAeHMOP EHT calculation failed")
        descriptors.update(CalcEHTDescriptors(result))
    if normalized in {"xtb", "both"}:
        result: XTBResult = XTBCalculator(
            conf_id=conf_id,
            backend=xtb_backend,
            **xtb_options,
        ).calculate(mol)
        descriptors.update(CalcXTBDescriptors(result, frontier_depth=frontier_depth,
                                                 occupation_threshold=occupation_threshold))
    return descriptors
