"""Result model shared by xTB backends. Units are explicit in field names."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping, Sequence

def _frontier_energies(
    energies: Sequence[float], occupations: Sequence[float], threshold: float = 1.0e-7
) -> tuple[float | None, float | None]:
    """Return HOMO and LUMO from xTB occupations.

    xTB can produce fractional occupations.  Any occupation above the small
    threshold is treated as occupied. At finite electronic temperature this
    convention can differ from xTB's reported HOMO/LUMO assignment.
    """

    if not energies or len(energies) != len(occupations):
        return None, None
    occupied = [energy for energy, occupation in zip(energies, occupations) if occupation > threshold]
    virtual = [energy for energy, occupation in zip(energies, occupations) if occupation <= threshold]
    return (occupied[-1] if occupied else None, virtual[0] if virtual else None)


@dataclass(frozen=True)
class XTBResult:
    """Normalized output from one gas-phase xTB single point.

    Energies in ``*_ev`` are electron-volts. ``total_energy_eh`` remains in
    Hartree because that is xTB's native total-energy unit. Coordinates use
    Angstrom, matching RDKit and XYZ files. Native orbital coefficients have
    axes [AO][MO]. With keep_overlap=True, the verified native extension
    supplies the same-calculation overlap and zero-based AO/shell/atom maps.
    ``version`` is the engine/package version when known, never the C ABI
    version (which is recorded separately in metadata).
    """

    method: str
    version: str | None
    total_energy_eh: float | None
    homo_lumo_gap_ev: float | None
    number_of_molecular_orbitals: int | None
    number_of_electrons: int | None
    number_of_unpaired_electrons: int | None
    orbital_energies_ev: tuple[float, ...] = ()
    orbital_occupations: tuple[float, ...] = ()
    partial_charges: tuple[float, ...] = ()
    dipole_au: tuple[float, ...] = ()
    wiberg_bond_orders: tuple[tuple[float, ...], ...] = ()
    raw_json: Mapping[str, Any] = field(default_factory=dict, repr=False)
    orbital_coefficients: tuple[tuple[float, ...], ...] = ()
    backend: str = "native"
    metadata: Mapping[str, Any] = field(default_factory=dict, repr=False)

    # Same-calculation native basis snapshot (opt-in).
    overlap_matrix: tuple[tuple[float, ...], ...] = ()
    ao_atom_indices: tuple[int, ...] = ()
    ao_shell_indices: tuple[int, ...] = ()
    shell_atom_indices: tuple[int, ...] = ()
    shell_angular_momenta: tuple[int, ...] = ()

    @property
    def homo_ev(self) -> float | None:
        homo, _ = _frontier_energies(self.orbital_energies_ev, self.orbital_occupations)
        return homo

    @property
    def lumo_ev(self) -> float | None:
        _, lumo = _frontier_energies(self.orbital_energies_ev, self.orbital_occupations)
        return lumo

    @property
    def gap_ev(self) -> float | None:
        homo, lumo = self.homo_ev, self.lumo_ev
        if homo is not None and lumo is not None:
            return lumo - homo
        return self.homo_lumo_gap_ev

    @property
    def fermi_energy_ev(self) -> float | None:
        """Use the occupied frontier energy as the discrete Fermi estimate."""

        return self.homo_ev
