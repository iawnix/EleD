"""Grimme Fractional Occupation Number Weighted Electron Density (FOD) descriptors.

Computes static-correlation / multireference diagnostics following Grimme &
Hansen (*Angew. Chem. Int. Ed.* 2015, 54, 12308) and matching the exact
``fermismear`` + ``fodenmak`` conventions of ``xtb --fod``:

- Separate Fermi-Dirac smearing of alpha and beta spin channels at electronic
  temperature ``T_el`` (default ``5000.0 K`` in :func:`RunFOD`):
  ``f_i^sigma = 1 / (exp((eps_i - E_F^sigma) / (k_B T_el)) + 1)``,
  with ``sum_i f_i^alpha = (N_e + N_unpaired) / 2`` and
  ``sum_i f_i^beta = (N_e - N_unpaired) / 2``.
- Grimme FOD orbital occupation weight:
  ``f_FOD,i^sigma = 1 - f_i^sigma`` for ``eps_i <= E_F^sigma`` and
  ``f_FOD,i^sigma = f_i^sigma`` for ``eps_i > E_F^sigma``,
  summed over spins as ``f_FOD,i = f_FOD,i^alpha + f_FOD,i^beta``.
- Molecular static-correlation index ``n_fod_e = sum_i f_FOD,i`` (in ``e``).
- Atomic FOD populations ``fod_pop_e`` (and ``s/p/d`` shell contributions) via
  Lowdin (default, matching ``xtb``) or Mulliken orbital populations when the
  same-calculation overlap basis is retained (``keepOverlap=True``).
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from numbers import Real
from typing import Any

import numpy as np

from .calculator import _calculate
from .electronic_descriptors import (
    _as_result,
    _has_population_basis,
    _orbital_basis,
    _orbitals,
)
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .populations import CalcOrbitalPopulations

# Match xtb's Boltzmann constant (kB in Hartree/K * autoev).
_KB_HARTREE_PER_K = 3.166808578545117e-6
_HARTREE_TO_EV_XTB = 27.21138505
BOLTZMANN_EV_PER_K = _KB_HARTREE_PER_K * _HARTREE_TO_EV_XTB

FOD_PROTOCOL = 'eled-fod-v1'

__all__ = [
    'BOLTZMANN_EV_PER_K',
    'CalcFODFromResult',
    'RunFOD',
]


def _validate_temperature(value: Any, name: str = 'electronicTemperature') -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real) or not math.isfinite(value) or value < 0.0:
        raise ValueError(f'{name} must be a finite nonnegative real number (kelvin)')
    return float(value)


def _validate_scheme(scheme: str) -> str:
    if scheme not in {'lowdin', 'mulliken'}:
        raise ValueError("scheme must be 'lowdin' or 'mulliken'")
    return scheme


def _fermi_channel(energies_ev: np.ndarray, n_spin: int, temperature_k: float) -> tuple[float, np.ndarray, np.ndarray]:
    """Solve Fermi-Dirac smearing and Grimme FOD weights for one spin channel."""
    nmo = int(len(energies_ev))
    if n_spin < 0 or n_spin > nmo:
        raise ValueError(f'spin electron count {n_spin} is out of range for {nmo} orbitals')
    if n_spin == 0:
        e_fermi = float(energies_ev[0] - 10.0 * max(BOLTZMANN_EV_PER_K * temperature_k, 1e-4))
        return e_fermi, np.zeros(nmo, dtype=float), np.zeros(nmo, dtype=float)
    if n_spin == nmo:
        e_fermi = float(energies_ev[-1] + 10.0 * max(BOLTZMANN_EV_PER_K * temperature_k, 1e-4))
        return e_fermi, np.ones(nmo, dtype=float), np.zeros(nmo, dtype=float)

    kt = BOLTZMANN_EV_PER_K * temperature_k
    if kt <= 1e-12:
        e_fermi = 0.5 * float(energies_ev[n_spin - 1] + energies_ev[n_spin])
        occ = np.zeros(nmo, dtype=float)
        occ[:n_spin] = 1.0
        return e_fermi, occ, np.zeros(nmo, dtype=float)

    def _occ(ef: float) -> np.ndarray:
        x = np.clip((energies_ev - ef) / kt, -500.0, 500.0)
        return 1.0 / (np.exp(x) + 1.0)

    span = max(20.0 * kt, 1.0)
    lo = float(energies_ev[0]) - span
    hi = float(energies_ev[-1]) + span
    while float(np.sum(_occ(lo))) > n_spin:
        lo -= span
        span *= 2.0
    span = max(20.0 * kt, 1.0)
    while float(np.sum(_occ(hi))) < n_spin:
        hi += span
        span *= 2.0

    ef = 0.5 * (lo + hi)
    for _ in range(120):
        ef = 0.5 * (lo + hi)
        f = _occ(ef)
        diff = float(np.sum(f)) - float(n_spin)
        if abs(diff) <= 1e-13 or (hi - lo) <= 1e-14 * max(1.0, abs(ef)):
            break
        if diff < 0.0:
            lo = ef
        else:
            hi = ef

    occ = _occ(ef)
    # Grimme fodenmak: below/at Fermi level -> hole (1 - f), above Fermi level -> particle (f).
    fod = np.where(energies_ev <= ef, 1.0 - occ, occ)
    fod = np.clip(fod, 0.0, 1.0)
    return float(ef), occ, fod


def _ao_orbital_weights(result, scheme: str) -> tuple[np.ndarray, float]:
    """Return [AO, MO] population weights matching Lowdin or Mulliken scheme."""
    basis = _orbital_basis(result)
    c, sc = basis.coefficients, basis.overlap @ basis.coefficients
    ao_mo = (basis.overlap_half @ c) ** 2 if scheme == 'lowdin' else c * sc
    return ao_mo / np.einsum('ij,ij->j', c, sc)[None, :], basis.orthonormality_error


def CalcFODFromResult(
    result: Any,
    *,
    scheme: str = 'lowdin',
    electronicTemperature: float | None = None,
) -> dict[str, Any]:
    """Compute Grimme FOD descriptors from an existing xTB single-point result.

    Never invokes xTB. When ``electronicTemperature`` is omitted, uses the
    electronic temperature recorded in ``result.metadata['calculation_settings']``.
    Atomic and ``s/p/d`` shell FOD populations are computed when same-calculation
    overlap and basis maps are retained (``keepOverlap=True``); otherwise
    ``atoms`` is empty and ``missing['atoms']`` records the reason rather than
    filling zeros.
    """
    res = _as_result(result)
    scheme = _validate_scheme(scheme)
    energies, occupations, _, _ = _orbitals(res)
    nmo = int(len(energies))
    if nmo == 0:
        raise ValueError('FOD analysis requires molecular orbital energies')

    calc_settings = res.metadata.get('calculation_settings') or {}
    recorded_temp = calc_settings.get('electronic_temperature')
    if electronicTemperature is None:
        if recorded_temp is None:
            raise ValueError(
                'electronicTemperature must be specified when calculation_settings '
                'does not record electronic_temperature'
            )
        temp_k = _validate_temperature(recorded_temp, 'recorded electronic_temperature')
    else:
        temp_k = _validate_temperature(electronicTemperature, 'electronicTemperature')

    self_consistent_temp = (
        recorded_temp is not None and abs(float(recorded_temp) - temp_k) <= 1e-6
    )

    n_e = res.numElectrons if res.numElectrons is not None else int(round(float(np.sum(occupations))))
    uhf = res.numUnpairedElectrons
    if uhf is None:
        uhf = res.metadata.get('uhf')
    if uhf is None:
        raise ValueError('FOD analysis requires numUnpairedElectrons spin provenance')
    n_e = int(n_e)
    uhf = int(uhf)
    if n_e < 0 or uhf < 0 or uhf > n_e or (n_e - uhf) % 2 != 0:
        raise ValueError(f'inconsistent electron count ({n_e}) and unpaired electrons ({uhf})')

    n_alpha = (n_e + uhf) // 2
    n_beta = (n_e - uhf) // 2
    ef_alpha, occ_alpha, fod_alpha = _fermi_channel(energies, n_alpha, temp_k)
    ef_beta, occ_beta, fod_beta = _fermi_channel(energies, n_beta, temp_k)

    fod_orb = fod_alpha + fod_beta
    n_fod_e = float(np.sum(fod_orb))
    n_fod_below_e = float(
        np.sum(np.where(energies <= ef_alpha, fod_alpha, 0.0))
        + np.sum(np.where(energies <= ef_beta, fod_beta, 0.0))
    )
    n_fod_above_e = float(
        np.sum(np.where(energies > ef_alpha, fod_alpha, 0.0))
        + np.sum(np.where(energies > ef_beta, fod_beta, 0.0))
    )

    missing: dict[str, str] = {}
    atoms_rows: list[dict[str, Any]] = []
    ortho_err: float | None = None
    atomic_available = False
    shell_available = False
    fod_sum_error: float | None = None

    if _has_population_basis(res):
        pop = CalcOrbitalPopulations(res, scheme=scheme)
        atom_pop = np.asarray(pop['atom_pop'], dtype=float)
        natoms = int(atom_pop.shape[0])
        atom_fod = atom_pop @ fod_orb
        ortho_err = float(pop['orthonormality_error'])
        fod_sum_error = float(np.sum(atom_fod) - n_fod_e)
        atomic_available = True

        has_shells = bool(res.aoShellMap) and bool(res.shellAngularMomenta)
        if has_shells:
            ao_mo, _ = _ao_orbital_weights(res, scheme)
            ao_fod = ao_mo @ fod_orb
            ao_atom = np.asarray(res.aoAtomMap, dtype=int)
            ao_shell = np.asarray(res.aoShellMap, dtype=int)
            shell_l = np.asarray(res.shellAngularMomenta, dtype=int)
            ao_l = shell_l[ao_shell]
            shell_available = True
            for a in range(natoms):
                mask_a = ao_atom == a
                f_s = float(np.sum(ao_fod[mask_a & (ao_l == 0)]))
                f_p = float(np.sum(ao_fod[mask_a & (ao_l == 1)]))
                f_d = float(np.sum(ao_fod[mask_a & (ao_l == 2)]))
                atoms_rows.append({
                    'atom_index': int(a),
                    'fod_pop_e': float(atom_fod[a]),
                    'fod_s_e': f_s,
                    'fod_p_e': f_p,
                    'fod_d_e': f_d,
                    'missing': {},
                })
        else:
            missing['atom_shells'] = 'shell_angular_momenta_not_retained'
            for a in range(natoms):
                atoms_rows.append({
                    'atom_index': int(a),
                    'fod_pop_e': float(atom_fod[a]),
                    'fod_s_e': None,
                    'fod_p_e': None,
                    'fod_d_e': None,
                    'missing': {
                        'fod_s_e': 'shell_angular_momenta_not_retained',
                        'fod_p_e': 'shell_angular_momenta_not_retained',
                        'fod_d_e': 'shell_angular_momenta_not_retained',
                    },
                })
    else:
        missing['atoms'] = 'overlap_basis_not_retained'

    payload = {
        'n_fod_e': n_fod_e,
        'n_fod_below_fermi_e': n_fod_below_e,
        'n_fod_above_fermi_e': n_fod_above_e,
        'fod_fermi_alpha_ev': ef_alpha,
        'fod_fermi_beta_ev': ef_beta,
        'electronic_temperature_k': temp_k,
        'scf_electronic_temperature_k': float(recorded_temp) if recorded_temp is not None else None,
        'self_consistent_temperature': self_consistent_temp,
        'orbital_fod_e': [float(v) for v in fod_orb],
        'atomic_fod_available': atomic_available,
        'shell_fod_available': shell_available,
        'fod_sum_error_e': fod_sum_error,
        'orthonormality_error': ortho_err,
        'atoms': atoms_rows,
        'scheme': scheme,
        'method': res.method,
        'backend': res.backend,
        'analysis_settings': {
            'protocol': FOD_PROTOCOL,
            'reference': 'Grimme & Hansen, Angew. Chem. Int. Ed. 2015, 54, 12308',
            'scheme': scheme,
            'electronic_temperature_k': temp_k,
            'boltzmann_ev_per_k': BOLTZMANN_EV_PER_K,
            'fod_definition': 'f_FOD,i^sigma = 1 - f_i^sigma (eps_i <= E_F^sigma) else f_i^sigma; summed over alpha and beta',
            'n_electrons': n_e,
            'num_unpaired_electrons': uhf,
            'n_alpha': n_alpha,
            'n_beta': n_beta,
        },
        'units': {
            'n_fod_e': 'e',
            'n_fod_below_fermi_e': 'e',
            'n_fod_above_fermi_e': 'e',
            'fod_fermi_alpha_ev': 'eV',
            'fod_fermi_beta_ev': 'eV',
            'electronic_temperature_k': 'K',
            'orbital_fod_e': 'e',
            'fod_pop_e': 'e',
            'fod_s_e': 'e',
            'fod_p_e': 'e',
            'fod_d_e': 'e',
            'fod_sum_error_e': 'e',
        },
        'invariance': (
            'invariant to rigid translation/rotation, AO permutation, MO phase, '
            'and degenerate orbital subspace rotation; preserves original RDKit atom ordering'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))


def RunFOD(
    mol: Any,
    *,
    electronicTemperature: float = 5000.0,
    scheme: str = 'lowdin',
    options: XTBOptions | None = None,
    interface: str = 'native',
) -> dict[str, Any]:
    """Run one self-consistent finite-temperature xTB single point and compute FOD.

    Defaults to ``electronicTemperature=5000.0 K`` and ``scheme='lowdin'`` to
    match Grimme's standard ``xtb --fod --etemp 5000`` protocol. Supports both
    ``interface='native'`` and ``interface='compiled'`` with identical outputs.
    """
    from rdkit import Chem

    temp_k = _validate_temperature(electronicTemperature, 'electronicTemperature')
    scheme = _validate_scheme(scheme)
    if interface not in ('native', 'compiled'):
        raise ValueError("interface must be 'native' or 'compiled'")
    snapshot = Chem.Mol(mol)
    base = options or XTBOptions()
    _xyz_from_mol(snapshot, base.confId)
    charge, unpaired = _effective_charge_uhf(snapshot, base.charge, base.numUnpairedElectrons)
    selected = replace(
        base,
        charge=charge,
        numUnpairedElectrons=unpaired,
        electronicTemperature=temp_k,
        keepOverlap=True,
    )
    result = _calculate(snapshot, selected, interface)

    report = CalcFODFromResult(result, scheme=scheme, electronicTemperature=temp_k)
    report['analysis_settings']['interface'] = interface
    return report
