"""Geometry optimization driven by gradients from the extended libxtb API."""

from __future__ import annotations

from dataclasses import dataclass
import math
from numbers import Integral, Real
from typing import Any, Iterable

import numpy as np

from .backends.native import ANGSTROM_TO_BOHR, _run_xtb_native
from .result import XTBResult

HARTREE_TO_KCAL_MOL = 627.5094740631


@dataclass(frozen=True)
class XTBOptimizationResult:
    """Owned optimized molecule, convergence report and final xTB result."""

    molecule: Any
    report: dict
    result: XTBResult | None = None

    @property
    def success(self) -> bool:
        return bool(self.report.get("converged", False))


def _positive(value, name, *, integer=False):
    if isinstance(value, bool) or not isinstance(value, (Integral if integer else Real)):
        raise ValueError(f"{name} must be a positive {'integer' if integer else 'real'}")
    if not math.isfinite(float(value)) or value <= 0:
        raise ValueError(f"{name} must be a positive {'integer' if integer else 'real'}")
    return int(value) if integer else float(value)


def _fixed_indices(values: Iterable[int], natoms: int) -> tuple[int, ...]:
    if values is None or isinstance(values, (str, bytes)):
        raise ValueError("fixed_atoms must be a sequence of atom indices")
    result = []
    for value in values:
        if isinstance(value, bool) or not isinstance(value, Integral) or not 0 <= value < natoms:
            raise ValueError("fixed_atoms contains an invalid atom index")
        result.append(int(value))
    return tuple(sorted(set(result)))


def _gradient(result: XTBResult, natoms: int) -> np.ndarray:
    values = result.metadata.get("gradient_hartree_per_bohr")
    if values is None or len(values) != 3 * natoms:
        raise ValueError("the selected libxtb build did not return a gradient")
    array = np.asarray(values, dtype=float).reshape(natoms, 3)
    if not np.isfinite(array).all():
        raise ValueError("xTB returned a nonfinite gradient")
    return array


def _lbfgs_direction(gradient: np.ndarray, pairs: list[tuple[np.ndarray, np.ndarray, float]]) -> np.ndarray:
    """Return a limited-memory quasi-Newton descent direction."""
    vector = gradient.reshape(-1).copy()
    alphas = []
    for s, y, rho in reversed(pairs):
        alpha = rho * float(np.dot(s, vector))
        alphas.append(alpha)
        vector -= alpha * y
    if pairs:
        s, y, _ = pairs[-1]
        curvature = float(np.dot(s, y))
        scale = curvature / max(float(np.dot(y, y)), 1e-30)
        scale = min(max(scale, 1e-3), 1e3)
    else:
        scale = 4.0
    vector *= scale
    for (s, y, rho), alpha in zip(pairs, reversed(alphas)):
        vector += s * (alpha - rho * float(np.dot(y, vector)))
    return -vector.reshape(gradient.shape)


def optimize_xtb(
    mol: Any,
    *,
    conf_id: int = -1,
    method: str = "GFN2-xTB",
    charge: int | None = None,
    uhf: int | None = None,
    xtb_library=None,
    max_iterations: int = 200,
    scf_max_iterations: int = 250,
    gradient_tolerance: float = 1.0e-4,
    energy_tolerance: float = 1.0e-8,
    step_tolerance: float = 1.0e-5,
    max_step: float = 0.20,
    fixed_atoms: Iterable[int] = (),
    accuracy: float = 1.0,
    electronic_temperature: float = 300.0,
) -> XTBOptimizationResult:
    """Optimize one RDKit conformer with an extended-libxtb gradient.

    The input molecule is copied. Coordinates are updated with bounded
    L-BFGS steps and an energy-backtracking line search; all engine calls
    use the same reviewed extended library as :class:`XTBCalculator`. ``max_iterations``
    counts geometry steps, while ``scf_max_iterations`` is passed to xTB for
    each energy/gradient evaluation. ``max_step`` and ``step_tolerance`` are in
    Angstrom, while gradient and energy tolerances use Hartree/Bohr and Hartree.
    """
    from rdkit import Chem

    iterations = _positive(max_iterations, "max_iterations", integer=True)
    scf_iterations = _positive(scf_max_iterations, "scf_max_iterations", integer=True)
    grad_tol = _positive(gradient_tolerance, "gradient_tolerance")
    energy_tol = _positive(energy_tolerance, "energy_tolerance")
    step_tol = _positive(step_tolerance, "step_tolerance")
    step_limit = _positive(max_step, "max_step")
    working = Chem.Mol(mol)
    # Validate the conformer before copying coordinates.  The native call
    # repeats this validation, but doing it here gives callers a deterministic
    # error before allocating an optimization state.
    from .input import _xyz_from_mol
    _xyz_from_mol(working, conf_id)
    natoms = working.GetNumAtoms()
    fixed = _fixed_indices(fixed_atoms, natoms)
    if len(fixed) == natoms:
        raise ValueError("fixed_atoms cannot contain every atom")
    coordinates = np.asarray(working.GetConformer(conf_id).GetPositions(), dtype=float)
    if not np.isfinite(coordinates).all():
        raise ValueError("coordinates must be finite")

    options = dict(conf_id=conf_id, method=method, charge=charge, uhf=uhf,
                   xtb_library=xtb_library, accuracy=accuracy,
                   electronic_temperature=electronic_temperature,
                   _return_gradient=True)
    current = _run_xtb_native(working, max_iterations=scf_iterations, **options)
    initial_energy = float(current.total_energy_hartree)
    history = []
    converged = False
    reason = "max_iterations"
    final_gradient = _gradient(current, natoms)
    movable = np.ones(natoms, dtype=bool)
    movable[list(fixed)] = False
    initial_masked = final_gradient[movable]
    initial_max_gradient = (float(np.linalg.norm(initial_masked, axis=1).max())
                            if len(initial_masked) else 0.0)
    pairs: list[tuple[np.ndarray, np.ndarray, float]] = []

    for iteration in range(iterations + 1):
        masked = final_gradient[movable]
        max_gradient = float(np.linalg.norm(masked, axis=1).max()) if len(masked) else 0.0
        if max_gradient <= grad_tol:
            converged, reason = True, "gradient_tolerance"
            break
        if iteration == iterations:
            break
        direction = np.zeros_like(final_gradient)
        direction[movable] = _lbfgs_direction(final_gradient[movable], pairs)
        if float(np.sum(direction[movable] * final_gradient[movable])) >= 0.0:
            direction[movable] = -final_gradient[movable]
        norms = np.linalg.norm(direction[movable], axis=1)
        direction_norm = float(norms.max()) if len(norms) else 0.0
        if not math.isfinite(direction_norm) or direction_norm <= 0.0:
            reason = "invalid_search_direction"
            break
        scale = min(1.0, step_limit * ANGSTROM_TO_BOHR / direction_norm)
        old_gradient = final_gradient.copy()
        accepted = None
        accepted_delta = None
        accepted_step = None
        for backtrack in range(14):
            displacement_bohr = direction * (scale * (0.5 ** backtrack))
            displacement_angstrom = displacement_bohr / ANGSTROM_TO_BOHR
            trial = Chem.Mol(working)
            trial_conf = trial.GetConformer(conf_id)
            for atom_index, delta in enumerate(displacement_angstrom):
                if movable[atom_index]:
                    point = coordinates[atom_index] + delta
                    trial_conf.SetAtomPosition(atom_index, point)
            trial_result = _run_xtb_native(trial, max_iterations=scf_iterations, **options)
            delta_energy = float(current.total_energy_hartree - trial_result.total_energy_hartree)
            if delta_energy > 0.0:
                accepted, accepted_delta, accepted_step = trial, delta_energy, displacement_angstrom
                current = trial_result
                break
        if accepted is None:
            reason = "line_search_failed"
            break
        working = accepted
        coordinates = np.asarray(working.GetConformer(conf_id).GetPositions(), dtype=float)
        final_gradient = _gradient(current, natoms)
        s = accepted_step[movable].reshape(-1) * ANGSTROM_TO_BOHR
        y = (final_gradient[movable] - old_gradient[movable]).reshape(-1)
        curvature = float(np.dot(s, y))
        if curvature > 1e-12 and np.isfinite(curvature):
            pairs.append((s, y, 1.0 / curvature))
            del pairs[:-7]
        max_step_observed = float(np.linalg.norm(accepted_step[movable], axis=1).max())
        final_masked = final_gradient[movable]
        max_gradient_after = (float(np.linalg.norm(final_masked, axis=1).max())
                              if len(final_masked) else 0.0)
        history.append({"iteration": iteration + 1, "energy_hartree": float(current.total_energy_hartree),
                        "energy_change_hartree": accepted_delta,
                        "max_gradient_before_hartree_per_bohr": max_gradient,
                        "max_gradient_hartree_per_bohr": max_gradient_after,
                        "max_displacement_angstrom": max_step_observed})
        # Energy/step stagnation alone is not sufficient for a stationary
        # point.  Check the new gradient immediately so a line-search plateau
        # cannot be reported as a successful optimization.
        if (max_gradient_after <= grad_tol and accepted_delta <= energy_tol
                and max_step_observed <= step_tol):
            converged, reason = True, "gradient_energy_and_step_tolerance"
            break

    final_gradient = _gradient(current, natoms)
    masked = final_gradient[movable]
    final_max_gradient = float(np.linalg.norm(masked, axis=1).max()) if len(masked) else 0.0
    if not converged and final_max_gradient <= grad_tol:
        converged, reason = True, "gradient_tolerance"
    report = {
        "schema_version": 1,
        "protocol_version": "1.0",
        "report_type": "xtb_geometry_optimization",
        "converged": converged,
        "status": "converged" if converged else "failed",
        "reason": reason,
        "method": current.method,
        "conformer_id": int(current.metadata.get("conformer_id", conf_id)),
        "n_atoms": natoms,
        "n_iterations": len(history),
        "max_iterations": iterations,
        "scf_max_iterations": scf_iterations,
        "fixed_atom_indices": list(fixed),
        "gradient_tolerance_hartree_per_bohr": grad_tol,
        "energy_tolerance_hartree": energy_tol,
        "step_tolerance_angstrom": step_tol,
        "max_step_angstrom": step_limit,
        "initial_energy_hartree": initial_energy,
        "final_energy_hartree": float(current.total_energy_hartree),
        "final_energy_kcal_mol": float(current.total_energy_hartree * HARTREE_TO_KCAL_MOL),
        "initial_max_gradient_hartree_per_bohr": initial_max_gradient,
        "final_max_gradient_hartree_per_bohr": final_max_gradient,
        "history": history,
        "engine": current.metadata.get("engine"),
        "geometry_sha256": current.metadata.get("geometry_sha256"),
    }
    return XTBOptimizationResult(working, report, current)
