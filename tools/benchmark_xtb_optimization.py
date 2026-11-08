#!/usr/bin/env python3
"""Compare :func:`eled.optimize_xtb` with the xTB CLI optimizer.

The two calculations start from the same RDKit-generated coordinates and use
the same GFN method, charge, spin, accuracy and SCF iteration limit.  The
comparison is a numerical regression check, not a claim that two optimizers
must visit the same path.  Final geometries are compared after optimal rigid
alignment because xTB may choose a different translation or rotation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile
import time

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem

from eled import optimize_xtb
from eled.input import _effective_charge_uhf, _geometry_hash, _xyz_from_mol
from eled.runtime import _find_xtb_library


DEFAULT_SMILES = ("O", "CCO", "c1ccccc1", "CC(=O)O")
_METHOD_NUMBERS = {"GFN0-xTB": "0", "GFN1-xTB": "1", "GFN2-xTB": "2"}
_ENERGY_RE = re.compile(r"\benergy:\s*([-+0-9.eE]+)")
_GNORM_RE = re.compile(r"\bgnorm:\s*([-+0-9.eE]+)")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _aligned_rmsd(reference: np.ndarray, candidate: np.ndarray) -> float:
    """Return Kabsch RMSD after translation and rotation alignment."""
    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    if reference.shape != candidate.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("coordinate arrays must have shape (n_atoms, 3)")
    ref = reference - reference.mean(axis=0)
    fit = candidate - candidate.mean(axis=0)
    covariance = fit.T @ ref
    left, _, right_t = np.linalg.svd(covariance)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    aligned = fit @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - ref) ** 2, axis=1))))


def _read_xyz(path: Path) -> tuple[list[str], np.ndarray, str]:
    lines = path.read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"invalid XYZ output: {path}")
    try:
        natoms = int(lines[0].strip())
    except ValueError as exc:
        raise ValueError(f"invalid atom count in {path}") from exc
    if len(lines) < natoms + 2:
        raise ValueError(f"truncated XYZ output: {path}")
    symbols, coordinates = [], []
    for line in lines[2:2 + natoms]:
        fields = line.split()
        if len(fields) < 4:
            raise ValueError(f"invalid XYZ atom row in {path}: {line!r}")
        symbols.append(fields[0])
        coordinates.append([float(value) for value in fields[1:4]])
    array = np.asarray(coordinates, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError(f"nonfinite XYZ coordinates in {path}")
    return symbols, array, lines[1]


def _cli_optimize(mol, executable: Path, method: str, charge: int, uhf: int,
                  cycles: int, accuracy: float, workdir: Path) -> dict:
    input_path = workdir / "input.xyz"
    input_path.write_text(_xyz_from_mol(mol, -1))
    command = [str(executable), str(input_path), "--gfn", _METHOD_NUMBERS[method],
               "--chrg", str(charge), "--uhf", str(uhf), "--acc", str(accuracy),
               "--opt", "normal", "--cycles", str(cycles), "--json"]
    started = time.perf_counter()
    completed = subprocess.run(command, cwd=workdir, capture_output=True, text=True)
    elapsed = time.perf_counter() - started
    if completed.returncode:
        raise RuntimeError(
            f"xTB CLI optimization failed with exit code {completed.returncode}: "
            f"{completed.stderr[-2000:]}")
    optimized_path = workdir / "xtbopt.xyz"
    json_path = workdir / "xtbout.json"
    if not optimized_path.is_file() or not json_path.is_file():
        raise RuntimeError("xTB CLI did not produce xtbopt.xyz and xtbout.json")
    symbols, coordinates, comment = _read_xyz(optimized_path)
    match_energy = _ENERGY_RE.search(comment)
    match_gradient = _GNORM_RE.search(comment)
    if not match_energy or not match_gradient:
        raise ValueError("xTB CLI optimized XYZ comment lacks energy or gnorm")
    summary = json.loads(json_path.read_text())
    return {
        "command": command,
        "elapsed_seconds": elapsed,
        "symbols": symbols,
        "coordinates_angstrom": coordinates.tolist(),
        "energy_hartree": float(match_energy.group(1)),
        "gradient_norm_cli": float(match_gradient.group(1)),
        "json_energy_hartree": float(summary["total energy"]),
        "version": summary.get("xtb version"),
    }


def _prepare(smiles: str, seed: int) -> Chem.Mol:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError(f"RDKit embedding failed: {smiles}")
    _xyz_from_mol(mol, -1)
    return mol


def run_case(smiles: str, *, executable: Path, library: str | None, method: str,
             seed: int, cycles: int, scf_max_iterations: int, accuracy: float,
             workdir: Path) -> dict:
    mol = _prepare(smiles, seed)
    charge, uhf = _effective_charge_uhf(mol, None, None)
    initial_geometry_sha256 = _geometry_hash(_xyz_from_mol(mol, -1))
    python_started = time.perf_counter()
    optimized = optimize_xtb(
        mol, method=method, charge=charge, uhf=uhf, xtb_library=library,
        max_iterations=cycles, scf_max_iterations=scf_max_iterations,
        accuracy=accuracy,
    )
    python_elapsed = time.perf_counter() - python_started
    cli = _cli_optimize(mol, executable, method, charge, uhf, cycles, accuracy, workdir)
    python_coordinates = np.asarray(optimized.molecule.GetConformer().GetPositions(), dtype=float)
    cli_coordinates = np.asarray(cli["coordinates_angstrom"], dtype=float)
    if len(cli_coordinates) != mol.GetNumAtoms():
        raise ValueError(f"CLI atom count changed for {smiles}")
    report = optimized.report
    return {
        "smiles": smiles,
        "seed": seed,
        "method": method,
        "charge": charge,
        "uhf": uhf,
        "n_atoms": mol.GetNumAtoms(),
        "initial_geometry_sha256": initial_geometry_sha256,
        "python": {
            "success": optimized.success,
            "reason": report["reason"],
            "n_iterations": report["n_iterations"],
            "energy_hartree": report["final_energy_hartree"],
            "gradient_max_hartree_per_bohr": report["final_max_gradient_hartree_per_bohr"],
            "elapsed_seconds": python_elapsed,
        },
        "cli": cli,
        "comparison": {
            "energy_difference_hartree": float(report["final_energy_hartree"] - cli["energy_hartree"]),
            "aligned_rmsd_angstrom": _aligned_rmsd(python_coordinates, cli_coordinates),
            "initial_energy_hartree": report["initial_energy_hartree"],
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--smiles", action="append", dest="smiles", default=None,
                        help="SMILES case; repeat this option to add cases")
    parser.add_argument("--method", choices=tuple(_METHOD_NUMBERS), default="GFN2-xTB")
    parser.add_argument("--library", type=Path, default=None,
                        help="extended libxtb; defaults to normal runtime discovery")
    parser.add_argument("--executable", type=Path, default=None,
                        help="xTB executable; defaults to the sibling of the selected library")
    parser.add_argument("--seed", type=int, default=17)
    parser.add_argument("--cycles", type=int, default=100)
    parser.add_argument("--scf-max-iterations", type=int, default=250)
    parser.add_argument("--accuracy", type=float, default=1.0)
    parser.add_argument("--check", action="store_true",
                        help="fail if Python does not converge or final results differ beyond thresholds")
    parser.add_argument("--max-abs-energy-difference", type=float, default=1.0e-5,
                        help="maximum Python-minus-CLI energy difference in Hartree")
    parser.add_argument("--max-aligned-rmsd", type=float, default=5.0e-2,
                        help="maximum rigid-aligned final-coordinate RMSD in Angstrom")
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError("output directory must be empty")
    if args.cycles < 1 or args.scf_max_iterations < 1:
        raise ValueError("cycles and scf-max-iterations must be positive")
    if args.max_abs_energy_difference <= 0 or args.max_aligned_rmsd <= 0:
        raise ValueError("comparison thresholds must be positive")
    library = str(_find_xtb_library(args.library))
    library_path = Path(library).resolve()
    executable = args.executable
    if executable is None:
        executable = library_path.parent.parent / "bin" / "xtb"
    executable = executable.resolve()
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise ValueError(f"xTB executable is not available: {executable}")
    cases = args.smiles or list(DEFAULT_SMILES)
    args.output.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    rows = []
    with tempfile.TemporaryDirectory(prefix="eled-optimization-", dir=args.output) as temporary:
        workroot = Path(temporary)
        for index, smiles in enumerate(cases):
            case_dir = workroot / str(index)
            case_dir.mkdir()
            rows.append(run_case(
                smiles, executable=executable, library=library_path, method=args.method,
                seed=args.seed + index, cycles=args.cycles,
                scf_max_iterations=args.scf_max_iterations, accuracy=args.accuracy,
                workdir=case_dir,
            ))
    summary = {
        "schema_version": 1,
        "protocol": "xtb-optimization-cli-comparison-v1",
        "rdkit_version": rdBase.rdkitVersion,
        "library": str(library_path),
        "library_sha256": _sha256(library_path),
        "executable": str(executable),
        "method": args.method,
        "cycles": args.cycles,
        "scf_max_iterations": args.scf_max_iterations,
        "accuracy": args.accuracy,
        "checks": {
            "enabled": args.check,
            "max_abs_energy_difference_hartree": args.max_abs_energy_difference,
            "max_aligned_rmsd_angstrom": args.max_aligned_rmsd,
        },
        "cases": rows,
        "n_cases": len(rows),
        "n_python_converged": sum(row["python"]["success"] for row in rows),
        "max_abs_energy_difference_hartree": max(
            abs(row["comparison"]["energy_difference_hartree"]) for row in rows),
        "max_aligned_rmsd_angstrom": max(row["comparison"]["aligned_rmsd_angstrom"] for row in rows),
        "elapsed_seconds": time.perf_counter() - started,
    }
    failures = []
    for row in rows:
        if not row["python"]["success"]:
            failures.append({"smiles": row["smiles"], "reason": "python_not_converged",
                             "detail": row["python"]["reason"]})
        if abs(row["comparison"]["energy_difference_hartree"]) > args.max_abs_energy_difference:
            failures.append({"smiles": row["smiles"], "reason": "energy_difference",
                             "value": row["comparison"]["energy_difference_hartree"]})
        if row["comparison"]["aligned_rmsd_angstrom"] > args.max_aligned_rmsd:
            failures.append({"smiles": row["smiles"], "reason": "coordinate_rmsd",
                             "value": row["comparison"]["aligned_rmsd_angstrom"]})
    summary["check_failures"] = failures
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print(json.dumps(summary, indent=2))
    if args.check and failures:
        raise SystemExit(f"optimization comparison failed: {len(failures)} check(s)")


if __name__ == "__main__":
    main()
