#!/usr/bin/env python3
"""Reproduce compiled rdXTBTools versus tblite comparisons on frozen conformers.

Run from a source checkout with the compiled extension and optional tblite.
Nonzero exit means calculation errors or differences outside diagnostic limits;
the report is still written in full. Timings are single observations, not a
performance benchmark (tblite also saves integrals, and engine warm-up differs).
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import time

import numpy as np
from rdkit import Chem, rdBase

from eled import rdXTBTools
from eled.backends.tblite import run_xtb_tblite
from eled.orbital_alignment import compare_orbitals, libxtb_to_tblite_rows


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASES = ROOT / "tests/data/orbital_alignment_cases.json"
# These are diagnostic limits, not guarantees for all molecules or versions.
LIMITS = {
    "energy_eh": 2e-7,
    "orbital_energies_ev": 1.5e-4,
    "occupations": 1e-6,
    "charges_e": 5e-5,
    "dipole_e_bohr": 2e-4,
    "wiberg": 5e-5,
    "reference_orthonormality_error": 1e-8,
    "candidate_orthonormality_error": 1e-8,
    "max_projector_distance": 1e-3,
    "density_orthonormal_max_abs": 1e-4,
    "density_trace_difference": 1e-6,
}


def molecule(case):
    mol = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
    numbers = [a.GetAtomicNum() for a in mol.GetAtoms()]
    if numbers != case["atomic_numbers"]:
        raise ValueError("fixture atom order differs from this RDKit version")
    positions = np.asarray(case["positions_angstrom"], dtype=float)
    if positions.shape != (len(numbers), 3) or not np.isfinite(positions).all():
        raise ValueError("invalid fixture coordinates")
    conf = Chem.Conformer(len(numbers))
    for i, xyz in enumerate(positions):
        conf.SetAtomPosition(i, xyz.tolist())
    mol.AddConformer(conf, assignId=True)
    if (Chem.GetFormalCharge(mol) != case["charge"]
            or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())):
        raise ValueError("fixtures must have the declared charge and be closed-shell")
    return mol


def difference(a, b):
    a, b = np.asarray(a), np.asarray(b)
    if a.shape != b.shape or not np.isfinite(a).all() or not np.isfinite(b).all():
        raise ValueError("incompatible or non-finite comparison arrays")
    return float(np.max(abs(a - b)))


def compare_case(case, method, library, degeneracy_ev, accuracy=1.0):
    mol = molecule(case)
    begin = time.perf_counter()
    ok, native = rdXTBTools.RunMol(mol, method=method, xtbLibrary=library,
                                  keepOrbitalCoefficients=True, accuracy=accuracy)
    native_seconds = time.perf_counter() - begin
    if not ok:
        raise RuntimeError(f"rdXTBTools: {native.errorMessage}")
    begin = time.perf_counter()
    other = run_xtb_tblite(mol, method=method, save_integrals=True, accuracy=accuracy)
    tblite_seconds = time.perf_counter() - begin
    if other.version != "0.8.0":
        raise ValueError("this diagnostic convention was reviewed for tblite 0.8.0")
    maps = other.raw_json
    rows, phases, labels = libxtb_to_tblite_rows(
        case["atomic_numbers"], native.parameterFile, maps["orbital_map"],
        maps["shell_map"], maps["angular_momenta"], native.numOrbitals)
    cn = native.GetOrbitalCoefficients()
    ct = np.asarray(other.orbital_coefficients)
    aligned = cn[rows, :] * phases[:, None]
    orbital = compare_orbitals(
        ct, aligned, maps["overlap_matrix"], other.orbital_energies_ev,
        native.GetOrbitalEnergies(), other.orbital_occupations,
        native.GetOrbitalOccupancies(), degeneracy_ev=degeneracy_ev)
    wbo_n = native.GetWibergBondOrders()
    wbo_t = np.asarray(other.wiberg_bond_orders)
    coordinates_bohr = np.asarray(case["positions_angstrom"]) * 1.8897261254578281
    distances = np.linalg.norm(coordinates_bohr[:, None, :] - coordinates_bohr[None, :, :], axis=2)
    wbo_delta = abs(wbo_n - wbo_t)
    wi, wj = np.unravel_index(wbo_delta.argmax(), wbo_delta.shape)
    differences = {
        "energy_eh": difference(native.totalEnergyHartree, other.total_energy_eh),
        "orbital_energies_ev": difference(native.GetOrbitalEnergies(), other.orbital_energies_ev),
        "occupations": difference(native.GetOrbitalOccupancies(), other.orbital_occupations),
        "charges_e": difference(native.GetAtomicCharges(), other.partial_charges),
        "dipole_e_bohr": difference(native.GetDipole(), other.dipole_au),
        "wiberg": difference(native.GetWibergBondOrders(), other.wiberg_bond_orders),
        "wiberg_within_10_bohr": float(wbo_delta[distances < 10].max()),
        "raw_coefficients_max_abs": difference(cn, ct),
        **{key: value for key, value in orbital.items()
           if isinstance(value, float) and key not in {"degeneracy_ev", "occupation_tolerance"}},
    }
    exceeded = [key for key, limit in LIMITS.items() if differences[key] > limit]
    return {
        "status": "outside_limits" if exceeded else "within_limits",
        "exceeded": exceeded,
        "differences": differences,
        "timings_seconds": {"rdXTBTools": native_seconds, "tblite_with_overlap": tblite_seconds},
        "engines": {
            "rdXTBTools": {"library": native.libraryPath, "sha256": native.librarySHA256,
                           "api_version": native.apiVersion, "source_commit": native.sourceCommit,
                           "parameter_file": native.parameterFile,
                           "parameter_sha256": native.parameterSHA256},
            "tblite": other.metadata["engine"],
        },
        "ao_mapping": {"libxtb_rows_in_tblite_order": rows.tolist(),
                       "row_phases": phases.tolist(), "labels": labels},
        "orbital_groups": orbital["groups"],
        "largest_wiberg_difference": {"atoms": [int(wi), int(wj)],
                                      "distance_bohr": float(distances[wi, wj]),
                                      "rdXTBTools": float(wbo_n[wi, wj]), "tblite": float(wbo_t[wi, wj])},
        "values": {
            "rdXTBTools": {"energy_eh": native.totalEnergyHartree,
                           "orbital_energies_ev": native.GetOrbitalEnergies().tolist(),
                           "occupations": native.GetOrbitalOccupancies().tolist(),
                           "charges_e": native.GetAtomicCharges().tolist(),
                           "dipole_e_bohr": native.GetDipole().tolist()},
            "tblite": {"energy_eh": other.total_energy_eh,
                       "orbital_energies_ev": other.orbital_energies_ev,
                       "occupations": other.orbital_occupations,
                       "charges_e": other.partial_charges, "dipole_e_bohr": other.dipole_au},
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cases", type=Path, default=DEFAULT_CASES)
    parser.add_argument("--output", type=Path, default=ROOT / "build/orbital-comparison")
    parser.add_argument("--library", default="")
    parser.add_argument("--method", choices=["GFN1-xTB", "GFN2-xTB"], action="append")
    parser.add_argument("--case", action="append", help="fixture name; repeat to select several")
    parser.add_argument("--degeneracy-ev", type=float, default=1e-3)
    parser.add_argument("--accuracy", type=float, default=1.0)
    args = parser.parse_args()
    data = args.cases.read_bytes()
    dataset = json.loads(data)
    cases = dataset["cases"]
    if args.case:
        unknown = set(args.case) - {c["name"] for c in cases}
        if unknown:
            parser.error(f"unknown case names: {sorted(unknown)}")
        cases = [c for c in cases if c["name"] in args.case]
    if not np.isfinite(args.degeneracy_ev) or args.degeneracy_ev < 0:
        parser.error("--degeneracy-ev must be finite and nonnegative")
    if not np.isfinite(args.accuracy) or args.accuracy <= 0:
        parser.error("--accuracy must be finite and positive")
    report = {
        "schema_version": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
        "implementation_sha256": {
            name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
            for name in ("tools/compare_orbitals.py", "src/eled/orbital_alignment.py",
                         "src/eled/backends/tblite.py", "cpp/rdXTBTools/XTBTools.cpp",
                         "cpp/rdXTBTools/Wrap.cpp")},
        "dataset": {"path": str(args.cases.resolve()), "sha256": hashlib.sha256(data).hexdigest(),
                    "metadata": {k: v for k, v in dataset.items() if k != "cases"}},
        "environment": {"python": platform.python_version(), "platform": platform.platform(),
                        "rdkit": rdBase.rdkitVersion, "numpy": np.__version__,
                        "tblite": importlib.metadata.version("tblite"),
                        "thread_environment": {k: os.environ.get(k) for k in
                                               ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")}},
        "settings": {"accuracy": args.accuracy, "max_iterations": 250, "temperature_kelvin": 300.0,
                     "uhf": 0, "solvent": None, "optimize": False,
                     "degeneracy_ev": args.degeneracy_ev, "occupation_tolerance": 1e-6,
                     "ao_convention": "libxtb-6.7.1-to-tblite-0.8.0-spd"},
        "limits": LIMITS, "cases": [],
    }
    for case in cases:
        for method in args.method or ["GFN1-xTB", "GFN2-xTB"]:
            entry = {"name": case["name"], "smiles": case["smiles"], "method": method,
                     "input": case}
            try:
                entry.update(compare_case(case, method, args.library, args.degeneracy_ev, args.accuracy))
            except Exception as exc:
                entry.update(status="error", error=f"{type(exc).__name__}: {exc}")
            report["cases"].append(entry)
            print(case["name"], method, entry["status"], entry.get("exceeded", entry.get("error", "")), flush=True)
    summary = {}
    for method in sorted({c["method"] for c in report["cases"]}):
        entries = [c for c in report["cases"] if c["method"] == method]
        measured = [c for c in entries if "differences" in c]
        keys = measured[0]["differences"] if measured else []
        summary[method] = {
            "counts": {status: sum(c["status"] == status for c in entries)
                       for status in ("within_limits", "outside_limits", "error")},
            "maxima": {key: max(c["differences"][key] for c in measured) for key in keys},
        }
    report["summary"] = summary
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output / "report.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    compact = {**report, "full_report_sha256": hashlib.sha256((args.output / "report.json").read_bytes()).hexdigest(),
               "cases": [{**{k: v for k, v in c.items()
                              if k not in {"input", "ao_mapping", "values", "orbital_groups"}},
                          "mixed_groups": [g for g in c.get("orbital_groups", []) if g["size"] > 1]}
                         for c in report["cases"]]}
    (args.output / "summary.json").write_text(json.dumps(compact, indent=2, allow_nan=False) + "\n")
    fields = ["name", "method", "status", "exceeded", "error", *LIMITS,
              "raw_coefficients_max_abs", "sign_aligned_max_abs", "subspace_aligned_max_abs"]
    with (args.output / "report.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for c in report["cases"]:
            writer.writerow({**c, **c.get("differences", {}), "exceeded": ";".join(c.get("exceeded", []))})
    print(json.dumps(summary, indent=2))
    return int(any(c["status"] != "within_limits" for c in report["cases"]))


if __name__ == "__main__":
    raise SystemExit(main())
