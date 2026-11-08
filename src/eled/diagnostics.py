"""Installation diagnostics, with explicit optional calculation checks."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import platform
import sys


def installation_info():
    """Report installed packages without importing optional calculation engines."""
    import rdkit
    import eled

    versions = {}
    for name in ("eled", "rdkit", "numpy"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = None
    spec = importlib.util.find_spec("eled.rdXTBTools")
    return {
        "schema_version": 1,
        "python": sys.executable,
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "package_path": str(Path(eled.__file__).resolve()),
        "rdkit_path": str(Path(rdkit.__file__).resolve()),
        "versions": versions,
        "compiled_module_path": spec.origin if spec else None,
        "compiled_module_import_checked": False,
    }


def check_backend(backend, method, library=None, keep_overlap=False, reactivity=False):
    """Calculate fixed water, returning actual engine identity and explicit units."""
    from rdkit import Chem

    mol = Chem.AddHs(Chem.MolFromSmiles("O"))
    conf = Chem.Conformer(3)
    for index, xyz in enumerate(((0., 0., 0.), (.9572, 0., 0.), (-.239987, .927297, 0.))):
        conf.SetAtomPosition(index, xyz)
    conf.Set3D(True)
    mol.AddConformer(conf)
    if backend == "compiled":
        try:
            from . import rdXTBTools
        except ImportError as exc:
            raise RuntimeError(
                "rdXTBTools is unavailable or incompatible. Build with "
                "ELED_BUILD_RDXTB=1 in the locked RDKit SDK environment; "
                "see docs/install.md. " + str(exc)
            ) from exc
        ok, result = rdXTBTools.RunMol(mol, method=method, xtbLibrary=library or "", keepOverlap=keep_overlap)
        if not ok:
            raise RuntimeError(result.errorMessage)
        energy = result.totalEnergyHartree
        charges = result.GetAtomicCharges()
        nao = result.numOrbitals
        engine = {"library": result.libraryPath, "sha256": result.librarySHA256,
                  "api_version": result.apiVersion, "parameter_file": result.parameterFile,
                  "parameter_sha256": result.parameterSHA256,
                  "source_commit": result.sourceCommit,
                  "native_gfn0_verified": result.gfn0StateFixVerified,
                  "compiled_rdkit_version": rdXTBTools.rdkitVersion,
                  "basis_export_verified": result.basisExportVerified,
                  "overlap_retained": result.overlapRetained}
    else:
        from . import XTBCalculator
        result = XTBCalculator(method=method, backend=backend,
                                xtb_library=library,
                                keep_overlap=keep_overlap).calculate(mol)
        energy = result.total_energy_hartree
        charges = result.partial_charges
        nao = result.number_of_molecular_orbitals
        engine = result.metadata["engine"]
    if (energy is None or not math.isfinite(energy) or energy >= 0
            or nao is None or nao < 1 or len(charges) != 3
            or not all(math.isfinite(float(q)) for q in charges)
            or abs(sum(charges)) > 1e-5):
        raise RuntimeError("water single point returned invalid energy, orbital count or charges")
    from .electronic_descriptors import CalcXTBDescriptors, CalcXTBAtomDescriptors, CalcXTBFragmentDescriptors
    descriptors = CalcXTBDescriptors(result)
    atoms = CalcXTBAtomDescriptors(result)
    fragments = CalcXTBFragmentDescriptors(result, {"oxygen": [0], "hydrogens": [1, 2]})
    if len(atoms) != 3 or abs(sum(a["q_e"] for a in atoms)) > 1e-5:
        raise RuntimeError("atomic descriptors failed the charge check")
    population_check = None
    if keep_overlap:
        from .populations import CalcOrbitalPopulations, CalcFrontierDescriptors
        populations = CalcOrbitalPopulations(result)
        frontier = CalcFrontierDescriptors(mol, result)
        electrons = sum(populations["mo_n_e"])
        if abs(electrons - descriptors["xtb_n_e"]) > 1e-6:
            raise RuntimeError("orbital populations failed the electron count check")
        population_check = {"electron_sum_e": electrons,
                            "orthonormality_error": populations["orthonormality_error"],
                            "hl_ovlp": frontier["hl_ovlp"]}
    from .reports import CalcDescriptorReport
    descriptor_report = CalcDescriptorReport(mol, result, include_orbitals=keep_overlap)
    if descriptor_report['quality']['status'] == 'error':
        raise RuntimeError('descriptor report failed quality checks')
    reactivity_report = None
    if reactivity:
        from .options import XTBOptions
        from .reactivity import RunReactivity
        reactivity_report = RunReactivity(mol, interface=backend,
            options=XTBOptions(method=method, xtb_library=library),
            population_scheme='lowdin' if keep_overlap else 'charges')
    return {"backend": backend, "method": method, "molecule": "water",
            "total_energy_hartree": energy, "n_orb": nao,
            "charge_sum_e": float(sum(charges)), "engine": engine,
            "descriptors": descriptors, "fragment_descriptors": fragments,
            "population_check": population_check,
            "descriptor_report": descriptor_report, "reactivity_report": reactivity_report}


def main(argv=None):
    parser = argparse.ArgumentParser(prog="eled")
    commands = parser.add_subparsers(dest="command", required=True)
    doctor = commands.add_parser("doctor", help="show installation identity and optionally test a backend")
    doctor.add_argument("--json", action="store_true", help="emit a machine-readable report")
    doctor.add_argument("--check", choices=("native", "compiled"))
    doctor.add_argument("--method", default="GFN2-xTB", choices=("GFN0-xTB", "GFN1-xTB", "GFN2-xTB"))
    doctor.add_argument("--library", help="explicit libxtb shared library for native/compiled checks")
    doctor.add_argument("--keep-overlap", action="store_true", help="verify native basis export and population descriptors")
    doctor.add_argument("--reactivity", action="store_true", help="also verify fixed-geometry GFN1/2 charge-state descriptors")
    args = parser.parse_args(argv)
    if args.reactivity and (not args.check or args.method == 'GFN0-xTB'):
        parser.error('--reactivity requires --check and GFN1/2')
    report = installation_info()
    report["status"] = "inventory"
    if args.check:
        try:
            report["check"] = check_backend(args.check, args.method, args.library, args.keep_overlap, args.reactivity)
            report["status"] = "ok"
            if args.check == "compiled":
                report["compiled_module_import_checked"] = True
        except Exception as exc:
            report["status"] = "error"
            report["error"] = f"{type(exc).__name__}: {exc}"
    if args.json:
        print(json.dumps(report, indent=2, allow_nan=False))
    else:
        print(f"eled {report['versions']['eled'] or '(source checkout)'}")
        print(f"Python: {report['python']}")
        print(f"Package: {report['package_path']}")
        print(f"RDKit: {report['versions']['rdkit']}")
        print(f"C++ module file: {report['compiled_module_path'] or 'not installed'}")
        if report["status"] == "ok":
            print(f"{args.check} / {args.method}: OK, water energy "
                  f"{report['check']['total_energy_hartree']:.10f} Hartree")
        elif report["status"] == "error":
            print(report["error"])
        else:
            print("Inventory only. Add --check native/compiled to run a water single point.")
    return int(report["status"] == "error")
