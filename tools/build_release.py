#!/usr/bin/env python3
"""Build source/pure-wheel artifacts, optionally a local compiled wheel.

Both wheels are built from separate extractions of the same source archive.
No artifact is uploaded. Run with the locked SDK and release-tool dependencies.
"""

import argparse
from email.parser import BytesParser
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
import zipfile


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def wheel_record(path, compiled):
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        extensions = [n for n in names if n.startswith("rdkit_eth/") and n.endswith(".so")]
        if bool(extensions) != compiled:
            raise RuntimeError(f"unexpected compiled module presence in {path}")
        if compiled and not any("rdXTBTools" in n for n in extensions):
            raise RuntimeError("compiled wheel is missing rdXTBTools")
        metadata = BytesParser().parsebytes(archive.read(next(n for n in names if n.endswith(".dist-info/METADATA"))))
        if "rdkit_eth/backends/tblite.py" in names:
            raise RuntimeError("wheel contains the removed tblite backend")
        if "rdkit_eth/backends/cli.py" in names:
            raise RuntimeError("wheel contains the removed CLI backend")
        if "rdkit_eth/__main__.py" not in names:
            raise RuntimeError("wheel is missing the installation diagnostic entry point")
        requirements = metadata.get_all("Requires-Dist", [])
        if any(r.startswith("tblite") for r in requirements):
            raise RuntimeError("tblite must not be a package dependency")
    return {"file": str(path), "sha256": digest(path), "compiled": compiled,
            "version": metadata["Version"], "requirements": requirements}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="new, empty artifact directory")
    parser.add_argument("--compiled", action="store_true")
    parser.add_argument("--xtb-prefix", type=Path, default=os.environ.get("RDXTB_XTB_PREFIX"))
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        parser.error("choose an empty output directory to preserve existing artifacts")
    if args.compiled and (args.xtb_prefix is None or not (args.xtb_prefix / "include/xtb.h").is_file()):
        parser.error("--compiled requires --xtb-prefix containing include/xtb.h")
    output.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env.pop("RDKIT_ETH_BUILD_RDXTB", None)
    subprocess.run([sys.executable, "-m", "build", "--sdist", "--no-isolation",
                    "--outdir", str(output), str(ROOT)], env=env, check=True)
    source, = output.glob("*.tar.gz")
    with tarfile.open(source) as archive:
        members = archive.getnames()
        for required in ("CMakeLists.txt", "tools/xtb-source.json", "locks/linux-64/rdxtb-sdk.txt",
                         "tests/data/orbital_alignment_cases.json"):
            if not any(n.endswith("/" + required) for n in members):
                raise RuntimeError(f"source archive is missing {required}")
        if any(n.endswith(".so") for n in members):
            raise RuntimeError("source archive contains a prebuilt extension")
    records = []
    for compiled in ([False, True] if args.compiled else [False]):
        destination = output / ("compiled" if compiled else "pure")
        destination.mkdir()
        with tempfile.TemporaryDirectory(prefix="rdkit-eth-build-") as temporary:
            with tarfile.open(source) as archive:
                archive.extractall(temporary, filter="data")
            checkout, = Path(temporary).iterdir()
            build_env = env.copy()
            if compiled:
                build_env.update(RDKIT_ETH_BUILD_RDXTB="1", RDXTB_XTB_PREFIX=str(args.xtb_prefix.resolve()),
                                 RDXTB_RDKIT_PREFIX=os.environ.get("RDXTB_RDKIT_PREFIX", sys.prefix))
            subprocess.run([sys.executable, "-m", "build", "--wheel", "--no-isolation",
                            "--outdir", str(destination), str(checkout)], env=build_env, check=True)
        wheel, = destination.glob("*.whl")
        records.append(wheel_record(wheel, compiled))
    manifest = {
        "schema_version": 1, "python": sys.executable,
        "source": {"file": str(source), "sha256": digest(source)}, "wheels": records,
        "build_packages": {name: importlib.metadata.version(name)
                           for name in ("rdkit", "numpy", "setuptools", "wheel", "build", "pyproject-hooks")},
        "lock_sha256": {str(p.relative_to(ROOT)): digest(p) for p in sorted((ROOT / "locks").rglob("*.txt"))},
        "xtb_source_lock_sha256": digest(ROOT / "tools/xtb-source.json"),
        "compiled_wheel_scope": "local Linux SDK/prefix; external libraries are not bundled",
    }
    if args.compiled:
        library = args.xtb_prefix.resolve() / "lib/libxtb.so"
        manifest["compiled_default_xtb"] = {"path": str(library), "sha256": digest(library)}
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"Artifacts and SHA256 manifest: {output}")


if __name__ == "__main__":
    main()
