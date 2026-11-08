#!/usr/bin/env python3
"""Install a wheel in a fresh venv and check it outside the source checkout.

The venv inherits only the active environment's third-party dependencies.
Package origin assertions prevent an editable checkout from satisfying the test.
"""

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import venv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("wheel", type=Path)
    parser.add_argument("--kind", required=True, choices=("pure", "compiled"))
    parser.add_argument("--check", action="append", choices=("native", "compiled"), default=[])
    parser.add_argument("--method", default="GFN2-xTB", choices=("GFN0-xTB", "GFN1-xTB", "GFN2-xTB"))
    parser.add_argument("--library")
    parser.add_argument("--keep-overlap", action="store_true")
    parser.add_argument("--reactivity", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.reactivity and (not args.check or args.method == 'GFN0-xTB'):
        parser.error('--reactivity requires --check and GFN1/2')
    wheel = args.wheel.resolve(strict=True)
    if args.library:
        args.library = str(Path(args.library).resolve())
    env = os.environ.copy()
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    with tempfile.TemporaryDirectory(prefix="eled-wheel-") as temporary:
        directory = Path(temporary)
        venv.EnvBuilder(with_pip=True, system_site_packages=True).create(directory / "venv")
        python = directory / "venv/bin/python"
        subprocess.run([str(python), "-I", "-m", "pip", "install", "--no-index", "--no-deps", "--ignore-installed",
                        str(wheel)], cwd=directory, env=env, check=True)
        subprocess.run([str(python), "-I", "-m", "pip", "check"], cwd=directory, env=env, check=True)

        def doctor(extra):
            completed = subprocess.run([str(python), "-I", "-m", "eled", "doctor", "--json", *extra],
                                       cwd=directory, env=env, capture_output=True, text=True)
            report = json.loads(completed.stdout)
            if completed.returncode:
                raise RuntimeError(json.dumps(report))
            if not Path(report["package_path"]).is_relative_to(directory / "venv"):
                raise RuntimeError("test imported a checkout instead of the installed wheel")
            return report

        inventory = doctor([])
        if bool(inventory["compiled_module_path"]) != (args.kind == "compiled"):
            raise RuntimeError("installed wheel does not have the requested module variant")
        # Check the installed console script too, from the neutral directory.
        subprocess.run([str(directory / "venv/bin/eled"), "doctor"], cwd=directory, env=env, check=True)
        checks = []
        for backend in args.check:
            extra = ["--check", backend, "--method", args.method]
            if args.keep_overlap:
                extra.append("--keep-overlap")
            if args.reactivity:
                extra.append("--reactivity")
            if args.library:
                extra.extend(["--library", args.library])
            checks.append(doctor(extra))
        result = {"wheel": str(wheel), "kind": args.kind, "inventory": inventory, "checks": checks,
                  "verification": "fresh venv, neutral working directory, isolated Python, package origin asserted"}
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Installed wheel verified: {wheel}")


if __name__ == "__main__":
    main()
