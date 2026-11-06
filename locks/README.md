# Locked environment snapshots (0.5)

These locks target **Linux x86-64 / CPython 3.11**. Conda explicit files pin
package URLs, versions, build variants and SHA256 hashes. They are snapshots
of the tested environments, not a new dependency solve on each installation.

| File | Purpose |
| --- | --- |
| `linux-64/runtime.txt` | Legacy RDKit + upstream xTB compatibility snapshot for packaging tests |
| `linux-64/rdxtb-sdk.txt` | Matched RDKit 2026.09.1 / Boost 1.90 / C++ SDK; no xTB package |
| `linux-64/xtb-toolchain.txt` | Fortran 16.2 / mctc-lib 0.7 toolchain for the pinned xTB source |
| `linux-64/release-tools-py311.txt` | Source/wheel build frontend and PEP 517 hooks |

Keep the SDK/runtime and source toolchain in separate prefixes. In particular,
the official xTB package requires an older mctc-lib series than the current
pinned source. The runtime lock may contain **native** tblite as a transitive
dependency of conda's xTB; it does not install a tblite Python API. The project-built libxtb disables
native tblite as well. Historical Python locks are archived under
`tools/diagnostics/archive/tblite-0.3/`.

The release-tools pip lock requires the matching conda environment (which supplies
NumPy and packaging). Install using `python -m pip install --require-hashes
--no-deps -r FILE`. All additional Python wheels are explicitly listed. ASE is not included.

Recreate a prefix with `conda create --prefix PREFIX --file FILE`. See
[the installation guide](../docs/install.md) for complete commands. YAML
environment files remain editable development specifications; the explicit
locks are the locked installation path.

To update a lock, first create and test the new environment, then export it
with `conda list --explicit --sha256 --prefix PREFIX`. Recreate it in a fresh
prefix and repeat source, installed-wheel, and numerical checks before
accepting the update. Do not just change a version string in an explicit file.
Pip wheel hashes must likewise come from the exact downloaded artifacts.

The upstream-fixed source revision, D4 file identity and basis-export patch are locked in
`tools/xtb-source.json`. Its installed binary manifest records the actual
library hash. The locks do not imply that compiler outputs are bit-for-bit
identical, that package hosting is permanent, or that other platforms/ABIs
have been tested.

The standard `tools/install_native.py` uses `rdxtb-sdk.txt` for its runtime
and `xtb-toolchain.txt` for its separate compiler prefix. It adds the built
extended libxtb directly to that runtime. `runtime.txt` is retained only as a
legacy upstream compatibility snapshot; it cannot run the public calculation
API because it has no project extensions.
