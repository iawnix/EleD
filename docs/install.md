# Installing rdkit-eth 0.4

The validated installation target is **Linux x86-64, Python 3.11**. There are
two profiles. Installing upstream RDKit alone does not add this package or
its optional `rdXTBTools` module. The commands below run from this source
checkout; no PyPI or conda-forge release of rdkit-eth is assumed.

## Python adapter with official xTB

This profile installs RDKit, the official xTB shared library and
its parameter files through conda, then installs this Python package:

```bash
conda create --yes --prefix "$PWD/.runtime-env" --file locks/linux-64/runtime.txt
conda activate "$PWD/.runtime-env"
python -m pip install --no-deps --no-build-isolation .
python -m rdkit_eth doctor --check native
```

Use `run_xtb(mol, method="GFN2-xTB")`; native is the default backend.
GFN1/2 work with the official library. GFN0 requires the verified patched
library described below. Orbital population analysis also requires the
extension below. This profile does not compile `rdXTBTools`. The conda xTB package includes an executable,
but this adapter does not launch it for calculations.

## C++/Python rdXTBTools with verified source xTB

This is the full development/research profile, including native GFN0 and same-calculation basis export. First
create the separate toolchain and RDKit SDK environments:

```bash
conda create --yes --prefix "$PWD/.xtb-build-env" --file locks/linux-64/xtb-toolchain.txt
conda create --yes --prefix "$PWD/.rdxtb-build-env" --file locks/linux-64/rdxtb-sdk.txt
conda run --no-capture-output --prefix "$PWD/.xtb-build-env" \
  python tools/build_patched_xtb.py "$PWD/.xtb-native" \
  --dependency-prefix "$PWD/.xtb-build-env" --jobs 4
```

The builder fetches the pinned source, applies the reviewed D4 and basis-export patches, runs
upstream tests, 12 method-switch checks and a native basis/lifecycle diagnostic, and writes the
binary capability manifest. It stops if verification fails. Then install the
compiled adapter using the matching SDK:

```bash
conda activate "$PWD/.rdxtb-build-env"
source .xtb-native/activate.sh
RDKIT_ETH_BUILD_RDXTB=1 RDXTB_XTB_PREFIX="$PWD/.xtb-native" \
  python -m pip install --no-deps --no-build-isolation .
python -m rdkit_eth doctor --check compiled
python -m rdkit_eth doctor --check compiled --method GFN0-xTB --keep-overlap
```

Keep `.rdxtb-build-env`, `.xtb-build-env` and `.xtb-native` in place: the local
binary uses their native libraries. Rebuild after moving prefixes or changing
the RDKit/Boost SDK. A pip RDKit wheel is not a substitute for the matching
C++ SDK. The compiled wheel is a **local environment-bound binary**, not a
portable manylinux distribution. macOS/Windows builds remain unsupported.

The SDK lock deliberately contains no official xTB package. Do not install
the current source's mctc-lib requirements into an environment containing the
official xTB package. To compile against the official release instead, create
the runtime profile separately and set `RDXTB_XTB_PREFIX="$PWD/.runtime-env"`.

## Native orbital populations

There is no tblite backend or installation extra in 0.4. The extended source
build sets `WITH_TBLITE=OFF`; it has no tblite native linkage. The official
conda xTB distribution may independently depend on a native tblite library.
That packaging choice does not add a Python backend to this adapter.

After selecting `.xtb-native/activate.sh`, both interfaces can retain basis data:

```python
from rdkit_eth import run_xtb, calc_orbital_populations, rdXTBTools
result = run_xtb(mol, keep_overlap=True)
populations = calc_orbital_populations(result)
ok, compiled = rdXTBTools.RunMol(mol, keepOverlap=True)
if not ok:
    raise RuntimeError(compiled.errorMessage)
populations = calc_orbital_populations(compiled)
```

## First calculation

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_xtb

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
result = run_xtb(mol, method="GFN2-xTB", backend="native")
print(result.total_energy_eh)       # Hartree
print(result.orbital_energies_ev)   # eV

# With the compiled profile:
from rdkit_eth import rdXTBTools
ok, result = rdXTBTools.RunMol(mol, method="GFN2-xTB")
if not ok:
    raise RuntimeError(result.errorMessage)
print(result.GetAtomicCharges())   # elementary charge, RDKit atom order
```

## Diagnose an installation

`python -m rdkit_eth doctor` (or `rdkit-eth doctor`) prints the interpreter,
package location, package versions and optional extension file location.
Inventory mode does not import optional engines and does not prove they work.
`--check native` or `--check compiled` runs a
fixed water single point. An unavailable/broken requested backend exits with
status 1; it does not silently select a different backend. `--json` preserves
the actual engine identity, parameter/library hashes where available, and
the method/energy units for installation reports.

Use `--library /absolute/path/libxtb.so` to diagnose a specific native engine. Set `XTBPATH` to its parameter directory if
necessary. Native GFN0 capability is determined by the reviewed manifest and
binary hashes, not by the string `6.7.1` or the C API version alone.

## Build and verify distributable files locally

In the locked SDK, install the release tooling and build both variants:

```bash
python -m pip install --no-deps --require-hashes -r locks/linux-64/release-tools-py311.txt
python tools/build_release.py --compiled --xtb-prefix "$PWD/.xtb-native" \
  --output build/release-0.4.0
```

Choose a new/empty output directory. This creates a source archive, a pure
Python wheel, a compiled wheel, and `manifest.json` containing SHA256 hashes,
build package versions and lock identities. Both wheels are built from
separate clean extractions of the same source archive. Omit `--compiled` to
build only the source archive and pure wheel in the runtime environment.
The script performs no upload.

Verify the compiled artifact using the SDK and selected engine:

```bash
python tools/check_wheel.py build/release-0.4.0/compiled/*.whl \
  --kind compiled --check compiled --check native \
  --output build/release-0.4.0/installed-compiled.json
```

Verify the pure wheel in the runtime environment without additional Python engines:

```bash
conda run --no-capture-output --prefix "$PWD/.runtime-env" \
  python tools/check_wheel.py build/release-0.4.0/pure/*.whl \
  --kind pure --check native \
  --output build/release-0.4.0/installed-pure.json
```

Each check creates a fresh venv using that environment's third-party
dependencies, installs the wheel, runs `pip check`, and executes from a
neutral directory with isolated Python. It asserts that the imported package
comes from the venv, preventing an editable checkout from masking a broken
wheel. The CI workflows repeat these checks with official and patched xTB;
the patched job also checks native GFN0 and the standalone C++ interface.

For an end-user install from a generated wheel, activate its matching locked
environment and run `python -m pip install --no-deps /absolute/path/to.whl`,
then run the appropriate `doctor --check` command.

The [descriptor interfaces](descriptors.md) are available after installation.
`doctor --check` also verifies scalar/atom/fragment extraction; the
`--keep-overlap` option retains native overlap and verifies population conservation.
