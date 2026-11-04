# Installing rdkit-eth 0.8

The standard profile uses **one extended libxtb for GFN0/1/2 and all current
descriptors**, with both ctypes and compiled `rdXTBTools`. It targets Linux
x86-64 / Python 3.11. Installing upstream RDKit alone does not enable it.

## Standard installation

Install conda, clone this repository, then run from its root:

```bash
python tools/install_native.py --jobs 6
conda activate "$PWD/.rdkit-eth/runtime"
python -m rdkit_eth doctor --check compiled --keep-overlap
python examples/electronic_descriptors.py --populations
```

The installer creates two dedicated prefixes from the supplied exact locks:

- `.rdkit-eth/runtime`: matched RDKit/Boost/Python SDK, the extended shared
  libxtb and parameters, and the installed `rdXTBTools` wheel.
- `.rdkit-eth/toolchain`: Fortran compiler, mctc-lib and numerical dependencies.

It fetches the pinned xTB revision containing upstream PR #1490, applies only
our basis-export patch, runs upstream and installed-library regressions, builds
pure/compiled wheels from clean source-archive extractions, installs the
compiled wheel and checks all three methods through both interfaces with
population analysis. Checks run outside the checkout, without `XTB_LIBRARY`
or `XTBPATH`. Both interfaces discover the same library in the runtime prefix.

Artifacts are retained in `.rdkit-eth/artifacts/build-*/`; `installation.json`
records wheel/source identity and the six installed calculation checks. No
upload is performed. Use `--directory /absolute/path` to select another root.
Reruns accept exactly matching conda environments and reviewed engine source;
changed environments or unrelated source edits are rejected without resetting
them. Use a new directory when changing the source/toolchain lock.

Keep **both** prefixes installed at their original paths: native libraries
use their RPATHs. This is a reproducible local source installation, not a
portable binary distribution. Windows/macOS and public package-index delivery
are not implemented. A plain `pip install .` still installs only the Python
adapter; it does not invoke a Fortran build. The old conda recipe remains a
basic official-library compatibility recipe, not the standard full profile.

If an older shell has `XTB_LIBRARY`, `XTBPATH` or `XTB_PREFIX` set, unset those
explicit overrides before using the new environment. The installer clears
these overrides for its own checks. Explicit library selection remains
available for diagnostics, but ordinary use needs only conda activation.

## First calculation

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_xtb, rdXTBTools, calc_orbital_populations

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
result = run_xtb(mol, method="GFN2-xTB", keep_overlap=True)
print(result.total_energy_eh)  # Hartree
print(calc_orbital_populations(result))

ok, compiled = rdXTBTools.RunMol(mol, method="GFN0-xTB", keepOverlap=True)
if not ok:
    raise RuntimeError(compiled.errorMessage)
print(calc_orbital_populations(compiled))
```

Omit `keep_overlap` / `keepOverlap` for basic descriptors to avoid retaining
the additional overlap matrix. Both APIs retain results in original atom order.
The source build disables tblite and CPCM-X, and needs neither Python tblite
nor ASE. GFN0/1/2 and radicals use the same engine, with spin-summed populations.

## Upgrade from 0.7

Rebuild/reinstall the adapter and restart Python; libxtb is unchanged. The C++
result layout now includes `failureKind`. Reactivity error handling distinguishes
runtime faults from calculation failures. See [persistent batches](persistent-batches.md)
for disk cache/checkpoint identity, explicit retry and storage limits.

## Upgrade from 0.6

Rebuild/reinstall the adapter with the standard installer and restart Python.
The native engine is unchanged. Reactivity now defaults to 300 K and reference
wavefunction reuse. For the former defaults, explicitly pass
`XTBOptions(electronic_temperature=0.)` and `initialization="independent"`.
Protocol version 2.0 and each state's initialization are recorded in reports.

## Upgrade from 0.5

The extended engine is unchanged in 0.6. Rebuild the adapter in the matched SDK
(the C++ result layout now includes numerical settings); restart Python sessions.
The standard installer can be rerun and verifies the new reports and GFN1/2
reactivity in addition to all-method population checks. Existing 0.5 artifacts
remain historical. See [new API conventions](descriptor-reports.md).

## Upgrade from 0.4

The 0.4 `.xtb-native` installation uses an older upstream revision plus our
D4 repair. Preserve it to reproduce old experiments. The standard installer
builds the current adapter in its own directory and uses the official GFN0 repair instead.
Previously verified 0.4 libraries are still recognized when selected explicitly.
Unknown libraries do not gain GFN0 permission just because they report 6.7.1.

## Manual build and official-library compatibility

For engine-only development, see [the source-build guide](xtb-patched-build.md).
To build the adapter against an already verified installation:

```bash
RDKIT_ETH_BUILD_RDXTB=1 RDXTB_XTB_PREFIX=/path/to/extended-xtb \
  python -m pip install --no-deps --no-build-isolation .
```

Run this in the matched RDKit SDK. The separate `locks/linux-64/runtime.txt`
retains an official xTB 6.7.1 compatibility environment, supporting basic
GFN1/2 but not our basis export or verified GFN0. It is used by compatibility
CI rather than the recommended installation. Its upstream xTB package may
link native tblite; this does not apply to our standard source build.

## Diagnostics and release artifacts

`python -m rdkit_eth doctor --check native` or `--check compiled` runs a water
calculation. Add `--method GFN0-xTB --keep-overlap` to check the full capability.
Add `--reactivity` with GFN1/2 to check the three charge states.
`--json` records the actual library, parameters and descriptor checks. Missing
capabilities fail explicitly. Inventory without `--check` does not prove the
engine works. Restart Python after replacing a loaded library.

To reproduce an independent artifact directory in the standard environment:

```bash
python tools/build_release.py --compiled --xtb-prefix "$CONDA_PREFIX" \
  --output build/release-0.8.0
python tools/check_wheel.py build/release-0.8.0/compiled/*.whl \
  --kind compiled --check compiled --check native --keep-overlap \
  --output build/release-0.8.0/installed.json
```

Choose a new/empty output directory. Wheels contain the adapter, not copies of
libxtb, RDKit or compiler libraries. `check_wheel.py` verifies imports from a
fresh venv outside the checkout. See [descriptor definitions](descriptors.md)
and [validation](validation.md).
