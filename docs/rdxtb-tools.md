# Compiled rdXTBTools interface

`rdkit_eth.rdXTBTools` is an optional Boost.Python extension backed by a C++
library. It accepts RDKit's native `ROMol`, calls the xTB C API directly, and
returns `(success, XTBResults)` in the style of `rdEHTTools.RunMol`. The C++
namespace is `RDKit::XTBTools`. It is independent of YAeHMOP and is not part of
the upstream RDKit distribution.

## Build and install

The supported development platform is Linux x86-64 with Python 3.11 and the
matched conda-forge RDKit 2026.09.1 / Boost 1.90 SDK. From the project root:

For the complete locked installation, use [install.md](install.md). The SDK
and xTB source toolchain use separate prefixes. First build/select xTB as
described there, then:

```bash
conda create --yes --prefix "$PWD/.rdxtb-build-env" --file locks/linux-64/rdxtb-sdk.txt
conda activate "$PWD/.rdxtb-build-env"
RDKIT_ETH_BUILD_RDXTB=1 RDXTB_XTB_PREFIX="$PWD/.xtb-native" python -m pip install --no-deps --no-build-isolation -e .
```

The environment installs RDKit, its C++ development headers/libraries, Boost,
NumPy, OpenSSL, a C++ compiler, CMake and Ninja. xTB is selected from its
separate prefix. An ordinary pip install
without `RDKIT_ETH_BUILD_RDXTB=1` remains a pure Python adapter install. A pip
RDKit wheel alone does not supply this C++ SDK. Use the same environment for
building and importing the extension; mixing RDKit or Boost builds can break
the C++ ABI even when their version numbers agree. Import checks the RDKit
version, but cannot certify all ABI combinations.

Always pass `RDXTB_XTB_PREFIX` to select the intended xTB installation.
Override with `RDXTB_XTB_PREFIX=/prefix/with/xtb`; the prefix must provide
`include/xtb.h` and `lib/libxtb.so`. `RDXTB_RDKIT_PREFIX` overrides the SDK
prefix. Runtime library selection, in order, is the `xtbLibrary` argument,
`XTB_LIBRARY`, the build's xTB path, then the system loader's `libxtb.so` if
no build path was configured. An explicit path that fails to load is an error.
The module supports the gas-phase calls shared by xTB C API major versions 1
and 2. This ABI version is not the xTB program version.

For native GFN0, first follow the [patched-build guide](xtb-patched-build.md),
then select its installation before launching Python:

```bash
source .xtb-native/activate.sh
```

The C++ core verifies the build manifest, reviewed patch digest, and selected
binary digest **before loading GFN0 parameters**. Official unpatched 6.7.1
can run GFN1/2; GFN0 raises `ValueError`. `rdXTBTools` always uses the native
library. The separate `run_xtb()` API defaults to ctypes/libxtb and also
requires a verified library for GFN0; neither interface falls back to a process.

Build a wheel in the same activated environment with:

```bash
RDKIT_ETH_BUILD_RDXTB=1 python -m pip wheel --no-deps --no-build-isolation -w dist .
```

These are local development binaries tied to the SDK and xTB installations;
they do not bundle RDKit, Boost, libxtb or its Fortran/numerical dependencies.
Keep those prefixes installed. The wheel is not a portable manylinux release.
macOS and Windows need platform loader support and have not been implemented.

## Python API and units

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdEHTTools
from rdkit_eth import rdXTBTools

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
ok, xtb = rdXTBTools.RunMol(mol, method="GFN2-xTB")
if not ok:
    raise RuntimeError(xtb.errorMessage)
print(xtb.totalEnergyHartree)
print(xtb.GetOrbitalEnergies())
print(xtb.GetWibergBondOrders())

eht_ok, eht = rdEHTTools.RunMol(mol)
if eht_ok:
    print(eht.totalEnergy)  # eV, EHT's model-specific energy
```

| Quantity | `rdXTBTools.XTBResults` | EHT distinction |
| --- | --- | --- |
| Total energy | `totalEnergyHartree`, Hartree | EHT `totalEnergy` is eV; definitions and reference zeros also differ |
| Orbital energies | `GetOrbitalEnergies()`, eV, `[MO]` | Both use eV, but different Hamiltonians and bases |
| Occupations | `GetOrbitalOccupancies()`, electrons, `[MO]` | xTB occupations can be fractional at finite electronic temperature |
| Partial charges | `GetAtomicCharges()`, elementary charge, `[atom]` | xTB's charge model; do not equate it with EHT's population analysis |
| Dipole | `GetDipole()`, e·Bohr, `[x,y,z]` | Origin-dependent for charged molecules; uses supplied coordinates |
| Bond orders | `GetWibergBondOrders()`, dimensionless, `[atom,atom]` | Not the EHT overlap population matrix |
| MO coefficients | `GetOrbitalCoefficients()`, `[AO,MO]` | Optional; xTB basis ordering, not EHT basis ordering |
| Overlap / maps | `GetOverlapMatrix()`, `GetAOAtomMap()`, `GetAOShellMap()`, `GetShellAtomMap()`, `GetShellAngularMomenta()` | Optional native basis data; not EHT matrices |
| Modeled electron count | `numElectrons`, rounded sum of occupations | Valence/model electrons; water has 8 |
| Total electron count | `numTotalElectrons`, sum of nuclear charges minus charge | Includes core electrons; water has 10 |

There is no EHT-compatible `totalEnergy`, `fermiEnergy`, Hamiltonian or
overlap-population alias. Converting Hartree to eV does not make total energies
equivalent. Orbital gaps are not excitation energies. With the verified native
extension, `keepOverlap=True` retains same-calculation S and C plus zero-based
AO/shell/atom maps. This enables the shared population helpers. Orbital signs
and degenerate orbital rotations remain arbitrary.

The call signature is:

```python
RunMol(mol, confId=-1, method="GFN2-xTB", charge=None,
       numUnpairedElectrons=None, keepOrbitalCoefficients=False,
       xtbLibrary="", parameterDirectory="", maxIterations=250,
       accuracy=1.0, electronicTemperature=300.0, keepOverlap=False)
```

`confId` is an RDKit conformer **ID**, not its position in the conformer list;
`-1` selects the first conformer. Input coordinates are Angstrom and must be
finite and marked as 3D. Use `Chem.AddHs` to expand all hydrogen counts into
atoms. Charge and unpaired electrons default to RDKit formal charges and
radical counts, with electron-count/parity validation. Electronic temperature
is in kelvin; accuracy is xTB's dimensionless numerical accuracy setting.
The scope is gas-phase single points, with no geometry optimization or solvent
option. These operations are outside the current package scope.

Parameter discovery checks `parameterDirectory` first, then colon-separated
`XTBPATH`, then the selected library's `share/xtb`. An explicit directory must
contain the requested parameter file. GFN1/2 can use compiled-in parameters
when no file is found; GFN0 needs a parameter file. Metadata records the actual
library path/hash, parameter path/hash (empty for built-in parameters), C API
version, selected method/conformer, charge/spin, and GFN0 verification/source
commit when available. Hashes identify files, not all transitive dependencies.
Restart the process after replacing a loaded library.

Invalid input/options raise `ValueError`; library and ABI errors raise
`RuntimeError`. xTB setup, SCF, or result-extraction errors return `False` with
`result.errorMessage`. Numerical getters on a failed result raise
`RuntimeError`. Missing optional coefficients raise `ValueError` with
instructions to rerun with `keepOrbitalCoefficients=True`.

Results own their storage. Each NumPy getter returns a copy; scalar metadata
is read-only in Python. Atoms retain the input molecule's ordering and the
input molecule is unchanged. The wrapper snapshots `ROMol` while holding the
GIL and releases the GIL during the calculation. C++ calls are serialized;
Python calls also share the existing ctypes adapter's lock. Independent xTB
bindings and separately linked copies of the static C++ library do not share
that lock. This does not establish general libxtb thread safety.

The low-level `XTBResults` object is distinct from the `XTBResult` dataclass.
Both work with the shared descriptor functions and `annotate_mol`. Batch/cache
APIs continue to produce the dataclass through ctypes.

## Shared descriptor extraction

```python
from rdkit_eth import calc_xtb_descriptors, calc_orbital_populations, annotate_mol
ok, xtb = rdXTBTools.RunMol(mol, keepOverlap=True)
if not ok:
    raise RuntimeError(xtb.errorMessage)
features = calc_xtb_descriptors(xtb)
populations = calc_orbital_populations(xtb)
annotated = annotate_mol(mol, xtb)
```

These helpers never rerun a calculation. The wrapper records provenance from
the molecule snapshot used by C++; mapping rejects subsequently changed atom
order, topology or geometry. Standalone C++ results leave `inputProvenanceJSON`
empty. `basisExportVerified` reports the verified ABI/manifest capability;
`overlapRetained` reports whether this result actually holds basis data.
Official libxtb still handles basic GFN1/2; requesting `keepOverlap=True` without
the reviewed extension raises `ValueError`. See [definitions](descriptors.md).

## Standalone C++ use

The C++ core needs no Python runtime. Configure and install it separately:

```bash
cmake -S . -B build/core -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DRDXTB_XTB_PREFIX="$PWD/.xtb-native" \
  -DRDXTB_BUILD_PYTHON=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/build/rdxtb-sdk"
cmake --build build/core --parallel 4
ctest --test-dir build/core --output-on-failure
cmake --install build/core
```

Use `-DRDXTB_XTB_PREFIX="$CONDA_PREFIX"` for the environment's release xTB.
An external CMake project can use:

```cmake
find_package(RDKitXTB CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE RDKitXTB::RDKitXTB)
```

Configure that project with both the SDK install prefix and conda environment
in `CMAKE_PREFIX_PATH`. Its C++ code can call:

```cpp
#include <XTBTools.h>

RDKit::XTBTools::XTBOptions options;
options.method = "GFN2-xTB";
RDKit::XTBTools::XTBResults result;
if (!RDKit::XTBTools::runMol(mol, result, options)) {
  throw std::runtime_error(result.errorMessage);
}
// result.totalEnergyHartree; result.orbitalEnergiesEV
```

Invalid input throws `std::invalid_argument`; library problems throw
`std::runtime_error`. Check the boolean before reading C++ numerical fields;
unlike the Python getters, public C++ fields do not guard failed results.
Matrices are owned row-major vectors, with axes documented in `XTBTools.h`.
Runtime loading uses `dlopen`; there is no link-time dependency on libxtb.

## Validation

```bash
source .xtb-native/activate.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 RDKIT_ETH_REQUIRE_RDXTB=1 \
  python -m pytest -q
python examples/rdxtb_tools.py
```

The required-module flag prevents missing extensions from silently skipping
the compiled tests. A present extension that fails to import always fails
collection. Tests compare all three methods on ethanol, water, benzene,
ammonium and methyl radical; verify units, coefficient layout, array ownership,
input validation, SCF failure/recovery, EHT coexistence, and mixed native
threads; and reject unverified GFN0 before initialization. The C++ test also
exercises multiple calls without Python. The
[compiled CI workflow](../.github/workflows/rdxtb-tools.yml) covers release
and patched dependencies. See [local validation evidence](validation.md).
