# Local validation — 2026-11-03


## Current scope: 0.4 native basis export

All production calculations use libxtb. The reviewed native extension exports
S and AO/shell/atom maps from the same successful calculation as C; tblite is
absent from the backend, dependencies and active comparison tests. Historical
comparison evidence is preserved below and in the diagnostic archive.

Local validation on Linux x86-64 / Python 3.11 / RDKit 2026.09.1:

| Check | Result |
| --- | --- |
| Extended libxtb upstream suite | 45 passed |
| Installed D4 method-switch sequences | 12 passed |
| C basis ABI and result lifecycle diagnostic | GFN0/1/2 passed; worst water CᵀSC error 1.33e-15 |
| Full adapter suite, extended library | 205 passed, 1 official-only test skipped |
| Full adapter suite, official conda xTB 6.7.1 | 122 passed, 84 capability-dependent tests skipped |
| Standalone C++ core | Passed against extended and official engines |

The 66 GFN1/2 frozen-conformer cases all satisfy the native overlap metric,
map/shell contracts and electron conservation. This replaces cross-engine
coefficient alignment as the production population validation. Additional
GFN0/1/2 cases cover ethanol, ammonium and methyl radical, compiled/ctypes
agreement, owned NumPy outputs, provenance checks and actual benzene
subspace invariance. Build tests verify patch idempotence and preserve tampered
source rather than overwriting it. The C diagnostic includes result copying,
invalid buffers, disabled retention, failed SCF and Hessian invalidation/recovery.
CPCM-X entry-point invalidation is implemented but not numerically exercised,
since this source build disables CPCM-X.

`ldd` on `.xtb-native/lib/libxtb.so` shows no tblite dependency. The locked
runtime environment has neither Python tblite nor ASE. Official conda xTB
may have a native tblite dependency imposed by its own distribution.

A [committed validation snapshot](validation/0.4.0/README.md) retains the
local test and installed-wheel records. Version 0.4 artifacts and full reports are recorded under
`dist/0.4.0/`; see its validation manifest for artifact hashes and installation
checks. The source and pure/compiled wheels are built locally; no upload or
remote CI execution is implied. Compiled wheels need the matching local SDK
and native-library prefixes. These checks establish numerical/API contracts,
not accuracy for excited-state prediction or a new electronic Hamiltonian.

## Historical 0.3 electronic descriptors

The descriptor layer now includes molecular scalars, atomic/bond/fragment
summaries, and opt-in same-calculation tblite population analysis. See the
[API definitions and validation contracts](descriptors.md). The local suites
passed **159 tests with the patched library** and **152 tests with 7 GFN0-only
skips with the official library**. New coverage includes nonorthogonal metric
handling, charge/electron conservation, signed Mulliken values, AO permutation,
orbital phase and degenerate rotations, coordinate transforms, fragment
accounting, compiled-result compatibility and overlap-aware cache identity.

On ten real GFN1/2 method/molecule combinations (water, ethanol, ammonium,
dimethyl sulfoxide, chloromethane), occupation-weighted Mulliken populations
recover tblite's atomic charges within 1e-8 e. Symmetric benzene tests actual
degenerate frontier subspaces. This checks numerical implementation and
invariance, not predictive accuracy for excited states. No native/tblite
cross-basis overlap substitution is used by the descriptor API.

Installed-wheel diagnostics now check scalar/atom/fragment extraction, and
explicit tblite checks verify population conservation. Reports for 0.3 live
with the new local artifacts. The 0.1/0.2 artifacts and earlier numerical
records are retained separately. Remote GitHub CI remains unexecuted.

## 0.2 validation after CLI removal

Version 0.2 removes the external calculation backend and automatic fallback.
The gas-phase libxtb calculation path is now the default; tblite remains an
explicit optional backend. Tests check execution with an empty `PATH` and
subprocess launches disabled, rejection of removed backends/options, missing
library failures, GFN0 verification, charge/electron conservation, mapping,
sequential batching and cache reuse. The compiled and tblite comparisons remain.
A fixed water reference retains a numerical check against an archived 0.1 CLI
calculation without requiring an executable at runtime.

The local patched-library suite passed **118 tests**. With the official xTB
library, **111 passed and 7 native-GFN0-only cases skipped**. Installation checks for
0.2 are recorded separately from the old 0.1 artifacts. Remote GitHub CI has
not been run. See [migration notes](migration-0.2.md) and [installation](install.md).

## Historical validation records (through 0.1)

The dated records below describe the implementations at the time of each
experiment. Their references to CLI, automatic fallback, optimization and
parallel subprocess batches do not describe the current API. Historical
numerical evidence and the archived 0.1 packages have been preserved.


Environment: Linux x86-64; Python 3.11.18; RDKit 2026.09.1; conda-forge
xTB 6.7.1 (`gfortran_ha02e54d_5`); C API version 10000. The isolated test
environment is at `/tmp/rdkit-xtb-check`; no existing user environment was
modified. CLI and shared library came from the same conda package.

Run from the project root:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /tmp/rdkit-xtb-check/bin/python -m pytest -q
```

Result: **28 tests passed**, including 9 real-engine integration cases.
The tests cover:

- GFN2 CLI/native comparisons for methanol, benzene, ammonium and methyl
  radical. Absolute tolerances: 2e-6 Hartree for total energy; 2e-5 eV for
  orbital energies; 2e-5 for atomic charges and existing-bond WBOs.
- GFN1 CLI/native agreement; GFN0 execution in a separate CLI process.
- Native coefficient `[AO][MO]` ordering against an independently generated
  xTB Molden file, comparing absolute values to allow orbital sign changes.
- ALPB water optimization; automatic backend selection and missing-library
  fallback; coordinate/conformer preservation and property mapping.
- Invalid H counts, inconsistent spin/electron parity, invalid coordinates,
  wrong structure/atom order, incomplete arrays, timeout and missing outputs.
- Parallel ordering, duplicate calculation reuse, failed-calculation retry,
  engine/parameter-aware cache keys, and independent result metadata copies.

## GFN0 native limitation

A reproducible sequential experiment found that a GFN0 call through the
6.7.1 C API changes subsequent GFN2 results in that process. For an ethanol
conformer embedded with seed 1, GFN2 after GFN0 differed from GFN2 before GFN0
by approximately -1.11e-4 Hartree in total energy and up to 7.64e-4 eV in
orbital energies. GFN1 did not produce this effect in the same check.
A follow-up pure-C experiment localized the cause to D4 initialization:
GFN0 overwrites module-level `refq` and `refh` arrays, and later GFN2
initialization reuses them. Loading GFN0 parameters alone triggers the change;
restoring just these two arrays restores the baseline. See the
[root-cause analysis and standalone reproducer](gfn0-native-state.md).
The result was reported upstream in [issue #1486](https://github.com/grimme-lab/xtb/issues/1486).

For an unverified library, the adapter rejects explicit native GFN0 calls and
routes auto GFN0 requests to the CLI. Verified patched builds now support
native GFN0 as described below. All native calls through this adapter are
serialized; a mutex alone does not resolve the GFN0 state issue.

## Patched dependency — 2026-11-02

The automated [source-build workflow](xtb-patched-build.md) was exercised from
a fresh Git fetch into `.xtb-rdkit/source` and repeated incrementally. It pins
upstream commit `a9040512dd7c15793472aa20f9647d36a266f7cd`, applies the reviewed
D4 patch, and installs to `.xtb-rdkit` using the separate `.xtb-build-env`
compiler/dependency environment. Both directories are ignored by Git.

- Upstream tests: **45/45 passed**.
- Installed-library C regression: **12/12 method-switch sequences passed**,
  with zero total-energy and orbital-energy differences.
- Python adapter suite against the patched installation: **33 passed**,
  including 9 real-engine integration cases.
- The same suite against official conda-forge 6.7.1: **33 passed**.
- The example workflow completed with the patched CLI.

The source commit uses C API 2.0 (`20000`), while release 6.7.1 uses 1.0
(`10000`). The 25 functions bound by this adapter have identical signatures
in both headers. The adapter now accepts both API major versions for these
gas-phase calls and still rejects unknown major versions. The changed
`xtb_setSolvent` function is not bound; solvation continues through the CLI.

The added source-preparation tests verify repeatable patch application,
rejection of a different commit, and preservation of unrelated source edits.
Build provenance, binary/parameter hashes, licenses, and regression results
are installed under `.xtb-rdkit/share/xtb/rdkit-eth/`.

## Native GFN0 enablement — 2026-11-02

The build helper now emits a versioned `native_gfn0` capability after the
installed-library regression succeeds. The adapter checks this record, the
reviewed D4 patch hash, and the actual library digest before selecting native
GFN0. Missing/malformed records and replaced binaries leave `auto` on CLI;
explicit native requests fail before initializing GFN0, including through the
direct `run_xtb_native` entry point. GFN1/2 do not require this capability.

Validation: **51 passed** against the patched installation; **45 passed,
6 skipped** against official 6.7.1. The six skipped cases require native GFN0.
The builder again passed 45 upstream tests and 12 method-switch checks.

Additional cases cover GFN0 CLI/native agreement for ethanol, water, benzene,
ammonium, and methyl radical (energy, orbitals, occupations, charges, and
bond orders); sequential batch caching without CLI calls; and
GFN2 → GFN0 → GFN1 → GFN0 on another molecule → GFN2. GFN2's final energy,
orbital energies and coefficients agree with baseline. Native calls remain
serialized; this does not establish general libxtb thread safety.

## Compiled rdXTBTools — 2026-11-02

The independent Boost.Python extension was compiled against the matching
conda-forge RDKit/librdkit-dev 2026.09.1, Boost 1.90, Python 3.11.18 and NumPy
SDK in `.rdxtb-build-env`, using GNU C++ 16.2. The C++ core calls libxtb
directly; a test disables the Python ctypes and CLI calculation entry points
while the compiled call still succeeds.

- Complete suite with patched xTB C API 20000: **88 passed**, including
  **37 compiled-interface tests**.
- Complete suite with official 6.7.1 C API 10000: **76 passed, 12 skipped**.
  The skips require a verified native-GFN0 installation.
- Standalone C++ CTest passed against both libraries. A separate core-only
  build using the release's API 1 header also compiled and passed CTest.
- C++ SDK installation and an external CMake consumer passed. The consumer
  links the exported `RDKitXTB::RDKitXTB` target and runs without Python.
- Opt-in pip editable installation and compiled wheel build passed. The wheel
  was installed to a separate target directory and executed GFN2 → GFN0 →
  GFN1 → GFN2, with identical initial/final GFN2 energies.
- The default pure Python wheel was also built and inspected: it contains no
  compiled extension, even when one exists in the source checkout.
- A source archive containing the C++ sources, CMake configuration and xTB
  build helpers was built, then successfully rebuilt into a compiled wheel.
- The EHT/xTB comparison example completed. This preserves separate total
  energy fields and bond-population semantics.

The numerical comparisons cover GFN0/1/2 on ethanol, water, benzene, ammonium
and methyl radical. CLI tolerances are 2e-6 Hartree for total energy and 2e-5
for orbital energies in eV, occupations, charges and existing-bond WBOs.
Compiled and ctypes results agree within 1e-10 Hartree for energies; dipoles
and WBOs within 1e-9; absolute orbital coefficients within 1e-8. The existing
Molden comparison independently checks the ctypes reference's `[AO,MO]`
layout. Modeled/total electron counts are distinguished (8/10 for water).

Other checks exercise atom/conformer identity and input preservation,
finite 3D coordinates and explicit-H requirements, charge/spin validation,
copy ownership, missing coefficients, failed SCF diagnostics and recovery,
shared serialization with the ctypes adapter, and rejection of missing or
mismatched GFN0 provenance before parameter loading. The standalone C++ test
also exercises four concurrent callers and failed-calculation cleanup.

To repeat the patched Python run from the project root:

```bash
source .xtb-rdkit/activate.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  RDKIT_ETH_REQUIRE_RDXTB=1 RDKIT_ETH_EXPECT_NATIVE_GFN0=1 \
  .rdxtb-build-env/bin/python -m pytest -q
```

For the release run, start a fresh Python process and set `XTB_LIBRARY`,
`XTB_PATH`, and `XTBPATH` to the matching official installation, with
`RDKIT_ETH_EXPECT_NATIVE_GFN0=0`. The [build/API guide](rdxtb-tools.md)
documents installation, C++ consumers and platform/ABI limitations.

## Optional tblite backend — 2026-11-02

The optional `tblite==0.8.0` wheel was installed separately with CFFI in the
same Python environment. It did not install ASE. The explicit
`backend="tblite"` path was tested on GFN1/2 water, ethanol, benzene and
ammonium against the patched libxtb C API. Results agree within `2e-7` Hartree
for total energies, `1.5e-4` eV for orbital energies, `5e-5` elementary charge
for charges, `2e-4` e·Bohr for dipoles and `5e-5` for Wiberg bond orders.
A separate exploratory **GFN2** run on water, ethanol, benzene, acetic acid
and ammonium found maxima of `4e-10` Hartree, `8.4e-5` eV, `2.4e-5` e,
`7.9e-5` e·Bohr and `1.7e-5`, respectively. These were not GFN1/GFN2
benchmark-wide maxima. The broader evidence below supersedes that small probe.

The comparison also verified that `backend="auto"` does not silently select
tblite, that GFN0 and spin-polarized results produce explicit errors, and that
tblite's orbital/shell/atom map metadata is retained. A separate environment
without tblite raises `TBLiteNotFoundError` with installation instructions.
The optional backend is excluded from the base dependency set; ASE is not
needed. These tolerances apply to the tested tblite 0.8.0 and patched xTB
build, not to all future versions.

With tblite installed, the complete suite is **98 passed** against the patched
engine and **86 passed, 12 skipped** against official 6.7.1. Without the
optional wheel, tblite-specific tests are skipped; engine discovery and
compiled-extension availability can produce additional skips.

## Direct compiled-module orbital comparison — 2026-11-03

[Orbital alignment](orbital-alignment.md) documents source-reviewed s/p/d
conversion, fixed conformers, S-metric sign/subspace/density comparisons and
per-case JSON/CSV reports. It uses `rdXTBTools.RunMol` directly. At default
accuracy 1 the local build has 51/66 cases within all diagnostic limits; at
accuracy 0.01 all 66 are within the same limits. Official conda 6.7.1 has
64/66 within limits at accuracy 0.01; GFN1 bromomethane and iodomethane remain
explicit exceptions. All three runs have zero calculation errors.

These cross-engine diagnostics are separate from regression-test pass counts.
The new tests cover rotated molecular coordinates, d-basis phases, invalid
maps, nonorthogonal metrics and degenerate/equal-occupation subspaces.

Complete suite: **117 passed** with the local source build; **105 passed,
12 skipped** with the official release. In the environment without tblite,
**73 passed, 24 skipped**; those skips also include compiled/engine-specific
checks. These runs used one OpenMP/OpenBLAS thread.

The follow-up [controlled source investigation](xtb-671-basis-differences.md)
reproduced the 64/66 outcome with release v6.7.1 rebuilt using the same
toolchain as current upstream. Separate on-site S/H0 and 5s/5p source edits
localized the Br and I exceptions, respectively; together they produce 66/66.
Current upstream without our D4 patch also produces 66/66 and passes 45
upstream tests. Four raw-tag test failures under this toolchain were separately
traced to the unassigned `stmp` in a gradient routine: applying upstream #1121
alone restores 38/38 tests. A water finite-difference check confirms the fix;
the installed conda binary also passes that derivative check. The research
builds and partial basis patches are not production replacements.

## 0.1 packaging and locked installation — 2026-11-03

The [installation guide](install.md) now uses separate SHA256-locked conda
prefixes for the official runtime, RDKit SDK, and xTB source toolchain.
All three were recreated in fresh directories. The fresh source build passed
45 upstream tests and 12 installed-library method-switch checks. The Python
suite including the new installation diagnostic checks passed **119 tests**.

`tools/build_release.py` builds pure and compiled wheels from independent
extractions of the same source archive. `tools/check_wheel.py` installs each
wheel into a fresh venv, runs outside the source checkout with isolated
Python, asserts the imported package origin, and runs `pip check`. Both
variants passed installed GFN2 checks; the compiled wheel also passed native
GFN0 checks. These checks initially ran without the optional tblite Python
package, confirming that it is not needed for the libxtb paths.

The generated artifact manifest records hashes and build/lock identities;
the installed-wheel reports record the actual loaded libraries and parameter
hashes. The C++ wheel remains local to its SDK/native-library prefixes, not a
portable binary distribution. CI now includes the same source-archive and
installed-wheel checks with separate official/patched engine jobs. Workflow
configuration was checked locally; remote GitHub CI has not been run in this
workspace, and no package has been published.

## Scope

This validates numerical agreement between adapters and an installed xTB
engine, not the scientific accuracy of xTB against experiment or excited-state
reference data. Tests do not establish portability across all xTB releases or
operating systems. macOS/Windows and the conda recipe remain unverified. The
patched source build is validated on Linux x86-64 as described in
[the patched-build guide](xtb-patched-build.md). The included Linux GitHub
Actions workflow has not been run remotely. Persistent caching and AO-to-atom
mapping are not implemented.
