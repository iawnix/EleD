# Independent GFN0 state-contamination verification

> Upstream update (2026-11-04): issue #1486 was fixed by merged PR #1490.
> The 0.4 source lock predates that repair; 0.5 uses the upstream fix. See [fix status and migration implications](upstream/xtb-gfn0-upstream-fix.md).

On 2026-11-01, the issue was reproduced with both the conda-forge xTB 6.7.1
library and independently compiled upstream source. A small D4 source patch
eliminated the observed GFN0-to-GFN2 contamination in both source builds.

This extends the [root-cause investigation](gfn0-native-state.md). Numerical
records, source/parameter hashes, and upstream test outcomes are preserved in
[the evidence JSON](../tools/diagnostics/evidence/gfn0-state-2026-11-01.json).

## Protocol

The [standalone C program](../tools/diagnostics/xtb_gfn0_state.c) uses fixed
geometries for ethanol and water. No RDKit or Python adapter is involved.
These tests do not use its optional `--inspect-restore` mode: they neither
inspect nor modify internal Fortran symbols.

Each sequence runs in a fresh, single-threaded process:

1. Create new handles, calculate GFN2, collect energy and orbital energies,
   and delete every result/calculator/molecule/environment handle.
2. With fresh handles, load GFN0, GFN1, or GFN2 parameters; optionally perform
   the intermediate single-point calculation; delete all handles again.
3. With another set of fresh handles, calculate GFN2 on the identical geometry
   and compare against step 1.

For each source build, the matrix contains two molecules × three intermediate
methods × two intermediate actions, giving 12 sequences before and 12 after
the patch. `--expect-clean` returns exit code 3 if the energy changes by more
than `1e-10 Hartree` or any orbital energy by more than `1e-8 eV`.

The environment fixes `OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, and
`MKL_NUM_THREADS=1`; `LD_LIBRARY_PATH`, `LD_PRELOAD`, `XTBPATH`, and `XTBHOME`
are unset. Each build uses its own source-tree parameter files. Runtime
linkage was checked for the `main` reproducer. The same reproducer executable
is used before and after rebuilding its linked library.

## Results

The following are GFN2 differences after **loading GFN0 parameters only**;
completing the intermediate GFN0 calculation gives identical differences.

| Library | Molecule | Δ total energy (Hartree) | Maximum absolute orbital-energy change (eV) |
|---|---|---:|---:|
| conda-forge 6.7.1 | Ethanol | -1.111703594478541e-4 | 7.635684228234822e-4 |
| Source-built v6.7.1 | Ethanol | -1.111703594478541e-4 | 7.635684228234822e-4 |
| Source-built v6.7.1 | Water | +4.229524108723126e-7 | 1.084464535799272e-5 |
| Source-built main `a9040512` | Ethanol | -1.111703626488492e-4 | 7.635450311748853e-4 |
| Source-built main `a9040512` | Water | +4.229524188659184e-7 | 1.084465219467583e-5 |
| Patched v6.7.1 and patched main | Both | 0 | 0 |

The conda-forge ethanol check was repeated in three fresh processes for GFN0
and three for the GFN1 control; results were identical between repeats.
In each unpatched source build, all four GFN0 sequences failed the invariance
check, while all eight GFN1/GFN2 controls passed. After patching, all 12
sequences passed, with computed differences exactly zero in this environment.
The patch also preserved each build's fresh GFN2 baseline and every recorded
intermediate GFN0/GFN1/GFN2 single-point energy.

Small differences between the absolute energies of `main` and `v6.7.1` are
not used as a bug criterion. Comparisons are within each build and geometry.

## Source intervention and upstream tests

The [candidate patch](../tools/diagnostics/xtb_d4_local_references.patch) changes
only `newD4Model` in `src/disp/dftd4.F90`. It removes writes to the module
reference arrays `refq` and `refh`, selects the method's charge tables into
model/local storage, and makes secondary reference charges local as well.
The original GFN2 reference tables therefore survive later GFN0 loading.

The exact same patch applies to both versions. This before/after intervention
supports the identified D4-table mutation as the cause of the measured
sequential contamination.

| Upstream CTest suite | Before patch | After patch |
|---|---:|---:|
| main `a9040512` | 45/45 passed | 45/45 passed |
| v6.7.1 | 34/38 passed | 34/38 passed |

The release build has pre-existing failures in `xtb/gfn1`, `xtb/gfn2`,
`xtb/hessian`, and `xtb/oniom` under this toolchain. Each fails with the same
expected and actual assertion values before and after the patch; these are
preserved in the evidence JSON. Their underlying cause was not investigated,
so this is not a claim that the release build passes its entire test suite.

This is a candidate upstream repair, not an assertion of complete libxtb
thread safety. Existing-calculator interleavings, concurrent calls, other
platforms, and other releases were not validated by this experiment. The
package initially retained GFN0 CLI isolation for unverified libraries; a later
[verified-build capability](xtb-patched-build.md) enabled native GFN0 for the
local patched library. Version 0.2 removes the CLI backend and rejects
unverified native GFN0. The installed conda-forge library was not patched.
The findings were submitted upstream as
[issue #1486](https://github.com/grimme-lab/xtb/issues/1486). The complete
[submitted report](upstream/xtb-gfn0-state-issue.md) is preserved locally;
no PR has been submitted.

## Build and reproduction details

- Platform: Linux x86-64.
- Sources: upstream tag `v6.7.1` and commit
  [`a9040512dd7c15793472aa20f9647d36a266f7cd`](https://github.com/grimme-lab/xtb/commit/a9040512dd7c15793472aa20f9647d36a266f7cd).
- Isolated toolchain: conda-forge GCC/GFortran 16.2.0, CMake 4.4.5,
  Ninja, mctc-lib 0.7.0, test-drive 0.6.1, OpenBLAS 0.3.34.
- Release builds with shared libraries enabled; optional tblite and CPCM-X
  disabled. `WITH_TESTS=ON` for main; the tag builds tests by default.
- CMake compatibility flag: `CMAKE_POLICY_VERSION_MINIMUM=3.5`.
  The release additionally uses `CMAKE_C_FLAGS=-std=c17`, because GCC 16's
  default C mode rejects the old source's function named `true`.
- C executable linkage uses
  `-Wl,-rpath-link,/tmp/xtb-state-build-env/lib` to resolve the isolated
  Fortran runtime dependencies.

The following reproduces the C test against the local source build. Set
`xtb_source` and `xtb_build` to the desired source/build pair when repeating
elsewhere. Exit 3 is the expected failure for an unpatched GFN0 sequence.

```bash
xtb_source=/tmp/xtb-state-verify-0cck0304/xtb-a9040512dd7c15793472aa20f9647d36a266f7cd
xtb_build=/tmp/xtb-state-verify-0cck0304/build-main

cc -std=c11 -Wall -Wextra -Werror \
  tools/diagnostics/xtb_gfn0_state.c \
  -I"$xtb_source/include" -L"$xtb_build" \
  -Wl,-rpath,"$xtb_build" \
  -Wl,-rpath-link,/tmp/xtb-state-build-env/lib \
  -lxtb -ldl -lm -o /tmp/xtb-gfn0-state-check

env -u LD_LIBRARY_PATH -u LD_PRELOAD -u XTBPATH -u XTBHOME \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  /tmp/xtb-gfn0-state-check "$xtb_source" 0 --expect-clean
```

Add `--water` for the second geometry, `--run-intermediate` for a full
intermediate calculation, and replace `0` with `1` or `2` for controls.
The temporary source/build directories now contain the patched version.
To repeat the full before/after experiment, begin with pristine upstream
source, record the baseline, apply the patch with `patch -p1`, rebuild the
library, then run the same commands again.
