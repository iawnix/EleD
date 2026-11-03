# Building the native xTB dependency (0.4)

> Upstream update (2026-11-04): issue #1486 was fixed by merged PR #1490.
> The current source lock predates that repair and still uses its verified
> local patch. See [fix status and migration implications](upstream/xtb-gfn0-upstream-fix.md).

This workflow builds the pinned upstream source with two reviewed patches: the
GFN0 D4 state fix and a project-local basis-export C ABI. The result supports
GFN0/1/2 single points and native orbital population descriptors through both
ctypes and `rdXTBTools`. It is not an upstream xTB or RDKit release.
The D4 report is [issue #1486](https://github.com/grimme-lab/xtb/issues/1486).

## Reproduce the build

From the checkout on Linux x86-64:

```bash
conda create --yes --prefix "$PWD/.xtb-build-env" --file locks/linux-64/xtb-toolchain.txt
conda run --no-capture-output --prefix "$PWD/.xtb-build-env" \
  python tools/build_patched_xtb.py "$PWD/.xtb-native" \
  --dependency-prefix "$PWD/.xtb-build-env" --jobs 4
source .xtb-native/activate.sh
```

The compiler environment is separate from the RDKit runtime/SDK. Keep both it
and `.xtb-native` installed: local binaries resolve compiler and numerical
libraries through their RPATHs. These binaries are not portable distributions.
The build disables `WITH_TBLITE` and `WITH_CPCMX`; xTB's own ALPB/GBSA remains
in the engine, while this adapter exposes only gas-phase single points.
The upstream executable is built for provenance and diagnostics; the adapter
does not launch it.

The [source lock](../tools/xtb-source.json) pins commit
`a9040512dd7c15793472aa20f9647d36a266f7cd`, the D4 patch, the basis-export patch
and all modified source files by SHA256. This upstream main commit reports
6.7.1 but differs from the release tag; see [earlier verification](gfn0-reverification.md).
The builder accepts a clean pinned checkout or an exactly matching patched
checkout. It refuses unrelated or mismatched edits and never resets them.
Reruns reuse the incremental build. An optional second positional argument
selects an existing checkout. Use a new prefix when changing the source lock
to preserve prior verified engines.

The build runs all upstream tests, 12 installed-library D4 method-switch
checks, and the native basis/lifecycle C diagnostic before writing capability
metadata. There is no skip-tests option. The installation includes:

```text
.xtb-native/
  source/                     # pinned source and visible local diff
  build/                      # upstream build and C diagnostics
  lib/libxtb.so
  include/xtb.h
  bin/xtb
  share/xtb/                  # matching GFN parameter files
  share/xtb/rdkit-eth/         # build.json, source.diff, patches, licenses, test log
  activate.sh
```

Activation selects `XTB_LIBRARY` and `XTBPATH` without changing `PATH`.
Restart Python after changing/replacing a loaded library. Build an adapter
using the matching RDKit SDK as described in [install.md](install.md).

## Basis extension ABI v1

Upstream `xtb_getAPIVersion()` is unchanged. Optional exported symbols are:

```c
int xtb_rdkit_getBasisAPIVersion(void);
void xtb_rdkit_setKeepOverlap(xtb_TEnvironment, xtb_TResults, bool);
void xtb_rdkit_getBasisSize(xtb_TEnvironment, xtb_TResults, int *sizes);
void xtb_rdkit_getOverlap(xtb_TEnvironment, xtb_TResults, int capacity, double *out);
void xtb_rdkit_getBasisMap(xtb_TEnvironment, xtb_TResults, int nao, int nshell,
                         int *ao_atom, int *ao_shell, int *shell_atom, int *shell_l);
```

Enable retention on the result before `xtb_singlepoint`. The actual overlap
and basis maps are captured in the GFN0 or GFN1/2 calculation and belong to
that result's wavefunction. S uses the same native spherical AO order as C.
The C buffer is column-major; both Python interfaces return `[AO, AO]` S and
`[AO, MO]` C. Size output has exactly three integers: natoms, nshell, nao.
Map indices are zero-based; shell angular momentum uses s=0, p=1, d=2.
For map buffers supply `nao` elements for AO maps and `nshell` for shell maps.
Overlap capacity must be exactly `nao * nao`; invalid/null buffers produce an
environment error without writing outputs. Callers remain responsible for
allocating the capacity they declare, as for other C buffer APIs.

Retention defaults to false and costs O(NAO²) additional storage. Copying a
result deep-copies its retained arrays. Changing retention or attempting a
new calculation clears the prior snapshot. Failed single points cannot expose
old basis data. Hessian and CPCM-X entry points invalidate and disable the
wavefunction snapshot; a subsequent successful single point restores it if
the result's retention flag is enabled. Hessian invalidation is exercised by
the C verifier; CPCM-X is disabled in this build and not numerically tested.

The adapter enables this capability only when the extension symbols/ABI,
reviewed patch hash, manifest capability and selected binary hash agree.
`capabilities.orbital_basis_api` is 1. GFN0 independently requires the verified
D4 fix and `capabilities.native_gfn0: true`. An engine version string cannot
establish either capability. Official C API 1/2 libraries still support basic
GFN1/2 calculations; requesting unavailable basis retention fails explicitly.

## Validation and scope

The native verifier checks GFN0/1/2 CᵀSC, shell maps, electron conservation,
GFN1/2 Mulliken charge reconstruction, result copying, disabled retention,
invalid buffers, SCF failure, Hessian invalidation and recovery. GFN0 engine
charges are EEQ, so they are not asserted equal to MO Mulliken charges.
The adapter tests cover 66 GFN1/2 frozen-conformer cases, radical/cation
calculations, compiled/ctypes agreement and degenerate-subspace invariance.
See [current validation](validation.md) for measured results.

This does not establish general libxtb thread safety. Adapter calls are
serialized, including the two Python interfaces; independent bindings are
outside that lock. Optimization, solvation and excited-state calculations
are not adapter features.

The source build is local and does not change the conda recipe's dependency
on official xTB. `pip install` installs the adapter, not a Fortran toolchain.
A portable distribution needs its own packaging and clean-environment checks.
xTB is LGPL-3.0-or-later; both license texts, exact patches, source diff and
commit are retained with the installation. Dependencies have their own licenses.
