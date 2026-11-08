# GFN0 C API state contamination: root cause and standalone reproducer

> Upstream update (2026-11-04): issue #1486 was fixed by merged PR #1490.
> The 0.4 source lock predates that repair; 0.5 uses the upstream fix. See [fix status and migration implications](upstream/xtb-gfn0-upstream-fix.md).

Investigated on 2026-11-01 using conda-forge xTB 6.7.1
(`gfortran_ha02e54d_5`, Linux x86-64) and upstream tag `v6.7.1` source.
This supplements [the original validation](validation.md). Version 0.2
rejects native GFN0 for unverified libraries; verified patched builds can
[enable native GFN0](xtb-patched-build.md). The former CLI isolation backend
has been removed.

A [second verification](gfn0-reverification.md) also compiled pristine upstream
`v6.7.1` and `main` (`a9040512`) and reproduced the issue through the public C API.
A local source patch removed the observed contamination in both builds.

## Root cause

The internal D4 initializer `newD4Model` overwrites module-level reference
charge arrays and fails to restore the GFN2 reference values when switching
back to GFN2.

1. `xtb_loadGFN0xTB` constructs an xTB calculator and reads its parameters.
   `readParam` calls `newD4Model(..., p_refq_goedecker)` for GFN0. GFN2 instead
   calls it with `p_refq_gfn2xtb`.
2. `refq(7,118)` and `refh(7,118)` are mutable Fortran module variables,
   initialized from `param_ref.fh`. They live at library/process scope, outside
   any calculator or environment object.
3. For GFN0's mode, `newD4Model` performs `refq = clsq`, `refh = clsh` and
   `secq = clss`. This destroys the original GFN2 reference tables in memory.
4. On subsequent GFN2 initialization, `secq` is reset to zero, but the first
   `select case(mode)` has no GFN2 branch restoring `refq`/`refh`. The second
   selection uses `dispm%q = refq` and `tmp_hq = refh` in its default branch.
   It consequently constructs the new GFN2 dispersion model from modified data.
5. Destroying calculators, molecules, results and environments deallocates
   those objects; it does not restore the module reference arrays. Serializing
   calls with a mutex cannot remove this sequential state dependence.

GFN2 uses charge-dependent, self-consistent dispersion. The dispersion model
contributes both to the energy and to the electronic potential (`addShift`),
so changing its reference tables can change orbital energies as well as the
total energy.

Source references, pinned to the inspected release:

- [GFN0/GFN2 parameter initialization](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/read_gfn_param.f90#L250-L306)
- [D4 table mutation and subsequent default branch](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/disp/dftd4.F90#L106-L181)
- [Module-level array declarations](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/disp/dftd4_parameters.f90#L53-L71)
- [Self-consistent dispersion potential and energy](https://github.com/grimme-lab/xtb/blob/v6.7.1/src/xtb/dispersion.f90#L109-L138)

## Independent verification

[xtb_gfn0_state.c](../tools/diagnostics/xtb_gfn0_state.c) uses only the xTB C API
for the baseline and contamination reproduction. It has fixed ethanol
coordinates and requires neither RDKit nor Python. Every invocation creates
and destroys fresh molecule, calculator, result and environment handles.

The optional `--inspect-restore` mode additionally inspects GNU Fortran internal
symbols. It snapshots and restores only `refq` and `refh` to test the mechanism.
This mode is diagnostic, compiler/version specific, and is not used by the
Python package or proposed as a production workaround.

```bash
cc -std=c11 -Wall -Wextra -Werror \
  tools/diagnostics/xtb_gfn0_state.c \
  -I/tmp/eled-check/include -L/tmp/eled-check/lib \
  -Wl,-rpath,/tmp/eled-check/lib \
  -lxtb -ldl -lm -o /tmp/xtb-state-repro

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /tmp/xtb-state-repro \
  /tmp/eled-check/share/xtb 0 --inspect-restore

# Negative control: load GFN1 instead of GFN0 in a fresh process.
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 /tmp/xtb-state-repro \
  /tmp/eled-check/share/xtb 1 --inspect-restore
```

Observed values:

| Sequence | GFN2 energy (Hartree) | Maximum orbital-energy difference from baseline (eV) |
|---|---:|---:|
| GFN2 in fresh process | -11.38860400183772 | 0 |
| Load GFN0 parameters only → delete all handles → new GFN2 | -11.38871517219716 | 0.0007635684228234822 |
| Restore only original `refq`, `refh` → new GFN2 | -11.38860400183772 | 0 |
| Fresh process: load GFN1 parameters → new GFN2 | -11.38860400183772 | 0 |

After GFN0 parameter loading, the maximum absolute changes in `refq` and
`refh` were 4.96895886239138 and 0.41407990519928, respectively. The control
GFN1 load changed neither table. Parameter loading alone is sufficient; a
GFN0 single-point calculation is not necessary to trigger the issue.

## Implications

The appropriate upstream repair is to keep the original reference tables
immutable and choose each method's tables in local/model-owned storage,
including local storage for secondary reference charges. Regression tests
should cover method-switch sequences and already-created calculators. A
[candidate patch](../tools/diagnostics/xtb_d4_local_references.patch) has now been
built and tested against both `v6.7.1` and the inspected `main`; see the
[second verification](gfn0-reverification.md) for results and limitations.

Version 0.1 isolated unverified GFN0 in a separate CLI process. Version 0.2
removes that backend and rejects GFN0 for unverified libraries. Other code calling the affected native GFN0
API directly in the same Python process could still alter libxtb state.
Recreating Python wrappers or calculator handles is insufficient. The findings
have been reported in [upstream issue #1486](https://github.com/grimme-lab/xtb/issues/1486);
releases other than 6.7.1 have not been tested here.

## Upstream awareness check — 2026-11-01

Public GitHub issues/PRs were searched in `grimme-lab/xtb` using `refq`,
`refh`, `newD4Model`, `GFN0 GFN2`, `GFN0 API`, `global`, and
`dispersion reference`; GFN0 reports in `grimme-lab/xtb-python` were also
checked before submitting this project's report. No matching report identifying
this specific D4 reference-table mutation was found. This does not establish that maintainers are privately
unaware of it or that no differently worded report exists.

There is explicit historical awareness of the broader global-state problem:

- [PR #102, Better C-API](https://github.com/grimme-lab/xtb/pull/102) says
  xTB uses many global variables and the then-current C API was not thread-safe.
- [Issue #98, maintainer discussion](https://github.com/grimme-lab/xtb/issues/98#issuecomment-574336920)
  discusses removing globals instead of merely introducing a mutex.
- [PR #157](https://github.com/grimme-lab/xtb/pull/157) and
  [PR #161](https://github.com/grimme-lab/xtb/pull/161) refactored global parameters.
- [Issue #186](https://github.com/grimme-lab/xtb/issues/186), despite its similar
  GFN0 parameter title, concerns locating the external GFN0 parameter file.
  It does not report the GFN0-to-GFN2 contamination reproduced here.

The inspected `main` commit was
[`a9040512dd7c15793472aa20f9647d36a266f7cd`](https://github.com/grimme-lab/xtb/commit/a9040512dd7c15793472aa20f9647d36a266f7cd)
(committed 2026-10-30). Its
[`newD4Model` implementation](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/disp/dftd4.F90#L106-L164)
still overwrites `refq`/`refh` for the Goedecker mode and reuses them in the
GFN2 default branch without restoration. The module-level declarations and
GFN0/GFN2 call sites also remain. Thus the inspected source still contains the
identified mechanism. The [second verification](gfn0-reverification.md) subsequently
built this exact commit and reproduced the contamination numerically, then
removed it with the local source patch.
This project subsequently submitted [issue #1486](https://github.com/grimme-lab/xtb/issues/1486)
with the complete public-C-API reproducer, numerical results, and candidate patch.
The [submitted issue body](upstream/xtb-gfn0-state-issue.md) is preserved locally.
No PR has been submitted.
