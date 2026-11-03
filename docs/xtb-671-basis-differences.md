# Why official xTB 6.7.1 met 64 of 66 comparison cases

> Historical 0.1–0.3 validation evidence. Version 0.4 has no tblite backend.
> Reproduce the commands below with the preserved `dist/0.3.0` source archive;
> comparison sources are archived in `tools/diagnostics/archive/tblite-0.3/`.
> Current population analysis uses the [native basis extension](xtb-patched-build.md).

The two GFN1 exceptions have distinct, source-level causes. They reproduce
when the official tag is built with the same compiler, flags and dependency
prefix as the newer source. Two isolated source edits remove the exceptions.
No change to the RDKit binding or AO permutation is needed.

The evidence is in
[`xtb-basis-cause-2026-11-03`](../tools/diagnostics/evidence/xtb-basis-cause-2026-11-03/).
[`investigation.json`](../tools/diagnostics/evidence/xtb-basis-cause-2026-11-03/investigation.json)
records commits, compiler settings, dependency hashes, upstream test outcomes
and SHA256 checksums of the archived scripts, patches, logs and reports.

## Bromomethane: on-site overlap was forced to one

The release's `src/xtb/hamiltonian.f90`, in `build_SDQH0`, assigns the on-site
overlap diagonal to exactly 1 and the on-site Hamiltonian to the shell's
self-energy. It computes the actual one-center integrals for multipoles but
does not use those overlap integrals to build S.

The contracted GFN1 4s/4p functions are not normalized to *exactly* one at
the precision of their tabulated expansion coefficients. For Br in the
frozen bromomethane conformer, tblite's diagonal overlap differs from one by
about `-9.54e-8` for 4s and `-2.08e-7` for each 4p component. This explains why
release xTB coefficients, orthonormal in the release's metric, fail a very
strict check in tblite's metric.

Upstream already changed this in
[#1475, Fix translational invariance for multipoles](https://github.com/grimme-lab/xtb/pull/1475),
commit `6d7cac2eb1172cb9c5b81f94e83942b48f173ac7`. It uses the contracted
on-site overlap and a consistent overlap-weighted on-site Hamiltonian, and
also makes the moment shifting consistent. The commit predates our diagnosis
but postdates release 6.7.1.

For causal isolation, the experiment backports **only the S/H0 portion** to
the release, keeping the old 5s/5p expansion and all other source unchanged.
The Br cross-metric orthonormality defect drops from `1.96e-7` to `2.07e-13`.
The full upstream patch is broader; our partial experimental patch is **not**
a production backport or a validation of gradients/multipoles.

## Iodomethane: a different 5s/5p radial expansion

Release 6.7.1 uses different six-Gaussian 5s and 5p expansions in
`src/slater.f90` than tblite 0.8.0. These are radial basis-function differences,
which cannot be removed by row permutation or orbital signs.

Upstream restored the original GFN1 shared-exponent 5s/5p expansion in
[#1476, Restore the original GFN1-xTB 5s/5p basis expansion](https://github.com/grimme-lab/xtb/pull/1476),
commit `59c8c0d9d4c1612c9201e09a66ce6039997559a2`. The commit changes the
hard-coded primitive exponents and contraction coefficients and adds heavy
element tests. It does not change the external `param_gfn1-xtb.txt` file.
This is why identical external parameter-file hashes did not imply identical
effective basis sets.

Replacing only `src/slater.f90` in the release with that upstream revision
reduces iodomethane's energy difference from `1.35e-6` to `4.48e-9 Hartree`
and its cross-metric defect from `1.89e-5` to `9.58e-10`, within the original
limits. Bromomethane's exception remains, as expected.

## Controlled experiments

All runs use the same frozen coordinates, single-thread settings, 300 K,
accuracy 0.01, compiled `rdXTBTools` and tblite 0.8.0. All source builds use
GNU Fortran 16.2 and the same compiler/linker flags and BLAS/LAPACK/mctc-lib
paths. Each variant is copied to a separate library file before the next
incremental rebuild and is loaded in a fresh Python process.

| Source variant | Br cross-metric defect | I cross-metric defect | I energy difference, Hartree |
| --- | ---: | ---: | ---: |
| Official v6.7.1 tag | 1.96e-7 | 1.89e-5 | 1.35e-6 |
| Release + on-site S/H0 only | 2.07e-13 | 1.89e-5 | 1.35e-6 |
| Release + 5s/5p expansion only | 1.96e-7 | 9.58e-10 | 4.48e-9 |
| Release + both isolated changes | 2.07e-13 | 4.86e-14 | 4.46e-9 |
| Current upstream commit, **without our D4 patch** | 1.18e-14 | 8.88e-15 | 1.78e-15 |

The baseline tag is `26b28010e805f7d1aeeef39813feb473e69cc4be`; the current
commit is `a9040512dd7c15793472aa20f9647d36a266f7cd`. The latter contains
both upstream fixes. Its two-case results reproduce the project build that
also has our D4 patch, establishing that the D4 patch is not needed to remove
these GFN1 exceptions.

The full 33-molecule, two-method comparison gives:

- Recompiled release tag: **64/66 within limits**, zero calculation errors.
- Release with both diagnostic changes: **66/66**, zero calculation errors.
- Unpatched current upstream source: **66/66**, zero calculation errors.

The remaining few `1e-9 Hartree` differences between the diagnostic release
variants and tblite do not need to be eliminated to explain the exceptions;
other source changes remain between the release and current upstream.

An additional localization probe reconstructs an effective S from the full
coefficient matrix as `inv(C).T @ inv(C)`. Its comparison with tblite places
the Br discrepancy in same-atom blocks (`2.08e-7`, interatomic difference
`1.19e-13`) and the I discrepancy mainly in interatomic overlaps (`1.94e-5`).
This reconstruction assumes a complete orthonormal MO basis; it is **not** an
overlap getter from the C API or an independent orthonormality test. The source
inspection and isolated-source interventions supply the causal evidence.

## Build/test limits and reproduction

The current unpatched source passes **45/45 upstream tests**. The raw release
tag under this newer toolchain initially passes **34/38**, with failures in
`xtb/gfn1`, `xtb/gfn2`, `xtb/hessian`, and `xtb/oniom`. The release build script
stopped before installation, and its build-directory library was used solely
for diagnostics. A follow-up isolated those four failures to a separate bug,
described below. The experiment does not certify the partial on-site patches
for production use.

The existing `.xtb-rdkit` dependency, its D4 verification manifest, and the
public backend result semantics were not changed. The project's current
pinned source already includes both upstream fixes, so an additional
production patch is unnecessary.

The archived `build_variant.py`, `run_ablation.py` and `probe_overlap.py` are
the scripts actually executed, not portable installers. To reproduce in the
same project/toolchain layout, first prepare dedicated research worktrees:

```bash
git clone --filter=blob:none --no-checkout https://github.com/grimme-lab/xtb.git \
  build/xtb-orbital-investigation/source
git -C build/xtb-orbital-investigation/source worktree add --detach \
  ../release-src v6.7.1
git -C build/xtb-orbital-investigation/source worktree add --detach \
  ../current-src a9040512dd7c15793472aa20f9647d36a266f7cd
```

Use a fresh research location/checkouts; the ablation script temporarily edits
only those dedicated release-source files and restores their original text
on successful completion. From the project root, run the archived
`build_variant.py release` and `build_variant.py current` with
`.xtb-build-env/bin/python`. The release command's upstream-test failure is
recorded above and must not be described as a successful production build.
Then run `run_ablation.py` with the same interpreter. It rebuilds the four
controlled variants, writes patches and library hashes, and invokes
`tools/compare_orbitals.py` using `.rdxtb-build-env/bin/python`.

For full comparisons, use `tools/compare_orbitals.py --accuracy 0.01`, point
`--library` at the frozen variant, and set `XTBPATH` to the research source
directory containing its parameter files. The archived build records contain
the exact executed commands and absolute paths. Full reports remain under
`build/xtb-orbital-investigation`; compact per-case reports are in the evidence
directory linked above.

## Separate gradient bug behind the four raw-tag test failures

In release `build_dSDQH0_noreset`, the expression
`g_xyz(ixyz) = g_xyz(ixyz)+stmp+dtmp+qtmp` reads `stmp` even though that local
variable has not been assigned in the routine. Its contribution depends on
compiler/runtime details. This is separate from the legitimate, assigned
`stmp` in the other gradient routine; replacing all occurrences would be
incorrect.

Upstream already removed the unused declaration, OpenMP-private entry and
addition in [#1121, Fix for GCC 14](https://github.com/grimme-lab/xtb/pull/1121),
commit `e39b5caab37b6b02cc63592a8acff14952ad5266`. Applying just these three
small edits to the raw release tag, with the same compiler settings, changes
the upstream outcome from **34/38 to 38/38 passed**. The separate Br/I orbital
comparison exceptions remain, demonstrating independence from the two basis
changes. A negative-control initialization change from upstream #1357 did
not fix these four failures; its failed test log is also retained.

To verify that this is an actual derivative problem rather than a stale test
reference, `check_gradient_fd.py` compares C API analytic gradients with
central finite differences of total energy for the upstream water geometry.
It uses fresh molecule/calculator/result handles, accuracy `0.0001`, and a
displacement of `0.0001 Bohr`, with no GFN0 calls.

| Binary | GFN1 max gradient error, Eh/Bohr | GFN2 max gradient error, Eh/Bohr |
| --- | ---: | ---: |
| Raw v6.7.1 tag, local GNU 16.2 build | 8.57e-2 | 8.79e-2 |
| Same tag/toolchain + upstream #1121 only | 1.32e-9 | 1.42e-9 |
| Installed conda 6.7.1 binary | 1.28e-9 | 1.42e-9 |

The installed conda package is `xtb-6.7.1-gfortran_ha02e54d_5`. Its tested
water gradients agree with finite differences. Thus the raw-tag build's
gradient failure must **not** be reported as an observed failure of that
conda binary, nor as an explanation for its two basis-comparison exceptions.
The precise reason this binary does not exhibit the undefined-variable
symptom was not established; one finite-difference probe is not a universal
gradient guarantee. The project source already contains the upstream fix.

The archived `check_gradient_stmp.py` reproduces the isolated patch/build/test
experiment in the research checkout and restores the source in a `finally`
block. `gradient-stmp.patch`, `gradient-stmp-ctest.log`, and the three
`*-gradient-fd.json` reports contain the exact changes and outcomes. None of
these diagnostic libraries replaces the project's installed dependency.
