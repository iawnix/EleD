# Orbital alignment: libxtb and tblite

> Historical 0.1–0.3 validation evidence. Version 0.4 has no tblite backend.
> Reproduce the commands below with the preserved `dist/0.3.0` source archive;
> comparison sources are archived in `tools/diagnostics/archive/tblite-0.3/`.
> Current population analysis uses the [native basis extension](xtb-patched-build.md).

The diagnostic in `tools/compare_orbitals.py` directly compares the compiled
`rdkit_eth.rdXTBTools.RunMol()` with tblite 0.8.0. It does not replace either
backend's returned coefficients. The experimental conversion helpers live in
`rdkit_eth.orbital_alignment`; tblite remains an optional dependency.

## Conventions verified from source

Rows are atomic **basis functions**, not atoms; columns are molecular orbitals.
Both codes build shells in atom order. The conversion reads libxtb's actual
parameter file (`ao=...`) and checks its shell counts and angular momenta
against tblite's zero-based orbital/shell/atom maps. This also handles repeated
s shells, such as GFN1 hydrogen's `1s2s`. It never infers a permutation by
fitting the computed coefficients.

| Shell | libxtb order | tblite order | libxtb indices for tblite order | Row phases |
| --- | --- | --- | --- | --- |
| s | s | s | `[0]` | `[+1]` |
| p | x, y, z | y, z, x | `[1,2,0]` | `[+1,+1,+1]` |
| d | x²−y², −z², xy, xz, yz | xy, yz, z², xz, x²−y² | `[2,4,1,3,0]` | `[+1,+1,−1,+1,+1]` |

Here z² denotes the normalized real solid harmonic proportional to
`2z²−x²−y²`, not a Cartesian z² function. libxtb includes the cross-term
normalization in its contracted d functions, so no additional sqrt(3) scaling
is applied to coefficients. f/g conversions are explicitly rejected.

Source references:

- [libxtb basis construction](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/xtb/basisset.f90): atom/shell loops and `set_d_function`.
- [libxtb angular transformation](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/intgrad.f90): Cartesian ordering and `dtrf2`.
- [libxtb coefficient getter](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/api/results.f90): returns the wavefunction coefficient matrix.
- [tblite v0.8.0 transformation](https://github.com/tblite/tblite/blob/v0.8.0/src/tblite/integral/trafo.f90): p permutation `[2,3,1]` in Fortran indexing and `dtrafo`.
- [tblite basis maps](https://github.com/tblite/tblite/blob/v0.8.0/src/tblite/basis/type.f90): `new_basis` shell/AO construction.

For a matrix `C` from libxtb, the diagnostic conversion is:

```python
from rdkit_eth.orbital_alignment import libxtb_to_tblite_rows

rows, phases, labels = libxtb_to_tblite_rows(
    atomic_numbers, native.parameterFile,
    tblite_result.raw_json["orbital_map"],
    tblite_result.raw_json["shell_map"],
    tblite_result.raw_json["angular_momenta"], native.numOrbitals,
)
C_in_tblite_basis = native.GetOrbitalCoefficients()[rows, :] * phases[:, None]
```

`labels` includes the atom, shell-on-atom, principal quantum number from the
libxtb parameter file, angular momentum and angular component. Matching shell
topology alone does **not** establish identical radial functions or integrals;
the metric checks below are also necessary. The C++ result itself still does
not expose a basis map or overlap matrix.

## What is checked

tblite is asked to retain its overlap matrix S. For both coefficient matrices
in tblite's AO convention, the report checks `C.T @ S @ C ≈ I`. This is a
cross-engine basis compatibility check for the libxtb matrix, not a claim
that libxtb is incorrectly normalized in its own native metric.

The report keeps raw, row-converted/sign-aligned, and subspace-aligned
coefficient differences separate. A whole MO column may change sign. For
near-degenerate orbitals, S-metric Procrustes alignment and principal-angle
projector distances compare the corresponding subspaces. Groups are admitted
only if their **full energy span in both engines** is at most `0.001 eV` and
their occupations agree within `1e-6`. These tolerances are recorded. Group
rotations are diagnostic, not new eigenvectors or evidence that slightly
different eigenvalues are exactly degenerate.

The occupied density is also compared in an orthonormal representation:
`D = (W @ C) @ diag(occupations) @ (W @ C).T`, with `W.T @ W = S`.
This removes arbitrary MO phases and rotations within equally occupied
subspaces. It is not an atom population descriptor. No core Hamiltonian is
treated as a converged SCF Hamiltonian, and no excited-state accuracy is claimed.

## Reproduce

Use the environment and compiled-extension installation in
[rdxtb-tools.md](rdxtb-tools.md), then from the project root:

```bash
source .xtb-rdkit/activate.sh
.rdxtb-build-env/bin/python -m pip install 'tblite==0.8.0'
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  .rdxtb-build-env/bin/python tools/compare_orbitals.py \
  --accuracy 0.01 --output build/orbital-comparison-patched-tight
```

Omit `--accuracy` for the default accuracy 1 run; use a different output
directory to preserve both. `--case NAME` and `--method GFN1-xTB` can select
subsets. Select another library with `--library /absolute/path/libxtb.so` and
set `XTBPATH` to its parameter directory. The conversion is source-reviewed
for the versions above; the C API version does not identify a program build.
The diagnostic rejects tblite versions other than 0.8.0 pending review.

The 33 fixed, unoptimized conformers in
[`tests/data/orbital_alignment_cases.json`](../tests/data/orbital_alignment_cases.json)
cover H/C/N/O/F/P/S/Cl/Br/I, neutral molecules, cations, anions and aromatics.
Coordinates and atom order are saved, so a rerun does not regenerate geometry.
This is an agreement fixture set, not a representative chemical accuracy dataset.

Each run writes `report.json`, `summary.json` and `report.csv`. They record
method, input geometry, numerical settings, code/fixture hashes, engine/library
and parameter hashes, differences, grouping and calculation errors (full
mapping and values are in `report.json`). An exceeded limit yields exit code 1
after writing the report; it is distinguished from an engine error. The
single-call timings include different work and warm-up costs and do not
establish a speed advantage. Library hashes are not a complete fingerprint
of every transitive BLAS/Fortran dependency.

## Observations — 2026-11-03

Stored per-case summaries and CSV files are under
[`tools/diagnostics/evidence/orbital-alignment-2026-11-03/`](../tools/diagnostics/evidence/orbital-alignment-2026-11-03/).

| Engine / accuracy | GFN1 within limits | GFN2 within limits | Calculation errors |
| --- | ---: | ---: | ---: |
| Local source build / 1 | 24/33 | 27/33 | 0 |
| Local source build / 0.01 | 33/33 | 33/33 | 0 |
| Official conda 6.7.1 / 0.01 | 31/33 | 33/33 | 0 |

For the local source build at accuracy 0.01:

| Maximum absolute difference / diagnostic | GFN1 | GFN2 |
| --- | ---: | ---: |
| Total energy, Hartree | 9.60e-14 | 1.69e-13 |
| Orbital energies, eV | 1.64e-6 | 3.23e-6 |
| Charges, e | 5.92e-7 | 7.51e-7 |
| Libxtb C in tblite overlap: orthonormality defect | 1.23e-14 | 3.67e-15 |
| Coefficients after sign alignment only | 0.275 | 0.845 |
| Coefficients after allowed subspace alignment | 2.50e-5 | 1.09e-5 |
| Subspace projector distance | 4.26e-5 | 3.08e-5 |
| Density in orthonormal representation | 2.58e-7 | 3.59e-7 |

Large sign-only differences can persist, or increase, with tighter SCF
convergence because nearly degenerate eigenvectors are not individually stable.
The much smaller subspace and density differences are the relevant evidence.

The official release has two remaining GFN1 exceptions: bromomethane and
iodomethane. Their libxtb coefficients have defects of about `1.96e-7` and
`1.89e-5` in tblite's overlap metric; iodomethane's total-energy difference is
`1.35e-6 Hartree`. The parameter-file hashes match those of the local build.
These differences are retained, not fitted away or accepted by loosening
limits. A subsequent [controlled source investigation](xtb-671-basis-differences.md)
localized Br's exception to the old on-site overlap treatment (upstream #1475)
and I's exception to the old 5s/5p radial expansion (upstream #1476). The local
source commit `a9040512...` already includes both upstream fixes. Rebuilding
that commit without our D4 patch also gives 66/66 within limits.

Wiberg matrices have another known semantic difference: libxtb's
[`get_wiberg`](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/scc_core.f90)
zeros pairs at distances >= 10 Bohr. The largest tight-run WBO differences
(`3.91e-5` GFN1, `3.77e-5` GFN2) occur for a phenol H/H pair 10.927 Bohr apart:
libxtb returns zero and tblite retains a small value. The report also records
the maximum difference restricted to distances < 10 Bohr. Tighter convergence
does not remove a definition/cutoff difference.

Regression tests additionally rotate molecular coordinates, check d-row
phases, exercise repeated s shells, distinguish equal-occupation degeneracy
from occupied/virtual mixing, and reject invalid maps and overlap metrics.
