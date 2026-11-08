# Development priorities

As of 2026-11-05, embedding-algorithm research is paused at the user's request.
Keep its source patches, negative results and frozen evidence for reproducibility.
Neither experimental patch is promoted to the package, and no new embedding
model or failure-only rescue integration is scheduled. The existing
`PrepareConformers` API remains input preparation using RDKit algorithms.

The active project is RDKit molecular input → extended libxtb → indexed,
unit-defined electronic descriptors through the Python and C++ interfaces.
The sole production calculation backend remains libxtb. Its Hamiltonians are
xTB methods; the project's contribution is the interface, population access,
descriptor analysis, provenance and reproducible workflows.

## First priority: complete descriptor validation

The historical generated-geometry QM8 study includes nine geometries subsequently
rejected by coordinate-derived stereo checks. It is not yet a fully validated
prediction benchmark. Complete this evidence before adding descriptor families:

1. Produce a separately named correction using the audited accepted historical
   geometries. Preserve original evidence, feature definitions, model settings
   and split memberships. Refit source and generated-geometry comparisons on the
   identical retained cohort, including the 2D baseline. Report excluded IDs and
   coverage against all 1,024 selected molecules. Filtering old predictions alone
   is insufficient because rejected rows may have been in training sets.
2. Evaluate the current checked `PrepareConformers` workflow separately with a
   frozen protocol. This changes geometry selection and retry budgets, so it must
   not be presented as merely correcting the old study. Record failures, costs,
   stereo checks, engine hashes and matched feature/model comparisons.
3. Report paired uncertainty and unfavorable results. Existing orbital-distribution
   features have not shown a stable additional benefit under scaffold splits;
   ground-state descriptors are not direct excitation-energy calculations.

Both QM8 evaluations are now complete; see [the 0.11 results](mainline-validation.md).
Model dependencies stay in the research environment rather than the installed runtime.

## Following priorities

- Extend scientific validation to independent chemical coverage and an additional
  relevant target family with a protocol fixed before inspecting test outcomes.
  New descriptor definitions should answer a specific analysis need and include
  units, invariance, missing-value semantics and interface parity checks.
- Improve installation and release portability from the current verified Linux
  x86-64 / Python 3.11 source-install profile. A plain pip install does not yet
  provide the full native runtime; broader platform support requires its own
  build and installed-package evidence.

Current functionality is described in [the README](../README.md), with detailed
[descriptor definitions](descriptors.md), [reactivity protocols](descriptor-reports.md),
[the EAS research plan](research-plan-eas-en.md), and [validation history](validation.md).
Historical embedding findings and the
reason for rejecting direct replacement remain in the
[nonregression report](embedding-nonregression.md).

## 0.11 delivery status

- Completed the two QM8 evaluations and a 256-molecule QM9 dipole comparison.
  Independent non-GDB coverage and IP/EA/reactivity calibration remain future
  validation, without implying that current descriptors are experimentally accurate.
- Added disposable libxtb process workers, per-task timeouts and parent-owned
  cache/checkpoint commits for single-point, descriptor and reactivity batches.
- Added successful-SCC diagnostic exports in the extended engine and both
  bindings. Failed-SCF traces and spin-resolved densities remain unavailable.
- Added offline Linux runtime bundle build/install tools with relocation checks.
  Public package-index delivery, Windows/macOS and a broad compatibility matrix
  remain outside the validated profile. The package still requires separate
  installation from upstream RDKit.
