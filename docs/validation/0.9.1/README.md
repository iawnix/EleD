# 0.9.1 validation and matched-geometry evidence

Recorded 2026-11-05 (Asia/Shanghai). This release fixes native libxtb discovery
from a venv inheriting the standard runtime, and adds research tooling for a
matched source/generated geometry comparison. The native engine, descriptor
formulas and physical calculation settings are unchanged.

| Check | Result |
| --- | --- |
| Development suite with extended libxtb | 359 passed, 4 skipped |
| Installed compiled wheel, isolated Python outside checkout | 359 passed, 4 skipped |
| Installed wheel with official libxtb 6.7.1 | 243 passed, 120 skipped |
| Optional benchmark and native-discovery tests in research venv | 23 passed |
| Standard installer, native/compiled GFN0/1/2 and populations | Passed |
| Final pure and compiled wheels, fresh isolated installation checks | Passed |
| Inherited venv with library/prefix overrides removed | Native and compiled population checks passed |
| Matched QM8 selected/common-complete cohort | 1,024 / 983; 41 embedding failures retained |
| Source/generated filtered split IDs | Identical |
| RDKit-only baseline predictions between geometry protocols | Identical within 1e-12 eV |

Three full-suite skips are optional sklearn tests, run separately above; the
other is an existing skip. Additional official-library skips reflect missing
capabilities. Tests check parent-prefix discovery and override priority,
deterministic geometry generation, MMFF convergence/selection, explicit failures,
invalid/frozen partitions, train/test group isolation, paired ID/target matching,
cluster bootstrap behavior and insufficient-group intervals. A test changes
only held-out targets and verifies that fitted predictions do not change.

See [the protocol and interpretation](../../qm8-geometry-study.md) for the
scientific scope and limitations. The data/model/split protocol is fixed before
execution. The source reference is the checksummed 0.9.0 pilot evidence; its
models are retrained on the same 983 successful geometries as the generated
features. No test-set tuning, molecule replacement or geometry fallback occurs.

Evidence files:

- `manifest.json`: source reference hashes, implementation/script hashes,
  selected IDs, geometry protocol, calculation settings and model/split policy.
- `geometry-preparation.jsonl.gz`: all 1,024 attempts, seeds, embedded IDs,
  optimization outcomes and selected coordinates for successful preparations.
- `generated-features.jsonl.gz`: 983 complete feature rows, input/provenance
  hashes, quality flags and timings. Source features remain in 0.9.0 evidence.
- `failures.json`: all 41 failed IDs, each `embedding_failed`; no xTB failures.
- `source-splits.json.gz` / `generated-splits.json.gz`: exact common partitions.
- `source-predictions.json.gz` / `generated-predictions.json.gz`: all 36 fits'
  held-out targets and predictions, with IDs/feature groups/splits/seeds.
- `summary.json`: individual fit MAE/RMSE, counts, software versions and costs.
- `paired-comparisons.json`: all 36 matched contrasts and conditional 95%
  intervals from 5,000 paired group-bootstrap replicates, without refitting.

On scaffold splits, mean MAE is 0.3433 eV for RDKit-only features, 0.2084 eV
with source-geometry basic xTB features and 0.2307 eV with generated-geometry
basic xTB features. Generated-geometry orbital distributions yield 0.2304 eV;
all three paired scaffold intervals for this incremental change include zero.
The basic feature benefit persists on this common pilot cohort. Additional
orbital features do not show stable incremental improvement. Geometry failures,
the reused test cohort and conditional intervals limit broader conclusions.

Generated-geometry featurization takes 113.3 s locally, including 91.9 s of
geometry preparation for all attempts; model fitting is recorded separately.
These timings exclude installation and source-reference geometry generation.

`inherited-native.json` and `inherited-compiled.json` were produced using
isolated Python in `build/benchmark-env` with `XTB_LIBRARY`, `XTBPATH`,
`XTB_PREFIX`, `CONDA_PREFIX` and `PYTHONPATH` removed. `installation.json` and
`artifact-manifest.json` identify tested local artifacts; `installed-pure.json`
and `installed-compiled.json` record final-wheel checks. `SHA256SUMS.json`
fingerprints the evidence. Source/wheel artifacts remain in `dist/0.9.1`.

The prior 0.9.0 commit `ae061cf` passed all three remote CI workflows:
41881848194, 41881848192 and 41881848177. These are not 0.9.1 CI results.
