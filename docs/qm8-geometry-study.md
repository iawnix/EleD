# Matched QM8 geometry study (0.9.1)

The [original pilot](qm8-pilot.md) used dataset-provided ground-state geometries.
This follow-up tests a practical input pipeline: canonical SMILES → ETKDGv3 →
MMFF94 → ground-state GFN2 descriptors → the same fixed regression model.
The target remains the original QM8 E1-CC2 vertical excitation energy at its
source geometry. These are approximate input features for predicting that target,
not new reference excitation energies at the generated geometries.

## Frozen comparison

We reuse the 1,024 selected IDs, original feature definitions, ExtraTrees settings
and six train/test partitions from 0.9.0. All source evidence is checksum-verified.
No model or geometry parameters are tuned on these test results.

For each selected isomeric SMILES, construct a fresh molecule with explicit
hydrogens and request eight ETKDGv3 conformers. The seed is a deterministic hash
of the protocol name and source ID; RMS pruning is 0.2 Å and one thread is used.
MMFF94 optimizes each embedded conformer for up to 500 iterations. Select the
lowest-energy **converged** MMFF conformer, with its ID breaking ties. Record
all statuses/energies and the selected coordinates. The workflow has no fallback,
no xTB geometry optimization and no energy-based selection using the target.
It uses source-derived connectivity/stereochemistry from SMILES, but no source
coordinates for generation. The reconstructed molecule must reproduce the
original 2D feature vector within numerical tolerance.

**41/1,024 molecules fail embedding** under this fixed protocol; none are
silently replaced or retried. All remaining 983 complete MMFF and xTB feature
extraction. These failures may depend on the representation and difficult
geometries in this cohort; they do not prove the molecules cannot be embedded
with other settings. Exact failed IDs and attempts are retained.

Both source and generated-geometry models are fitted again on the **same 983
successful molecules**. Original split memberships are retained; removed rows
are simply deleted from both sides. No new random split is drawn. Consequently,
the source scores below differ from the earlier 1,024-molecule pilot and must
not be compared directly to that table as though the cohorts were unchanged.

GFN2 settings stay at 300 K electronic temperature, accuracy 1 and independent
initial guesses. The current libxtb binary and parameter-file hashes must match
the reference feature evidence, preventing an engine change from being confused
with a geometry change. Absolute installation paths need not match. This study
does not average conformers or implement Boltzmann weighting.

## Results on the common cohort

Test MAE in eV, mean ± sample standard deviation over three fixed seeds:

| Geometry / features | Random split | Scaffold split |
| --- | ---: | ---: |
| RDKit 2D + Morgan (same in both protocols) | 0.3519 ± 0.0075 | 0.3433 ± 0.0233 |
| Source geometry + basic xTB | 0.2138 ± 0.0083 | 0.2084 ± 0.0226 |
| Source geometry + basic xTB + orbital distributions | 0.2048 ± 0.0098 | 0.2070 ± 0.0263 |
| Generated geometry + basic xTB | 0.2386 ± 0.0123 | 0.2307 ± 0.0221 |
| Generated geometry + basic xTB + orbital distributions | 0.2317 ± 0.0137 | 0.2304 ± 0.0255 |

Basic xTB features still help on the generated geometries in this pilot. Their
average MAE is about 0.022–0.025 eV worse than with source geometries. Source
geometry quality therefore matters, but does not explain the entire observed
benefit over 2D features in this common cohort. This conclusion is conditional
on successful geometry generation; the roughly 4% failure coverage is part of
the practical result.

Additional orbital distribution features show only small mean improvements.
For **generated geometry, scaffold splits**, the candidate-minus-basic MAE
differences and exploratory paired 95% intervals are:

| Seed | MAE difference (eV) | 95% interval (eV) |
| --- | ---: | ---: |
| 17 | +0.00220 | [−0.00609, +0.01093] |
| 29 | +0.00107 | [−0.00570, +0.01156] |
| 43 | −0.00415 | [−0.01567, +0.00311] |

All three intervals cross zero. We cannot claim a stable additional benefit
for these orbital features under this protocol. Every contrast, including the
other splits and unfavorable outcomes, is archived rather than selecting only
positive results.

## What the intervals mean

Predictions are paired by molecule ID and checked for identical targets. For
each fitted model pair and each split, compute the difference in absolute error
on each test molecule. Resample whole test groups with replacement 5,000 times:
chemical identity for random splits, Bemis–Murcko scaffold for scaffold splits.
All acyclic molecules remain in the single empty-scaffold group. Each sampled
group contributes all its molecules; MAE remains molecule-weighted. Report the
2.5th/97.5th percentiles of candidate-minus-reference MAE; negative favors the
candidate. A single group has no reported interval, and fewer than 20 groups
are flagged. Here scaffold test sets contain 67–69 groups.

These intervals are **conditional on the already trained models and current
test cohort**. Models are not refitted inside the bootstrap. The procedure does
not include training-set/model-selection uncertainty, geometry failures,
reference-label error or multiple-comparison adjustment. The seeds reuse
molecules, so intervals are never pooled across seeds as independent samples.
This is an exploratory robustness check, not an independent external validation
or proof of general statistical significance.

## Cost and reproduction

In the recorded local run, generation/MMFF across all attempts took 91.9 s,
xTB on successful rows 16.9 s and descriptor analysis 3.3 s. Overall
featurization took 113.3 s. Model fitting/prediction is recorded separately for
all 36 fits. Native/BLAS computation uses one thread and each forest two workers.
These numbers exclude runtime installation and creation of source geometries.

Install 0.9.1 first: it fixes native-library discovery from a venv inheriting
the standard runtime. Keep model dependencies outside that locked runtime:

```bash
.rdkit-eth/runtime/bin/python -m venv --system-site-packages build/benchmark-env
build/benchmark-env/bin/python -m pip install -r tools/benchmark-requirements.txt
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  build/benchmark-env/bin/python -I tools/benchmark_qm8_geometry.py \
  --output build/qm8-geometry-new-run
```

The default reference is `docs/validation/0.9.0`; `--reference` accepts an
equivalent evidence directory with the same named files and checksum manifest.
A different native binary/parameter hash requires regenerating the source
reference features with that engine before comparing. The tool intentionally
refuses to mix engines or overwrite an existing output directory.

[Archived evidence](validation/0.9.1/README.md) includes preparation outcomes and
selected coordinates, generated features, failure IDs, both sets of predictions,
exact filtered splits, all paired intervals, settings/hashes and timings. It
also records the isolated-venv fix and package installation checks. The model
dependencies are still research tools, not package runtime requirements.
