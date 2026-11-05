# Preserve successful embedding behavior

**Status (2026-11-05): this algorithm route is paused.** Future integration
ideas below are historical proposals, not scheduled work. Active priorities
are recorded [separately](development-priorities.md).

**Neither current ring prototype passes the conservative quality-retention gate
for direct replacement of RDKit behavior.** Expanded tests preserve baseline
success counts, but some previously successful inputs converge to higher-energy
MMFF94 minima. Recovery of other failures does not offset these regressions.

## Acceptance contract

Compare identical input, seed, preset, initialization and attempt budget, with
the same post-embedding optimization:

1. Every baseline initial coordinate/stereo success must remain successful.
2. Every baseline converged-MMFF94/stereo success must remain successful.
3. For paired valid optimized structures, flag an MMFF94 energy increase above
   **0.001 kcal/mol**. This is a predeclared conservative numerical allowance,
   not a universal chemical accuracy threshold.
4. Non-target controls must retain their exact probe events. Also require
   coverage of **baseline successes inside the patch's target scope**. If this
   stratum is empty, zero observed losses is inconclusive, not a pass.
5. Runtime, diversity and downstream descriptor accuracy require separate
   evaluation. Finite tests cannot certify all chemistry or all quality metrics.

The [gate](../tools/embedding_research/nonregression_gate.py) returns
`reject_observed_regression` for any observed success/energy regression, even
with many recovered failures. It reports missing critical coverage separately
and never promotes finite test results to `production_ready`.

## Expanded cohort and protocol

The checksum-verified QM8 archive has 21,786 records and 20,309 eligible unique
identities under the existing pilot's filtering. The topology predicate matches
43; three belong to the original development cohort. The remaining **40 matching
identities** are frozen before generation. Twelve have scaffolds absent from the
original cohort, representing **six distinct new scaffolds**. This is additional
coverage within the same dataset, not an external validation dataset.

Forty nonmatching controls are selected in the existing stable hash order,
excluding development identities. The main experiment comprises 40 matching
molecules × 3 presets × 3 initialization modes × 5 seeds = 1,800 paired settings
with UFF bounds, plus 40 controls × ETKDGv4/DG × 5 seeds = 200 settings.

A supplementary experiment tests the same 40 matching identities with ETKDGv4/DG
and five seeds in three strata: specified-stereo SMILES with MMFF bounds,
unspecified-stereo SMILES with UFF bounds, and unspecified-stereo SMILES with MMFF
bounds. These 600 settings use **the same 40 identities**, not 120 new molecules.
Removing input stereo annotations changes the specification of the task;
`enforceChirality` remains enabled. No source geometry initializes generation.
Each embedding call retains 100 attempts; post-optimization always uses MMFF94
with at most 500 iterations. The binaries are unchanged from the prior study.

## Results

The 2,600 paired settings require 7,800 executable runs. They contain **443
baseline initial successes** and **433 baseline post-optimization successes**.
Both candidates preserve all of these success flags in this experiment.

| Stratum | Baseline initial successes | Lost initial/post successes, either candidate | Higher optimized energy: planarity-only / combined |
|---|---:|---:|---:|
| Specified matching inputs, UFF, full grid | 8/1,800 | 0 / 0 | 1 / 1 |
| Nonmatching controls | 190/200 | 0 / 0 | 0 / 0 |
| Specified matching inputs, MMFF, ETKDGv4/DG | 0/200 | No affected baseline-success coverage | 0 / 0 |
| Unspecified stereo, UFF | 130/200 | 0 / 0 | 5 / 17 |
| Unspecified stereo, MMFF | 115/200 | 0 / 0 | 10 / 10 |

Planarity-only has **16** flagged quality comparisons and combined has **28**.
These are repeated variant/settings comparisons, not distinct molecule counts.
The specified-MMFF stratum is inconclusive for retention because the original
generator had no passing result there. Candidate successes cannot supply that
missing coverage.

All 200 nonmatching control trials have identical probe events in both
candidates: 190 initial successes and 180 post-optimization successes. The ten
initial failures and ten post-optimization failures are also retained and are
not counted as new regressions.

## Concrete counterexamples

| Source ID and input | Settings | Variant | Baseline → candidate optimized MMFF94 energy |
|---|---|---|---:|
| 11734, `CC[C@@]12CC(=O)N1C2` | srETKDGv3, random, seed 17, UFF bounds | Combined | 16.6617 → 18.0638 kcal/mol |
| 11734, same input | srETKDGv3, internal, seed 43, UFF bounds | Planarity-only | 16.6617 → 18.0638 kcal/mol |
| 11686, unspecified `N#CC12CC(=O)N1C2` | ETKDGv4, DG, seed 17, MMFF bounds | Both | 11.4502 → 67.6257 kcal/mol |

All listed structures pass specified-stereo checks and optimizations report
convergence. The largest increase, **56.1756 kcal/mol**, was rerun independently
with identical results. Its maximum bonded-distance/covalent-radii-sum ratio
changes from about 1.022 to 1.190; this descriptive ratio is not a new validated
chemical-validity threshold.

Higher energy can reflect a different local minimum or, for unspecified inputs,
a different allowed stereo outcome. It does not alone prove an invalid
structure. It does demonstrate that convergence and stereo flags are
insufficient to claim preserved quality under the declared metric. The
1.4021 kcal/mol examples retain explicitly specified stereo and use the same
molecule/force field.

Wall times include subprocess startup and trace output, with rotated execution
order. The supplementary run overlapped part of the main experiment. These
timings cannot certify production speed retention; no performance gate passed.

## Integration decision

Keep both direct-replacement patches experimental. Do not change the default
generator or propose them as a proven regression-free upstream replacement.
Earlier positive development results remain valid within their scope; the
expanded retention test provides counterexamples to the stronger claim.

For future compatibility-oriented integration, preserve a successful baseline
result unchanged and enter the new constraint solver only after a defined
baseline failure. This preserves old successful output by design but adds work
on rescued failures. Such an integration policy is neither an equal-budget
speed improvement nor algorithmic novelty by itself. It is **proposed**, not an
implemented public API here. Constraint-model research can continue with these
regression cases as mandatory controls. Package 0.10.0 behavior is unchanged.

## Reproduction

Build the research executables using [the prototype instructions](../tools/embedding_research/ring-prototype.md), then:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build/install-validation/sdk/bin/python tools/embedding_research/check_ring_nonregression.py \
  --data build/qm8 --baseline build/ring-angle-control/embedding_probe \
  --planarity-only build/ring-planarity-only/embedding_probe \
  --combined build/ring-angle-and-planarity/embedding_probe \
  --output build/ring-nonregression-v1

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build/install-validation/sdk/bin/python tools/embedding_research/check_ring_input_variants.py \
  --cohort build/ring-nonregression-v1/cohort.json \
  --baseline build/ring-angle-control/embedding_probe \
  --planarity-only build/ring-planarity-only/embedding_probe \
  --combined build/ring-angle-and-planarity/embedding_probe \
  --output build/ring-input-retention-v1
```

Unexpected execution errors abort; observed regressions are recorded while the
remaining frozen cohort completes. Read `summary.json` verdicts: successful
process completion does not mean a candidate passed. Outputs are incremental,
without resumption. [Frozen evidence](research/ring-retention-v1/README.md)
includes outcomes, optimized coordinates/energies, regression traces and tests.
