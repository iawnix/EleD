# Initial embedding algorithm investigation

The first source-level study supports investigating **local constraint modeling
for strained ring systems**, alongside sampling. It does not establish a new
algorithm or an improvement over RDKit. In particular, existing RDKit options
already address part of the observed difficulty and must be included as baselines.

A subsequent [local-ring source-patch experiment](local-ring-prototype.md) tests
angle correction and targeted nitrogen planarity removal separately. It records
both a negative angle-only result and improvements/regressions in the factorial
ablation; the patches remain research-only.

## What was implemented

A standalone diagnostic executable instruments RDKit 2026.09.1's embedding
implementation. It captures the actual bounds matrix, chiral/tetrahedral sets,
torsion atom quadruplets and rejected intermediate coordinates. A Python audit
replays distance, signed-volume and center-containment predicates with atom-level
details, then applies the same predicates to mapped QM8 source geometries.

The executable never reads source coordinates. Reference geometries are loaded
only after generation finishes. Logging changes no coordinates, constraints or
random numbers; the installed RDKit is used as an independent parity check.
This research tooling adds no production dependency or public API.

## Frozen protocol and results

The cohort is the existing 81-case preparation audit: 41 original embedding
failures (22 subsequently accepted by `prepare_conformers`, 19 still rejected)
and 40 controls. The latter 19 are **workflow failures**, not all failures of
initial embedding: IDs 17445, 12123 and 17312 can embed, but previous MMFF94
relaxation produced coordinates rejected by the specified-stereo check.

Each independent run requests one conformer, seed 42, one thread and at most 100
embedding attempts. Chirality remains enforced. Presets are ETKDGv3, ETKDGv4 and
srETKDGv3; initial coordinates use DG or the existing internal-coordinate mode.
Both existing `embedForceField=UFF` and `embedForceField=MMFF` settings are tested.
**MMFF here parameterizes embedding bounds; it is not post-embedding MMFF94
optimization.** There is no force-field relaxation or xTB calculation afterward.

Counts below mean returned coordinates pass the package's finite-coordinate,
collision and specified-stereo checks, before any relaxation:

| Preset / bounds | Difficult group, DG / internal (of 19) | Previously recovered, DG / internal (of 22) | Controls, DG / internal (of 40) |
|---|---:|---:|---:|
| ETKDGv3 / UFF | 0 / 0 | 0 / 1 | 39 / 36 |
| ETKDGv3 / MMFF | 4 / 1 | 9 / 6 | 37 / 36 |
| ETKDGv4 / UFF | 3 / 3 | 18 / 19 | 40 / 40 |
| ETKDGv4 / MMFF | 9 / 9 | 20 / 20 | 39 / 39 |
| srETKDGv3 / UFF | 3 / 3 | 9 / 9 | 39 / 37 |
| srETKDGv3 / MMFF | 7 / 4 | 12 / 11 | 39 / 38 |

Across all tested settings, 10 of the difficult 19 have a passing initial
geometry. This is a union across unequal total work, **not** a new method's
success rate, and does not establish end-to-end recovery after relaxation.
Equal attempt limits also do not imply equal optimizer work or wall time.

All **972** probe/library comparisons passed. **2,628** captured snapshots passed
the geometric replay checks. All 81 mapped source geometries passed the package
stereo check; no specified center had a wrong affine chiral-volume sign.
One molecule, ID 3250, triggers the same `bad direction in linearSearch`
optimizer exception in the probe and installed RDKit under all four srETKDGv3
settings. These are counted as exceptions, not ordinary embedding failures.

## Concrete constraint findings

`FINAL_CHIRAL_BOUNDS` is misleading if read as a wrong-handedness diagnosis. The
implementation checks **pair distances within the union of specified chiral
centers and their neighbors**, allowing a deviation of 10% of the upper distance
bound. Signed chiral-volume checks are separate.

With the default UFF-derived bounds, source geometries of **7 of the difficult
19** fail at least one audited geometric predicate:

| Predicate | Source IDs |
|---|---|
| Final chiral-neighborhood distance bounds | 6596, 8594, 8596, 18847, 19217, 19394 |
| Specified center containment/margin | 8594, 19214 |
| Unspecified tetrahedral-center volume/containment | 8596 |

These categories overlap. ID 18185 in the previously recovered group also fails
the source distance check. None of the 40 controls does under UFF bounds.
Source geometries are useful reference structures, not guarantees of a global
minimum or universal chemical correctness. Rejection of a reference geometry
does **not** prove that the constraint set has no feasible solution, nor that
this one predicate caused all failed generation attempts.

Three cases share a useful mechanism hypothesis: a nitrogen assigned SP2 lies in
three- and four-membered rings. `BoundsMatrixBuilder.cpp::set13Bounds()` assigns
ring angles first, then distributes the remaining planar angle around SP2 atoms.
The resulting distance bounds differ substantially from the source structures:

| ID | Canonical explicit-H atom path | Source angle | Source end-to-end distance | UFF distance interval |
|---|---|---:|---:|---:|
| 19217 | C(3)–N(5)–C(6) | 105.05° | 2.3398 Å | 2.8518–2.8923 Å |
| 19394 | C(1)–N(2)–C(3) | 106.02° | 2.3639 Å | 2.8518–2.8923 Å |
| 18847 | C(1)–N(7)–C(5) | 108.62° | 2.4260 Å | 2.8518–2.8923 Å |

For ID 19217, existing MMFF bounds change this interval to 2.3833–2.4662 Å and
remove the source's final distance-check failure. However, ETKDGv4 still fails
to generate it in this protocol: the distance discrepancy alone does not explain
everything. MMFF also introduces reference distance rejections elsewhere and
reduces ETKDGv4 control coverage from 40/40 to 39/40. Simply switching all molecules
to MMFF bounds is not a complete solution or a new research contribution.

Center-containment failures also need careful interpretation. At the source
centers of IDs 8594 and 19214, the center is on the correct side of every ligand
face, but one unnormalized face-volume margin is below the fixed 0.1 threshold.
For example, ID 8594's margin magnitude is about 0.0154. Its four-ligand affine
volume is 8.4375, above the chiral minimum and of the requested sign. Calling this
simply a chirality inversion would be incorrect. ID 8596 additionally has an
unspecified center outside one ligand face and a normalized volume of about
0.0341, below the existing fused-small-ring threshold of 0.125.

## Prior-art boundary

This is a targeted review of implementation, publication metadata and upstream
issues, not an exhaustive novelty review. Full papers have not all been reviewed.

- [Riniker and Landrum, 2015](https://doi.org/10.1021/acs.jcim.5b00654): ETKDG's
  use of existing conformational knowledge is established prior art.
- [Wang, Witek, Landrum and Riniker, 2020](https://doi.org/10.1021/acs.jcim.0c00025):
  small-ring and macrocycle adaptations are established prior art.
- [Liberti, Lavor and Maculan, 2008](https://doi.org/10.1111/j.1475-3995.2007.00622.x):
  branch-and-prune for molecular distance geometry already exists. Local
  feasibility and combinatorial reconstruction are not novel by themselves.
- [RDKit 2026.09.1 release notes](https://github.com/rdkit/rdkit/blob/fece8caa860bdf6c9c82bdb9253f22a0cf87150c/ReleaseNotes.md):
  ETKDGv4 and internal-coordinate initialization are available as beta features.
  The implementation also already has signed-volume constraints, fused-ring
  tolerance flags, MMFF bounds and Z-matrix chirality correction.
- [PR #9740](https://github.com/rdkit/rdkit/pull/9740), open when checked:
  user-provided chiral-volume constraints are already being developed upstream.
- [Issue #9760](https://github.com/rdkit/rdkit/issues/9760), open when checked:
  distorted-center 3D stereo assignment can depend on neighbor ordering. A future
  method needs permutation-aware geometric validation beyond one assigned tag.
- [Issue #8715](https://github.com/rdkit/rdkit/issues/8715), closed 2026-10-16:
  a torsional-gradient issue was fixed after the pinned source tag. The SDK's
  source/patch provenance must be checked and newer fixed baselines tested before
  attributing optimizer problems to a novel algorithmic limitation. This study
  does not establish that issue as the cause of ID 3250's exception.

## Proposed algorithm experiment

The candidate is **ring-system-aware local constraint generation and repair**.
It must change the constraints/optimization inside embedding to qualify as this
research direction; another preset retry policy would not test the hypothesis.

1. Identify coupled fused/bridged ring subgraphs from the molecular graph. Build
   local bond-length and angle constraints jointly, retaining their dependencies
   instead of assigning every angle independently from a planar/regular-ring
   template. A local Gram-matrix or equivalent coordinate formulation can test
   geometric compatibility; this mathematical ingredient itself is established.
2. Preserve specified handedness and connectivity. Treat empirical angle ranges
   and geometric-quality margins separately from stereo signs. Diagnose which
   terms are violated and distinguish a local sampling failure from evidence of
   an inappropriate prior. A failed optimizer is not an infeasibility proof.
3. Reconstruct the implicated ring subgraph, retaining compatible external
   geometry; enlarge the subgraph if ring closure or boundary constraints require
   it. Introduce torsional preferences after obtaining a feasible local geometry,
   with subsequent constrained refinement and independent stereo/geometry checks.
4. First test two separate ablations: local constraint modeling with unchanged
   global sampling, then local re-embedding with unchanged constraints. Only then
   combine them, so benefits and regressions are attributable.

The current evidence most directly motivates the first ablation. It does not yet
demonstrate that local reconstruction will improve these molecules. A useful
first target is the small-ring nitrogen class, but neither its three source
angles nor source coordinates may become molecule-specific generation inputs.

Before evaluating novelty/performance, use an independent ring-scaffold holdout,
multiple seeds, matched time/optimization budgets and stronger current RDKit
baselines including MMFF bounds, random coordinates and internal coordinates.
Report geometry validity, stereo invariance under atom/bond reorderings, strain,
diversity and downstream relaxation stability separately. The selected QM8 hard
cases are diagnostic development data, not a generalization benchmark.

## Evidence and reproduction

[Research tools](../tools/embedding_research/README.md) include build/run commands.
[Frozen evidence](research/embedding-constraints-v1/README.md) contains complete
constraints, capped failure snapshots, source mappings, summaries, source and
binary provenance. The package remains version 0.10.0; its defaults and historical
QM8 prediction results are unchanged.
