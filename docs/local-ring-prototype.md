# Local ring-angle and nitrogen-planarity prototype

The first source-patch experiment identifies a useful intervention: the default
SP2 nitrogen planarity term obstructs embedding of three selected strained ring
systems. A local angle correction alone did not recover them. Removing the
targeted planarity term, with all specified chirality checks retained, did.
Combining the two changes improves some presets and worsens another; this is
an experimental ablation, not a new default or a validated general algorithm.

**Subsequent [success-retention testing](embedding-nonregression.md) rejects both
direct-replacement candidates under the declared MMFF94 quality gate.** Success
flags are retained in that larger test, but optimized energies worsen on some
previously successful inputs. The positive results below do not establish
regression-free behavior.

## Source changes

The baseline is RDKit 2026.09.1, source commit
`fece8caa860bdf6c9c82bdb9253f22a0cf87150c`, with the same SDK and instrumentation
as the [initial investigation](embedding-algorithm-research.md).

Two independent patches are supplied:

1. [`local-ring-angle.patch`](../tools/embedding_research/local-ring-angle.patch)
   changes `BoundsMatrixBuilder.cpp::set13Bounds()`. At neutral, nonaromatic,
   degree-three SP2 nitrogen in both three- and four-membered rings, it replaces
   a computed planar-completion angle **greater than 180 degrees** with the
   existing MMFF angle parameter, when available. MMFF parameterization operates
   on a copy. Missing parameters preserve the original rule. Bond parameters,
   molecular hybridization, stereochemical tags, sampling and convergence
   thresholds are unchanged. Downstream distances and triangle smoothing are
   recomputed by the original builder; this is not a matrix-entry overwrite.
2. [`ring-n-planarity-ablation.patch`](../tools/embedding_research/ring-n-planarity-ablation.patch)
   changes `Embedder.cpp::initETKDG()` by removing the improper term at nitrogen
   centers matching the same topology predicate. This second ablation does not
   depend on successful MMFF parameterization or the angle-greater-than-180 test.
   It leaves other impropers, experimental torsions, signed chiral volumes and
   final geometric/stereo checks intact. It changes the term list consumed by
   both embedding optimization and the corresponding knowledge-term check.

The second change is an explanatory ablation. **Deleting a planarity term is not
yet a replacement physical model for nitrogen pyramidalization.** The two trigger
conditions would need to be reconciled before designing a coupled production
model. In particular, unsupported MMFF chemistry can retain the old angle while
the independent planarity ablation still removes the improper.

The executable compiles the bounds-builder and embedding translation units,
linking other components from the SDK. A separately compiled, unmodified bounds
control matched the original probe exactly on all 1,024 screen molecules. The
installed RDKit and the package API are not modified.

## Protocol

The fixed 1,024-molecule QM8 cohort is screened using ETKDGv4, DG initialization,
seed 42, UFF embedding bounds and at most 100 embedding attempts. The topology
predicate selects exactly IDs **19217, 19394 and 18847**. These were already
identified during diagnosis; this is development data, not a new holdout.

Those three molecules are evaluated more deeply using:

- ETKDGv3, ETKDGv4 and srETKDGv3;
- DG, random-coordinate and internal-coordinate initialization;
- seeds 17, 29, 42, 43 and 101;
- one requested conformer and at most 100 attempts per call;
- baseline, angle-only, planarity-only and combined variants, plus the existing
  global MMFF-bounds option as a separate baseline.

This gives **135 trials per variant on three molecules**, not 135 molecules.
All successful target geometries are checked independently and then optimized
with the same MMFF94 protocol (500 iterations), followed by another stereo check.
The 1,024-molecule screen measures initial embedding only, without post-MMFF.

No source coordinates enter generation. In the angle experiment, saved source
coordinates are read only afterward for diagnostic constraint replay. The
planarity ablation does not read source coordinates at all. No thresholds were
tuned to the source angles and no source geometry was used as a template.

## Results

| Variant | Valid initial geometries, target trials | Converged MMFF94 + valid stereo, target trials | Valid initial geometries, 1,024 screen |
|---|---:|---:|---:|
| Original UFF bounds | 0/135 | 0/135 | 996/1,024 |
| Local angle correction only | 0/135 | 0/135 | 996/1,024 |
| Existing global MMFF bounds | 0/135 | 0/135 | Not screened here |
| Targeted nitrogen planarity removal only | 77/135 | 77/135 | 999/1,024 |
| Local angle correction + targeted planarity removal | 109/135 | 109/135 | 999/1,024 |

All three patched variants leave the **1,021 non-target screen molecules**
unchanged: angle-only is compared using complete probe events; the two planarity
variants are compared using all recorded bounds hashes, outcomes, failure
counters, accepted coordinates, coordinate checks and snapshot counts. No
screen regressions were observed in this single-seed protocol.

The two-factor experiment is more informative than the aggregate success count:

| Preset | Planarity removal only | Combined |
|---|---:|---:|
| ETKDGv3 | 15/45 | 25/45 |
| ETKDGv4 | **45/45** | **39/45** |
| srETKDGv3 | 17/45 | 45/45 |

The combined change gains 38 paired trials over planarity removal alone, but
loses six. All six losses are ETKDGv4 trials for ID 19394. Consequently, the
combined variant is **not uniformly better**, and its 109/135 aggregate is not
sufficient reason to replace the ETKDGv4 behavior.

## What this establishes

The original UFF-derived end-to-end bound for the C(3)–N(5)–C(6) path of ID 19217
is 2.8518–2.8923 Angstrom. The local-angle patch changes it to 2.4567–2.5410
Angstrom. The three source structures then pass the audited final distance,
chiral-volume and center-containment predicates for the local-angle variant,
but generation still fails without changing the nitrogen improper.

The original experimental-torsion setup adds an out-of-plane/improper term to
three-coordinate SP2 C, N and O. For these strained nitrogen environments, the
paired removal experiment demonstrates that this term and its associated
knowledge check are material obstacles under the tested protocol. It does not
prove the same mechanism explains every difficult ring, or that all such
nitrogens should have an unconstrained out-of-plane degree of freedom.

The generation improvements are achieved inside the RDKit embedding process.
They are not the result of turning off specified chirality, increasing the
attempt limit, changing post-embedding optimization, or supplying known 3D
coordinates. All 186 successful trials across the two planarity variants also
converge under the common MMFF94 protocol and pass the subsequent specified-stereo
check. These are repeated trials, not 186 distinct molecules.

## Decision and remaining work

Keep both patches experimental. The strongest next upstream candidate is a
carefully justified treatment of **nonplanar strained-ring nitrogen**, with the
local-angle change evaluated separately. A production model could replace the
inappropriate planar prior with a ring-aware pyramidalization prior, rather than
simply omit the term. That proposal remains to be implemented and validated.

Before proposing a default change upstream, reproduce the cases on a newer
complete RDKit build, test additional independent ring scaffolds and negative
controls, and investigate the six combined-variant regressions. Atom/bond-order
robustness, strain/energy distributions, conformer diversity and matched wall-time
budgets also remain open. The present checks do not establish universal chemical
validity or absence of regressions outside the tested cohort. Existing RDKit
MMFF typing and bounds are reused; no algorithmic novelty is claimed.

## Reproduction and verification

See [commands and tool descriptions](../tools/embedding_research/ring-prototype.md)
and [frozen evidence](research/local-ring-prototype-v1/README.md).
The focused test run passed **46 tests, with one existing integration skip**.
The native tests exercise local scope, unchanged ordinary/single-ring molecules,
unchanged global MMFF behavior, missing-MMFF fallback, downstream bounds
propagation and retained chirality checks in the successful combined case.

The public package remains version 0.10.0 with its previous runtime behavior.
