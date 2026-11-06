# Research plan: RDKit-native xTB density descriptors for EAS selectivity

## Scientific question

On an atom-labelled electrophilic aromatic substitution (EAS) task with strict
scaffold and similarity-cluster separation, do shell-resolved density features
from one xTB calculation add stable site-ranking information beyond RDKit
structure features and a GFN1 partial-charge baseline?

The project does not claim to invent Mulliken, Löwdin or Fukui definitions. Its
contribution is a reproducible RDKit/libxtb interface, an atom-level benchmark,
and a reliability analysis that separates chemical signal from conformer and
numerical sensitivity.

## Pre-registered feature groups

- **A:** RDKit atom, bond and local-topology features;
- **B:** A plus GFN1 charge, shell and Wiberg features;
- **C:** corresponding GFN2 features;
- **D:** C plus Mulliken/Löwdin s/p/d populations and p/s ratios;
- **E:** D plus Mayer, frontier-bond and π features;
- **F:** D plus Fukui, local softness and local electrophilicity;
- **G:** D/F plus conformer mean, standard deviation, minimum and maximum.

A--D define the main conclusion. E--G are exploratory extensions and are reported
separately. Missing values remain explicit masks with reasons; they are never
replaced by zero.

## Frozen calculation protocol

Use the public RegioML/RegioSQM EAS atom labels and the selection rules in
`protocol.json`. RDKit generates and checks the conformers. GFN2 supplies the
predeclared geometry, while GFN1 and GFN2 single points use identical coordinates.
The extended `libxtb` build retains `C`, `S` and AO/shell maps for density analysis.

The primary electronic temperature is 300 K. A second 1000 K calculation is a
deliberately amplified occupation-sensitivity diagnostic; it is not a nuclear
temperature, thermal conformer ensemble or excited-state calculation.

## Models and evaluation

Use regularized logistic regression and gradient-boosted trees. Keep molecules
sharing a Bemis--Murcko scaffold or Morgan Tanimoto cluster in the same partition.
Fit scaling, feature selection and imputation inside training folds only. The
untouched test set is evaluated with molecule-weighted macro ROC-AUC, mean
reciprocal rank, top-1 and top-2 accuracy, with PR-AUC and calibration as
secondary measures.

For every candidate minus GFN1 baseline comparison, use 2,000 molecule-level
bootstrap replicates and report the difference with a 95% confidence interval.
A positive result requires at least one primary metric whose interval is entirely
above zero without significant deterioration in the other primary metrics. If
this condition is not met, report no stable incremental benefit.

## Required technical evidence

Before model fitting, archive AO permutation, orbital-phase, degenerate-rotation
and rigid-motion invariance tests; `trace(P S)` and shell-to-atom conservation;
ctypes/compiled parity for GFN0/1/2; coverage and failure reasons; engine and
parameter hashes; and conformer-level atom/bond index checks.

The final package contains the frozen protocol and hash, selected-record manifest,
feature files, failure manifest, training scripts, locked environment, test logs,
models and bootstrap results. A negative result is retained as a valid outcome.
