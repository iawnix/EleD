# Research and case-study plan for EleD

## Positioning

EleD is a software contribution first. Its purpose is to give RDKit users a
stable way to calculate, retain and reuse xTB electronic-structure information
at molecular, atom, bond and fragment levels. The research cases demonstrate
that this interface can be inserted into ordinary machine-learning workflows;
they do not redefine EleD as a new predictive model or claim a new quantum
chemical population analysis.

The central question is:

> Can a researcher go from an RDKit molecule to indexed electronic features,
> quality diagnostics and a small ML/GNN experiment without writing a custom
> xTB parser or losing atom and bond identity?

## What is already part of EleD

The core package owns the calculator, extended `libxtb` runtime contract,
descriptor extraction, provenance, missing-value semantics, conformer analysis
and RDKit-native visualisation. Core tests cover units, conservation,
invariance, native/compiled parity and failure reporting. Model training and
dataset-specific code stay outside the importable `eled` package and outside
the core CI path.

## System boundary and work products

The project is organized into four layers with separate release criteria:

| Layer | Repository location | Output | Acceptance gate |
| --- | --- | --- | --- |
| Descriptor engine | `src/eled/`, `cpp/rdXTBTools/` | `XTBResult` and indexed descriptor reports | units, conservation, invariance, failure semantics and native/compiled parity |
| Feature archive | `tools/run_*_features.py`, protocol JSON | versioned JSON/JSONL with atom and bond identity | protocol hash, coordinate hash, coverage and explicit missing reasons |
| Application case | external research workspace | model-ready archive and study report | frozen data split, untouched test set and reproducible metrics |
| Research artifact | `docs/`, release archive | report, environment lock and checksums | a reader can reproduce the case without changing the core API |

The data flow is one-way:

```text
RDKit Mol + conformer
        │
        ▼
extended libxtb single point / optional optimization
        │
        ▼
XTBResult (C, S, AO maps, energies, charges, provenance)
        │
        ├── descriptor reports and RDKit annotations
        └── versioned feature archive ──► tabular model / GNN case
```

This boundary makes the software contribution independently useful. A dataset
or model can be replaced while the descriptor API, report schema and numerical
tests remain stable.

## Software validation questions

The cases measure four practical properties of the library:

| Question | Evidence |
| --- | --- |
| Can a user start from an RDKit molecule? | one documented calculator call and a report with original atom/bond indices |
| Can electronic features be reused? | one archive consumed by a tabular model and a graph adapter without recalculation |
| Is the extra information affordable? | engine-call count, wall time and storage size reported per feature group |
| Can a result be trusted or diagnosed? | units, invariance checks, coverage, missing reasons and protocol hashes |

The primary usability result is therefore a short, reproducible path from
`ROMol` to model-ready features. Predictive scores are evidence for the case
workflow; archive reuse, provenance and calculation cost are evidence for EleD.

## Case studies

### Case A: molecular tabular model

Use a small public molecular property set with fixed 3D structures or a frozen
SMILES-to-conformer preparation rule. Generate a versioned EleD feature table
and fit a compact model such as ridge regression, random forest or gradient
boosting. The case should report:

- the target, source version and license;
- the exact molecule and geometry manifest;
- RDKit-only features and EleD feature groups;
- coverage, failures, wall time and feature-generation cost;
- held-out performance with uncertainty and a short reproducible script.

The RDKit-only model is a sanity control for feature usefulness. It is not a
claim that EleD must beat a particular GFN1 method, and no GFN1 superiority
criterion is required.

### Case B: atom-level EAS model

Use a versioned, publicly available RegioML/RegioSQM-style EAS export when its
atom mapping and license can be reproduced. Keep this as an application case,
not as the definition of the EleD API. Generate node features from EleD and
edge features from the indexed bond descriptors, then run either a small atom
classifier/ranker or a message-passing GNN. The case should show that the same
feature archive can feed a table model and a graph model.

The primary report is an honest capability report: label coverage, descriptor
coverage, failure reasons, runtime, molecule-level split, and held-out ranking
metrics. A structure-only model is an optional sanity control. Comparisons to
GFN1, published RegioSQM features or other systems are optional ablations and
must not become the acceptance condition for EleD.

### Optional Case C: graph-level transfer

After Cases A and B are stable, use a second small public task to test whether
the same `XTBResult`-derived node and edge schema can be consumed by a graph
model without changing the descriptor API. This case is about schema reuse and
developer ergonomics, not about claiming state-of-the-art accuracy.

## Frozen case protocol

Each case receives its own machine-readable protocol. It records the data source
and version, molecule selection, conformer preparation, xTB method and
temperature, feature groups, split seed, model seed, metrics, and output
artifacts. The feature protocol and case protocol are hashed into every
archive. A case cannot silently change geometry, atom indexing or missing-value
handling after results are inspected.

The protocol has three checkpoints:

1. **Before calculation:** freeze source identity, selection rules, geometry
   settings, feature groups, seeds and split rules.
2. **After feature generation:** freeze the manifest, archive hash, coverage and
   failure table before opening the test set.
3. **After modeling:** publish model configuration, held-out predictions,
   uncertainty and environment information together with the archive hash.

The common electronic protocol is:

- extended `libxtb` only;
- GFN2 for the primary geometry or single-point features unless the case says
  otherwise;
- 300 K for the main analysis and 1000 K only as an electronic-occupation
  sensitivity diagnostic;
- `keep_overlap=True` for density, shell, frontier and π features;
- explicit failure rows and missing reasons; no zero imputation during feature
  generation;
- all train-time imputation, scaling and feature selection fitted inside the
  training partition.

## External model study boundary

EleD publishes the feature archive contract and the provenance needed to use it.
Model code, hyperparameter files, datasets, trained weights and statistical
results are maintained in a separate research workspace. They are not imported
by `eled`, included in the EleD test suite or shipped in the EleD repository.

An external study must consume a frozen archive, keep molecule identities and
atom/bond indices, freeze its test partition before fitting, and record its
complete environment and model configuration. The study can use tabular models
or a GNN and can report held-out metrics, uncertainty, coverage and runtime.
Those results demonstrate archive usability; they do not change the EleD API.

The public protocol records the electronic method, geometry, feature groups,
missing-value semantics and provenance fields. The external study owns model
selection, hyperparameters, calibration and statistical reporting.

## Deliverables

Each completed case produces:

1. a source and license record;
2. a selected-record manifest with atom/bond indexing;
3. a protocol JSON and SHA256;
4. an EleD feature archive with provenance and failure rows;
5. a short script that trains the model from the archive;
6. held-out predictions, metrics, uncertainty and coverage tables;
7. an environment lock and reproduction log.

The repository contains descriptor fixtures and feature-archive contract tests.
Source data, model code, trained models and statistical results belong to the
separate research workspace or an external archive with checksums.

## Milestones

| Milestone | Work | Exit condition |
| --- | --- | --- |
| M0: API freeze | descriptor catalog, units, result capabilities and missing reasons | catalog and core tests pass; no case code is required |
| M1: archive freeze | feature exporter, protocol hash, atom/bond mapping and failure rows | a 20–50 molecule dry run is reproducible |
| M2: external tabular case | hand off a frozen archive to the research workspace | external study records split, model and held-out metrics |
| M3: external atom-level case | hand off an atom-labelled archive | external study reports coverage and atom-level metrics |
| M4: external graph case | consume the same node/edge archive | GNN work proceeds without modifying EleD |
| M5: release | bilingual case report, lock file, checksums and failure analysis | a clean environment reproduces the published artifacts |

M0–M1 are software deliverables. M2–M5 are independent demonstrations that
can be published as case studies when their data and licenses are available.

Current repository status:

| Milestone | Status |
| --- | --- |
| M0: API freeze | descriptor catalog, report schema and validation suite are in place |
| M1: archive freeze | EAS protocol, exporter and dry-run contract are in place; no public source archive is checked in |
| M2: external tabular case | EleD handoff format is ready; study data and manifest remain external |
| M3: external atom-level case | feature archive interface is ready; model workspace and licensed export remain external |
| M4: external graph case | node/edge schema is specified; implementation remains outside EleD |
| M5: release | follows the first completed case |

## Execution order

1. Freeze the case protocol and source adapter before calculating features.
2. Run a 20–50 molecule dry run to inspect atom mapping, coverage and runtime.
3. Freeze the complete manifest and split; do not tune the selection after seeing
   test results.
4. Generate the feature archive and run descriptor-level validation first.
5. Hand the frozen archive to the separate research workspace for model fitting.
6. Publish the external case report, including failures and negative results,
   with the archive hash and environment lock.

EleD does not contain the EAS model adapter or training utilities. A future
public release can link to an independently versioned study repository after its
data, license and results are ready.
