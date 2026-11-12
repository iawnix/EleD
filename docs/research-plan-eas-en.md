# EAS case study protocol

This page is the EAS application of the broader [EleD research and case-study
plan](research-plan-eled.md). EAS is a useful atom-labelled demonstration of
EleD; it does not define the core package and it is not a claim that EleD is an
EAS predictor.

## Question

Can a researcher use one indexed EleD feature archive to train a small atom-level
classifier/ranker or graph model for electrophilic aromatic substitution (EAS)
site labels?

The case measures interface usefulness, descriptor coverage and reproducibility.
It does not require a GFN1-versus-EleD superiority result. A structure-only
feature group (A) may be reported as a sanity control; GFN1 and other published
feature recipes are optional ablations.

## Data

Use a versioned public RegioML/RegioSQM-style EAS export only after its source,
license, atom mapping and release identity have been recorded. `protocol.json`
freezes the allowed elements, event declaration, label indexing, deduplication and
leakage rules. The source adapter must reject unmapped or ambiguous labels and
must preserve every excluded row with a reason.

If the selected records are too few for an independent scaffold test set, label
the case as a pilot. Do not turn a pilot into a generalization claim.

## Electronic protocol

- Prepare at most ten ETKDGv3 conformers using the frozen seed and retry policy.
- Optimize the primary coordinates with extended-libxtb GFN2-xTB.
- Run GFN2/300 K as the primary feature calculation.
- Run GFN1/300 K on the identical coordinates as an optional method diagnostic.
- Run GFN2/1000 K as an electronic-occupation sensitivity diagnostic.
- Set `keep_overlap=True` for density, shell, frontier and π features.
- Keep engine/parameter hashes, coordinate hashes, coverage and failure rows.
- Never convert a failed or undefined descriptor into zero during feature export.

## Feature groups

| Group | Contents | Role |
| --- | --- | --- |
| A | RDKit atom, bond and local-topology features | sanity control |
| B | GFN1 charge and basic electronic features | optional method ablation |
| C | GFN2 charge and basic electronic features | optional method ablation |
| D | GFN2 plus Mulliken/Löwdin shell and density features | primary EleD case |
| E | D plus Mayer, frontier-bond and π features | exploratory extension |
| F | D plus Fukui and local conceptual-DFT features | exploratory extension |
| G | D plus conformer statistics and reliability fields | robustness extension |

Groups are reported separately. A single best concatenated score is not the case
conclusion.

## External model study

The EleD repository does not contain model training code, hyperparameter files,
datasets, trained weights or statistical results. An external research workspace
consumes the frozen feature archive and owns model selection, preprocessing,
calibration and uncertainty estimation. It must freeze molecule-level splits,
keep atom and bond indices, record its environment and report held-out metrics
without using test labels for decisions.

The public EAS protocol evaluates the archive contract: descriptor coverage,
missingness, provenance, runtime and reuse by downstream tabular or graph code.
Feature groups and any model families are declared before the external study;
EleD itself does not select a winning model or claim predictive superiority.

## Artifacts

The public archive consists of the frozen input manifest, protocol hash, feature
JSON/JSONL archive, failure report and reproduction log. External model outputs,
source exports and trained models are kept in a separate research workspace with
checksums; they are not EleD release artifacts.

`tools/prepare_eas_dataset.py` and `tools/run_eas_features.py` implement the
public preparation and feature stages. Model stages are external. No real
RegioML/RegioSQM archive is currently stored in this repository.
