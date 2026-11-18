# EleD multi-model application plan

This plan turns the EleD descriptor archive into a reproducible application
study. The external research workspace owns data adapters, model code, model
weights and numerical results. The EleD repository owns the descriptor API,
feature archive contract and descriptor-level validation.

## Research question

Can one RDKit-indexed EleD archive be consumed by conventional tabular models,
GCN, GAT, GIN and Chemprop for molecular-property regression and atom-level
reactivity prediction without recalculating electronic structure?

The study is an interface and capability evaluation. It is not a claim that a
particular descriptor family or neural architecture is state of the art.

Many existing studies already add tight-binding electronic features to machine
learning. EleD's contribution is the software handoff: stable RDKit atom/bond
identity, explicit units and missing reasons, one archive reusable by different
models, and measured calculation cost. The evidence must therefore cover both
low-data usefulness and large-data operational scaling.

## First two cases and data sources

| Case | Unit | Primary source | Target | Role |
| --- | --- | --- | --- | --- |
| Reactivity | atom | RegioSQM20 public export plus ESNUEL/ESNUEL-ML ERDA archive after rights review | EAS reactive atom set; MCA and MAA | first atom case, covering small and large data |
| Optical spectroscopy | molecule | Joung *et al.*, Scientific Data 2020; Figshare v3 `10.6084/m9.figshare.12045567.v3` | absorption and emission maxima | first molecular case |

GNN-SOM, SET-LSF and other reaction-condition datasets are deferred until these
two cases have stable archive and model contracts. They are extensions, not
requirements for the first results.

Every source receives a manifest containing source URL, commit or release,
license, file SHA256, row count, structure parse status, target coverage and
the exact filtering and duplicate policy. A mirror distributed with a model
package is never treated as the authoritative data source.

## Frozen feature groups

All models consume the same archive and all train-time transformations are fit
inside the training fold.

* **A, structure:** RDKit atom and bond features, local topology and optional
  Morgan fingerprints.
* **B/C, core electronic:** GFN1 or GFN2 charges, HOMO, LUMO, gap, dipole,
  spin, total charge and SCF diagnostics.
* **D, population:** Mulliken/Löwdin shell populations, s/p/d populations,
  density summaries and conformer statistics.
* **E, bond/frontier:** Wiberg, Mayer, HOMO/LUMO frontier overlap and pi bond
  order.
* **F, reactivity:** fixed-geometry Fukui indices, local softness,
  electrophilicity and charge-difference features.
* **G, reliability:** conformer mean, standard deviation, range, valid-count
  and failure indicators.

Null values remain null with an explicit missing reason during export. A model
adapter may impute and add a missingness indicator, but only from its fitting
partition. Zero is a valid physical value and must not be used as a missing
sentinel.

## Graph representation

The graph adapter preserves RDKit atom indices. Nodes contain structure and
EleD atom features; edges contain RDKit bond features and incident EleD bond
descriptors. Molecular quantities are carried as a global feature vector.
Solvent and reaction-condition features are separate global inputs so that a
model can be evaluated with or without environmental information.

For atom-level tasks, the graph produces one score or regression value per
node. For molecular optics, node embeddings are pooled with mean, sum or
attention pooling and concatenated with the global vector. The graph adapter
is external and must not add a dependency to the EleD runtime package.

## Model matrix

### Conventional models

Use Ridge and ElasticNet for a linear reference, ExtraTrees and Random Forest
for robust nonlinear baselines, and LightGBM/XGBoost or HistGradientBoosting
for the primary tabular comparison. For atom labels, the matrix contains one
row per atom and evaluation remains molecule-weighted. Molecular tabular
models use fixed aggregations of atom and bond descriptors: mean, sum, maximum,
minimum, standard deviation, aromatic-subset statistics and reactivity
quantiles.

### Graph models

The first implementation uses the following small models:

| Model | Initial configuration |
| --- | --- |
| GCN | 4 layers, hidden 128, residual connections, dropout 0.1 |
| GAT | 3–4 layers, hidden 128, 4 or 8 heads, dropout 0.1 |
| GIN | 4 layers, hidden 128, two-layer MLP, learnable epsilon |
| Chemprop | D-MPNN depth 3, hidden 300, FFN depth 2–3, dropout 0.1 |

The common optimizer is AdamW with learning rate `1e-3` for GCN/GAT/GIN and
`1e-4` to `5e-4` for Chemprop, weight decay `1e-5`, batch size 32–64, gradient
clipping at 5, maximum 200 epochs and early stopping after 30 validation epochs.
Hidden sizes 64/128/256, layers 2–6, dropout 0/0.1/0.2/0.3 and learning rates
`1e-4`/`3e-4`/`1e-3` are the declared search space. Use 20–50 trials per model
and task, with five predeclared random seeds.

Chemprop uses an external custom featurizer to add EleD atom and bond fields.
Its default molecular head is suitable for optics; atom-level tasks require a
node-embedding head and a molecule-preserving ranking loss.

## Task-specific objectives

For absorption and emission, train separate heads or a shared backbone with
two masked regression heads. Report MAE and RMSE in nm and eV, Spearman
correlation, calibration coverage and runtime. `e_abs_ev` is photon energy,
not an experimental HOMO–LUMO gap.

For EAS and SOM, retain multiple positive atoms. Train weighted BCE or focal
loss (`gamma=2`) and optionally a listwise ranking loss. Report molecule-level
top-1/top-3/top-5 hit rate, MRR, nDCG, PR-AUC and ROC-AUC. Global atom accuracy
is only a secondary diagnostic because molecule sizes and positive rates vary.

For MCA/MAA and other local reactivity targets, report atom-level and
molecule-level MAE, RMSE, (R^2), Spearman correlation and errors stratified
by element, aromaticity and functional group.

## Splits and leakage control

The primary split is 80/10/10 by Bemis–Murcko scaffold with a Morgan radius-2,
1024-bit cluster constraint at Tanimoto 0.7. All records, conformers and atoms
from a molecule stay in one partition. The test manifest is sealed before
feature selection or hyperparameter tuning.

Optical data also receive a solvent-transfer split: a chromophore may occur in
training, but a held-out solvent condition is evaluated separately. The model
must receive explicit solvent/environment features for this task. Duplicate
chromophore/solvent records retain their source DOI and are resolved by a
frozen policy before splitting.

## Model selection and uncertainty

Hyperparameters, imputation, scaling and feature selection are chosen inside
the training partition using deterministic grouped folds. The validation set
is reported as a development check; the final model is refit on train plus
validation and evaluates the test set once. Report every declared seed, mean,
standard deviation and molecule-bootstrap 95% interval. Calibration is
descriptive and is never selected on the test labels.

## Small- and large-sample evidence

The independent unit is always the molecule or reaction record. Counting every
atom as an independent sample would overstate evidence for atom-level tasks.

The low-data study subsamples the same frozen source at 50, 100, 250 and 500
molecules (or the largest available balanced subsets for atom labels). It uses
repeated scaffold splits and reports the learning curve, mean and standard
deviation across seeds, and the number of examples needed to reach a fixed
error or ranking target. Ridge/ElasticNet, a shallow tree model and a small GIN
or Chemprop model are sufficient for this study; large neural models are not
automatically preferred when they overfit.

The large-data study uses the complete optical archive and an ESNUEL-scale or
other openly licensed reaction archive. It measures throughput, peak memory,
parallel scaling, feature coverage, archive size, cache reuse and model fit
time. The expensive xTB stage is run once and its cost is amortized across all
tabular and graph models. A successful result shows that adding models does not
silently trigger new electronic-structure calculations.

The primary software comparisons are RDKit-only features, a direct one-off
xTB script or published descriptor export when reproducible, and the EleD
archive. Predictive improvement is one result; time to a valid model, mapping
failures, reproducibility and archive reuse are equally important results.

## EleD usability measurements

Every run also records descriptor coverage, missing reasons, conformer and SCF
failure rates, xTB engine calls, wall time, archive size, cache hits and native/
compiled numerical parity. A successful case must show that the same archive
feeds every declared model without an additional electronic-structure call.

## Execution milestones

1. Freeze the two source manifests, protocol hashes and split IDs.
2. Run 20–50 molecule/record dry runs and validate atom/bond mapping.
3. Generate the RegioSQM20, ESNUEL and optical archives and run descriptor checks.
4. Run A/B/C/D tabular baselines and archive candidate fold scores.
5. Add GCN, GAT, GIN and Chemprop adapters using the same graph schema.
6. Complete the small/large reactivity study and the optical regression study.
7. Add feature-group, GFN1/GFN2, conformer and solvent ablations.
8. Publish metrics, cost, coverage, failures, negative results and all hashes.

Model code, source data, trained weights and result tables remain in the
ignored external research workspace. EleD changes are limited to descriptor
definitions, archive contracts and validation needed to make the handoff
reproducible.
