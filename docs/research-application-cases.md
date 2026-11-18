# Application case research

This page records the candidate application cases for EleD and the decisions
behind them. Dataset files, model code, trained weights and numerical results
remain in the external research workspace.

## What an EleD case must show

An appropriate case has four properties:

1. labels are attached to atoms or bonds so RDKit identity is part of the
   problem;
2. the source, version and data rights can be recorded;
3. a molecule-level split can prevent close analogues from crossing the test
   boundary; and
4. the same EleD archive can feed a small tabular model and a graph model while
   reusing the recorded electronic-structure calculation.

The result demonstrates a usable electronic-descriptor interface. EleD supplies
the descriptor archive and provenance; downstream research supplies the task
models and compares feature groups under a declared protocol.

## Why the optical database is a strong EleD case

The optical database is a better molecular-level companion to the atom-level
EAS case. The paper by Joung *et al.* describes 20,236 literature data points,
7,016 chromophores and measurements in solvents or solid state. Freeze the
current version-3 Figshare file for a new study:
`10.6084/m9.figshare.12045567.v3`, licensed CC BY 4.0. Record the exact
download URL, file SHA256 and row count for this release alongside the rounded
paper statistics and the DeepMPP copy.

The DeepMPP repository provides model-training code and a convenient
`_Data/optical.csv` example. The external study cites and freezes the Figshare
record, while the training implementation stays outside EleD. The package was
marked deprecated in July 2026 in favour of D4CMPP2; the underlying optical data
paper remains the scientific source.

The targets should be named precisely:

| Target | Definition | Unit | Interpretation |
| --- | --- | --- | --- |
| `lambda_abs_nm` | first reported absorption maximum | nm | experimental optical transition; solvent/solid-state dependent |
| `lambda_em_nm` | reported emission maximum | nm | experimental fluorescence/emission observable; solvent and environment dependent |
| `e_abs_ev` | `1239.841984 / lambda_abs_nm` | eV | photon energy corresponding to the absorption maximum |
| `e_em_ev` | `1239.841984 / lambda_em_nm` | eV | photon energy corresponding to the emission maximum |
| `stokes_ev` | `e_abs_ev - e_em_ev` when both exist | eV | derived Stokes shift in energy |
| `xtb_gap_ev` | EleD `LUMO - HOMO` | eV | ground-state xTB orbital gap used as an input descriptor; the optical target is defined by the measured photon energy |

`e_abs_ev` names the photon energy of the absorption maximum. The first allowed
optical transition can differ from the orbital gap because of orbital relaxation,
exciton binding, selection rules, vibronic structure and the environment.
EleD supplies ground-state electronic features for downstream prediction of
experimental optical properties; excited-state and solvent-shifted emission
models belong to the external study.

The current Figshare v3 CSV contains 20,836 rows in the downloaded file, with
17,703 absorption values, 18,851 emission values and 16,422 rows containing
both. It contains repeated chromophore/solvent records from different papers;
61 chromophore/solvent groups have conflicting absorption values and 77 have
conflicting emission values. Treat these as literature measurement records.
Retain the reference DOI, define a duplicate policy before splitting, and report
the number of excluded or aggregated records.

The same chromophore appears in multiple solvents, which is scientifically
valuable and statistically dangerous. A random row split leaks a molecule's
structure into the test set. The primary evaluation should split by canonical
chromophore or scaffold. A secondary solvent-transfer experiment can hold out
solvent conditions for chromophores that remain in training, but it must be
reported as a different task. Solvent structure and environment features must
be included explicitly when that task is studied.

The optical case should compare three feature groups in the external research
workspace:

1. RDKit structure and solvent features;
2. the same features plus EleD core descriptors such as `xtb_homo_ev`,
   `xtb_lumo_ev`, `xtb_gap_ev`, dipole and charge summaries; and
3. the full retained electronic feature archive, with shell, frontier and
   π/reaction descriptors added only when their coverage is sufficient.

The report should use molecule/scaffold-held-out MAE and RMSE in nm and eV,
calibration or interval coverage where uncertainty is available, coverage and
runtime. Emission and absorption should be evaluated separately. The optical
database is therefore a strong first molecular case, while RegioSQM20 remains
the complementary atom-level case.

## Candidates reviewed on 2026-11-13

| Case | Labels and scale | Source and rights | EleD fit | Decision |
| --- | --- | --- | --- | --- |
| **Experimental chromophore optics (Joung DB for chromophore)** | Experimental absorption and emission maxima, lifetime, quantum yield and bandwidths; the current Figshare v3 file contains 20,836 rows with 16,422 rows having both absorption and emission maxima | The database paper and Figshare record define the source; the Deep4Chem/DeepMPP copy is a software-distributed mirror of that record | Excellent for molecular electronic descriptors, solvent-aware features and a direct usability case; optical properties remain experimental excited-state targets | **Primary molecular-property case** |
| **RegioSQM20** | EAS reaction-site labels; 535 molecules and 604 raw labelled atom entries in the supplied `compounds.smiles` export | Code and supplied repository carry an MIT license; the article is CC BY 4.0. The rights of the exported reaction data still need to be recorded before redistribution. | Excellent: atom labels match charges, Fukui, shell, π and frontier descriptors; the input preserves RDKit atom order | **Primary atom-level case** |
| **RegioML** | EAS atom scores; published model and descriptor code, with the full training export hosted in ERDA | GitHub code is MIT; the ERDA data package and its redistribution terms must be checked separately | Very good for a related-work and interoperability check; weak as the first independent dataset because the data archive is external | **Reference / optional replication** |
| **GNN-SOM** | Atom-SOM and bond-SOM labels from KEGG RPAIR; 21,023 molecules in the paper | KEGG RPAIR was downloaded under an academic license; the GitHub repository supplies inference code and weights, while the curated training data remains in the controlled research workspace | Good for testing EleD node and edge archives; license and reaction-centre reconstruction are part of the case protocol | **Second-stage case** |
| **ESNUEL / ESNUEL-ML** | QM-derived methyl cation and methyl anion affinities for 47,921 and 47,440 molecules, respectively; atom-level regression targets | MIT code; data and raw QM archives are distributed through ERDA and must be checked before use | Strong for a large numerical regression and archive-throughput study; expensive and less directly tied to experimental selectivity | **Optional scale case** |
| **SET-LSF** | Minisci/P450 late-stage functionalisation labels; an open literature subset plus a proprietary full set | GitHub code and non-proprietary subset are MIT; the full Pfizer set is restricted | Chemically interesting and atom-labelled, but mixed mechanisms, reaction-condition features and a small open subset complicate the first case | **Future optional case** |

## Why RegioSQM20 is the first case

The RegioSQM20 paper describes 535 literature EAS bromination reactions. It
reports atom-level observed sites and compares quantum and ML approaches. The
supporting repository snapshot reviewed here is
`jensengroup/SI_RegioSQM20@6f9539594b9a07babd76db347ca6d2f550c3ad97`.
Its `RegioSQM20/compounds.smiles` file has one record per molecule in the form
`record_id SMILES site_indices`; comma-separated indices represent multiple
observed sites. Parsing the source snapshot with RDKit 2026.09.1 gives 535
source records and 604 raw positive atom indices. Under EleD's current frozen
allowed-element filter, 533 records and 602 positive labels are selected; the
two excluded records contain Ge and Se. There are no out-of-range labels in the
source export. The associated SDF is a 2-D structure file; the adapter keeps
the SMILES export and its row identity together with the reaction labels.

This is a useful test for EleD because the target has the same granularity as
the feature report. A feature row can retain the original RDKit atom index,
the label, every missing-value reason and the calculation provenance. A single
frozen archive can then feed an atom-ranking table model or a node/edge GNN.
The paper itself used xTB/GFN1, proton affinities, tautomer enumeration and a
solvent model. These calculations form related-work references and optional
ablations; the EleD case uses its own declared descriptor protocol.

The raw export needs explicit handling before calculation:

- preserve comma-separated labels as a set of atom indices;
- verify that labels remain aromatic carbon atoms after SMILES sanitisation;
- retain records with multiple observed sites as a multi-label target;
- record the exact source commit, file hash and the interpretation of its
  labels; and
- keep the repository code license separate from the rights of the underlying
  reaction observations.

RegioSQM20 is still used as a reference point in later reaction-selectivity
and quantum/ML work. For example, the 2024 ESNUEL paper and several recent
site-selectivity studies cite it. This makes the case recognisable while
leaving EleD's contribution clear: a reusable RDKit-indexed electronic feature
archive for downstream models.

## What the other candidates add

**RegioML** is the closest historical comparison. Its public code calculates
GFN1-xTB CM5 charges and a sorted graph-charge-shell vector, then applies
LightGBM. EleD can reproduce the same input category through its public API
while leaving the model and parser in their original research project. The full
data package is linked from ERDA. Confirm its data terms and atom mapping before
freezing it for a primary study.

**GNN-SOM** broadens the target to both atoms and bonds. The paper constructs
the labels from KEGG RPAIR reaction-centre tags and reports 21,023 molecules,
with ten 80/10/10 splits. The GitHub snapshot
`HassounLab/GNN-SOM@dbe88c5aea3cb45f416d4d11d8bc882b72ce3c57` contains the
inference notebook, model states and configuration. The curated training set
is generated and stored in the controlled research workspace. KEGG's academic
license defines the feature-archive release scope.

**ESNUEL** is a useful scale and target-quality study. ESNUEL-ML reports
650,857 atomic MCA values from 47,921 molecules and 534,119 atomic MAA values
from 47,440 molecules, using GFN1-xTB CM5-derived local descriptors and a
higher-level QM target. The public code is MIT. The data and raw QM output are external ERDA archives; the raw calculations are about 29 GB compressed
and 355 GB unpacked. It is appropriate after the EAS archive contract is stable
as a large-scale follow-up case.

**SET-LSF** is a possible radical-reactivity extension. Its public repository
contains an MIT-licensed non-proprietary literature subset and code for
Minisci/P450 site labels, while the full Pfizer set is restricted. Mixed reaction
mechanisms and condition vectors make it a later test of the general electronic
schema. The external study can evaluate EleD
features together with a reaction-condition vector.

## Frozen first-case plan

The external research workspace should freeze the following before any xTB
calculation:

1. the RegioSQM20 source commit and SHA256 of `compounds.smiles`;
2. the parser rule for comma-separated zero-based atom indices;
3. the molecule selection and duplicate policy;
4. the ETKDGv3/conformer and extended-libxtb GFN2/300 K calculation settings;
5. the feature groups and missing-value semantics from `protocol.json`; and
6. a scaffold-plus-Morgan-cluster molecule split, with the test partition
   sealed before model selection.

The first dry run should use 20–50 molecules and report label mapping,
descriptor coverage, failed conformers, runtime and archive size. After that
manifest is frozen, calculate the complete archive and hand it to the external
model project. The external project fits the compact ranking/classification model
and small GNN; model implementations remain in that research workspace.

The case report should include structure-only, charge-only and full EleD feature
groups, together with coverage and cost. It should report failures and negative
results, and present comparisons to GFN1, RegioML and RegioSQM20 as context or
ablation results.

## Source records

- Optical database paper: [Joung *et al.*, Scientific Data (2020)](https://doi.org/10.1038/s41597-020-00634-8), “Experimental database of optical properties of organic compounds”.
- Optical database record: [DB for chromophore, Figshare v3](https://doi.org/10.6084/m9.figshare.12045567.v3), CC BY 4.0.
- JACS Au application paper: [Joung *et al.*, JACS Au (2021)](https://doi.org/10.1021/jacsau.1c00035), “Deep Learning Optical Spectroscopy Based on Experimental Database: Potential Applications to Molecular Design”.
- Deep4Chem/DeepMPP software: [ACCA-KU/DeepMPP](https://github.com/ACCA-KU/DeepMPP). Its bundled optical CSV is a convenience mirror of the frozen Figshare record.

- RegioSQM20 paper: [Ree, Göller and Jensen (2021)](https://doi.org/10.1186/s13321-021-00490-7); article text is CC BY 4.0.
- RegioSQM20 code: [jensengroup/RegioSQM20](https://github.com/jensengroup/RegioSQM20).
- RegioSQM20 supplied data: [jensengroup/SI_RegioSQM20](https://github.com/jensengroup/SI_RegioSQM20).
- RegioML paper: [Ree, Göller and Jensen (2022)](https://doi.org/10.1039/D1DD00032B).
- RegioML code: [jensengroup/RegioML](https://github.com/jensengroup/RegioML).
- GNN-SOM paper: [Porokhin et al. (2023)](https://doi.org/10.1093/bioinformatics/btad089).
- GNN-SOM code: [HassounLab/GNN-SOM](https://github.com/HassounLab/GNN-SOM).
- ESNUEL workflow: [jensengroup/ESNUEL](https://github.com/jensengroup/ESNUEL).
- ESNUEL-ML: [jensengroup/ESNUEL_ML](https://github.com/jensengroup/ESNUEL_ML).
- SET-LSF: [emmaking-smith/SET_LSF_CODE](https://github.com/emmaking-smith/SET_LSF_CODE).
