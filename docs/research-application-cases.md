# Application case research

This page records the candidate application cases for EleD. It is a research
decision record, not a model benchmark. Dataset files, model code, trained
weights and numerical results remain in the external research workspace.

## What an EleD case must show

An appropriate case has four properties:

1. labels are attached to atoms or bonds so RDKit identity is part of the
   problem;
2. the source, version and data rights can be recorded;
3. a molecule-level split can prevent close analogues from crossing the test
   boundary; and
4. the same EleD archive can be consumed by a small tabular model and a graph
   model without recalculating the electronic structure.

The result is a demonstration of a usable electronic-descriptor interface. It
does not make EleD a reaction predictor and does not require an EleD feature
group to outperform a published model.

## Why the optical database is a strong EleD case

The optical database is a better molecular-level companion to the atom-level
EAS case. The paper by Joung *et al.* describes 20,236 literature data points,
7,016 chromophores and measurements in solvents or solid state. The current
version-3 Figshare file is the source to freeze for a new study:
`10.6084/m9.figshare.12045567.v3`, licensed CC BY 4.0. Its metadata and file
row count have changed between releases, so the study must record the exact
download URL, file SHA256 and row count rather than relying on the rounded
numbers in a paper or on the copy shipped by DeepMPP.

The DeepMPP repository is a model-training package, not an independent
database authority. Its `_Data/optical.csv` is convenient for a software
example, but the external study should cite and freeze the Figshare record and
keep the training implementation outside EleD. The package itself was marked
deprecated in July 2026 in favour of D4CMPP2; that software status does not
invalidate the underlying optical data paper.

The targets should be named precisely:

| Target | Definition | Unit | Interpretation |
| --- | --- | --- | --- |
| `lambda_abs_nm` | first reported absorption maximum | nm | experimental optical transition; solvent/solid-state dependent |
| `lambda_em_nm` | reported emission maximum | nm | experimental fluorescence/emission observable; solvent and environment dependent |
| `e_abs_ev` | `1239.841984 / lambda_abs_nm` | eV | photon energy corresponding to the absorption maximum |
| `e_em_ev` | `1239.841984 / lambda_em_nm` | eV | photon energy corresponding to the emission maximum |
| `stokes_ev` | `e_abs_ev - e_em_ev` when both exist | eV | derived Stokes shift in energy |
| `xtb_gap_ev` | EleD `LUMO - HOMO` | eV | ground-state xTB orbital gap, an input descriptor rather than the experimental optical gap |

`e_abs_ev` must not be labelled as a HOMO–LUMO gap. The first allowed optical
transition can differ from the orbital gap because of orbital relaxation,
exciton binding, selection rules, vibronic structure and the environment.
Likewise, EleD currently does not solve an excited state or a solvent-shifted
emission state. The scientifically useful claim is therefore that EleD makes
ground-state electronic features easy to calculate and archive for downstream
prediction of experimental optical properties.

The current Figshare v3 CSV contains 20,836 rows in the downloaded file, with
17,703 absorption values, 18,851 emission values and 16,422 rows containing
both. It contains repeated chromophore/solvent records from different papers;
61 chromophore/solvent groups have conflicting absorption values and 77 have
conflicting emission values. These are literature measurement records, not
duplicates that may be silently averaged. The external study must retain the
reference DOI, define a duplicate policy before splitting, and report the
number of excluded or aggregated records.

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
| **Experimental chromophore optics (Deep4Chem/Joung DB)** | Experimental absorption and emission maxima, lifetime, quantum yield and bandwidths; the current Figshare v3 file contains 20,836 rows with 16,422 rows having both absorption and emission maxima | Original database article and the Figshare record are CC BY 4.0; the DeepMPP copy is a software-distributed mirror and must not replace the frozen source record | Excellent for molecular electronic descriptors, solvent-aware features and a direct usability case; optical properties remain experimental excited-state targets | **Primary molecular-property case** |
| **RegioSQM20** | EAS reaction-site labels; 535 molecules and 604 raw labelled atom entries in the supplied `compounds.smiles` export | Code and supplied repository carry an MIT license; the article is CC BY 4.0. The rights of the exported reaction data still need to be recorded before redistribution. | Excellent: atom labels match charges, Fukui, shell, π and frontier descriptors; the input preserves RDKit atom order | **Primary atom-level case** |
| **RegioML** | EAS atom scores; published model and descriptor code, with the full training export hosted in ERDA | GitHub code is MIT; the ERDA data package and its redistribution terms must be checked separately | Very good for a related-work and interoperability check; weak as the first independent dataset because the data archive is external | **Reference / optional replication** |
| **GNN-SOM** | Atom-SOM and bond-SOM labels from KEGG RPAIR; 21,023 molecules in the paper | KEGG RPAIR was downloaded under an academic license; the GitHub repository contains inference code and weights, not the curated training data | Good for testing EleD node and edge archives, but the license and reaction-centre reconstruction add risk | **Second-stage case** |
| **ESNUEL / ESNUEL-ML** | QM-derived methyl cation and methyl anion affinities for 47,921 and 47,440 molecules, respectively; atom-level regression targets | MIT code; data and raw QM archives are distributed through ERDA and must be checked before use | Strong for a large numerical regression and archive-throughput study; expensive and less directly tied to experimental selectivity | **Optional scale case** |
| **SET-LSF** | Minisci/P450 late-stage functionalisation labels; an open literature subset plus a proprietary full set | GitHub code and non-proprietary subset are MIT; the full Pfizer set is restricted | Chemically interesting and atom-labelled, but mixed mechanisms, reaction-condition features and a small open subset complicate the first case | **Future optional case** |

## Why RegioSQM20 is the first case

The RegioSQM20 paper describes 535 literature EAS bromination reactions. It
reports atom-level observed sites and compares quantum and ML approaches. The
supporting repository snapshot reviewed here is
`jensengroup/SI_RegioSQM20@6f9539594b9a07babd76db347ca6d2f550c3ad97`.
Its `RegioSQM20/compounds.smiles` file has one record per molecule in the form
`record_id SMILES site_indices`; comma-separated indices represent multiple
observed sites. Parsing with RDKit 2026.09.1 gives 535 valid records, 604 raw
positive atom indices, and no out-of-range labels. The associated SDF is a 2-D
structure file and does not contain the reaction labels, so the adapter must
keep the SMILES export and its row identity together.

This is a useful test for EleD because the target has the same granularity as
the feature report. A feature row can retain the original RDKit atom index,
the label, every missing-value reason and the calculation provenance. A single
frozen archive can then feed an atom-ranking table model or a node/edge GNN.
The paper itself used xTB/GFN1, proton affinities, tautomer enumeration and a
solvent model. Those calculations are related work and optional ablations;
they are not a requirement that EleD reproduce the RegioSQM20 predictor.

The raw export needs explicit handling before calculation:

- preserve comma-separated labels as a set of atom indices;
- verify that labels remain aromatic carbon atoms after SMILES sanitisation;
- retain records with multiple observed sites rather than turning them into a
  single-class target;
- record the exact source commit, file hash and the interpretation of its
  labels; and
- keep the repository code license separate from the rights of the underlying
  reaction observations.

RegioSQM20 is still used as a reference point in later reaction-selectivity
and quantum/ML work. For example, the 2024 ESNUEL paper and several recent
site-selectivity studies cite it. This makes the case recognisable while
leaving EleD's contribution clear: a reusable RDKit-indexed electronic feature
archive, not another RegioSQM implementation.

## What the other candidates add

**RegioML** is the closest historical comparison. Its public code calculates
GFN1-xTB CM5 charges and a sorted graph-charge-shell vector, then applies
LightGBM. EleD can reproduce the *kind* of input through its public API without
copying the model or its parser. The full data package is linked from ERDA,
so it should not be the first frozen source until its data terms and atom
mapping are confirmed.

**GNN-SOM** broadens the target to both atoms and bonds. The paper constructs
the labels from KEGG RPAIR reaction-centre tags and reports 21,023 molecules,
with ten 80/10/10 splits. The GitHub snapshot
`HassounLab/GNN-SOM@dbe88c5aea3cb45f416d4d11d8bc882b72ce3c57` contains the
inference notebook, model states and configuration, but not the curated
training set. KEGG's academic license means the feature archive should be
generated in a controlled research workspace and released only when the terms
allow it.

**ESNUEL** is a useful scale and target-quality study. ESNUEL-ML reports
650,857 atomic MCA values from 47,921 molecules and 534,119 atomic MAA values
from 47,440 molecules, using GFN1-xTB CM5-derived local descriptors and a
higher-level QM target. The public code is MIT, but the data and raw QM output
are external ERDA archives; the raw calculations are about 29 GB compressed
and 355 GB unpacked. It is therefore appropriate after the EAS archive
contract is stable, not as the first end-to-end case.

**SET-LSF** is a possible radical-reactivity extension. Its public repository
contains an MIT-licensed non-proprietary literature subset and code for
Minisci/P450 site labels, while the full Pfizer set is restricted. Mixed
reaction mechanisms and condition vectors make it a poor first test of a
general electronic descriptor schema, but it could later test whether EleD
features remain useful when a reaction-condition vector is supplied by the
external study.

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
model project. The external project may fit a compact ranking/classification
model and a small GNN, but EleD will not contain either implementation.

The case report should include structure-only, charge-only and full EleD
feature groups, together with coverage and cost. It should report failures and
negative results, and should treat any comparison to GFN1/RegioML/RegioSQM20 as
context or ablation rather than an acceptance criterion for EleD.

## Source records

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
