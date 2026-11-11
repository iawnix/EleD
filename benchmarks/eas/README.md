# EAS benchmark layer

This directory is the model benchmark boundary for the EAS study. It consumes
feature archives produced by `tools/run_eas_features.py` and never changes the
RDKit/xTB descriptor implementation.

`model_protocol.json` freezes the model families, feature groups, preprocessing,
random seeds, metrics and bootstrap acceptance rule. A future training command
must verify both the feature archive's `protocol_sha256` and this model
protocol before fitting anything.

The current utilities are split into three steps:

- `pipeline.py` validates the archive and expands one A–G atom feature group;
- `modeling.py` computes molecule-weighted ranking metrics and paired bootstrap;
- `train.py` optionally uses scikit-learn for training-fold tuning and untouched
  test evaluation.

The command-line entry point is `tools/run_eas_benchmark.py`. It requires the
separate packages in `tools/benchmark-requirements.txt`; the core
`eled` installation does not acquire those machine-learning dependencies.

The benchmark layer is intentionally separate from the public `eled`
calculation API. Descriptor unit tests and numerical invariance tests therefore
remain valid even when the EAS dataset, model, or statistical evaluation is
replaced.

Bond features in groups E/F read their corresponding report sections:

| Feature | Report source | Unit |
| --- | --- | --- |
| `gfn2_wbo` | `bonds[].wbo` | dimensionless |
| `gfn2_mbo` | `mayer_bonds.bonds[].mbo` | dimensionless |
| `gfn2_homo_mul_ovlp_e` | HOMO rows in `frontier_bonds.bonds` | e |
| `gfn2_lumo_mul_ovlp_e` | LUMO rows in `frontier_bonds.bonds` | e |
| `gfn2_pi_bo` | `pi.bonds[].pi_bo` | dimensionless |

Each atom feature is the arithmetic mean over available incident
RDKit-bond/conformer pairs. Failed or error-quality conformers are excluded;
missing pairs are omitted and an empty set yields `None`. Zero values and
frontier overlap signs are retained. HOMO and LUMO are separate columns.
These are unweighted benchmark summaries; they do not use the report's optional
Boltzmann block. The per-conformer descriptor reports retain full coverage data.
