# EAS benchmark layer

This directory is the model benchmark boundary for the EAS study. It consumes
feature archives produced by `tools/run_eas_features.py` and never changes the
RDKit/xTB descriptor implementation.

`model_protocol.json` freezes the model families, feature groups, preprocessing,
random seeds, metrics and bootstrap acceptance rule. A future training command
must verify both the feature archive's `protocol_sha256` and this model
protocol before fitting anything.

The benchmark layer is intentionally separate from the public `rdkit_eth`
calculation API. Descriptor unit tests and numerical invariance tests therefore
remain valid even when the EAS dataset, model, or statistical evaluation is
replaced.

