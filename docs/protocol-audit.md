# EAS protocol audit

This note records the executable boundary of `protocol.json`.  It is kept next
to the protocol so a feature archive cannot be mistaken for a completed model
benchmark. Model choices and acceptance rules live separately in
[`benchmarks/eas/model_protocol.json`](../benchmarks/eas/model_protocol.json).

## Checks now enforced

- Dataset manifests carry the canonical SHA256 of `protocol.json`.
- The allowed element set, optional EAS event type, reaction/scaffold
  de-duplication, formal charge, and explicit open-shell spin fields are
  validated before splitting.
- Primary splits keep connected Bemis–Murcko/Morgan leakage groups together.
  Each record also receives a deterministic, leakage-safe five-fold `cv_fold`.
- Feature generation rejects a nonmatching protocol hash and uses exactly ten
  ETKDGv3 conformers at most.  ETKDGv4 and `srETKDGv3` are not silent fallbacks.
- GFN2 optimization is followed by GFN2/300 K, GFN1/300 K, and GFN2/1000 K
  analyses on the same optimized coordinate hashes.  The local Fukui and
  conceptual-DFT reports are generated from fixed-geometry charge states.
- Failed conformers and failed charge states remain explicit records.  No
  failed value is converted to zero.  The output includes conformer, method,
  and electronic-temperature sensitivity diagnostics.

## Items deliberately outside the feature runner

The repository does not yet ship a training/evaluation command for the RegioML
or RegioSQM source data, nor a flattened A–G model matrix. The feature runner
stores the complete RDKit/xTB reports from which that matrix must be constructed.
Nested hyperparameter selection,
molecule-level bootstrap confidence intervals, calibration, and the A–G
ablation tables still require a separate benchmark stage.  Until that stage is
archived, no predictive gain or failure should be inferred from a feature file.

The source adapter also cannot infer an EAS event from an arbitrary SMILES-only
row.  Rows without `event_type`/`reaction_class` are accepted only because the
caller is declaring that the upstream public-data export has already applied
the EAS event filter; this declaration is recorded in the protocol itself.

