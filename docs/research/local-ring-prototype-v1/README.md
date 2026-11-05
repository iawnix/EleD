# Frozen local-ring prototype evidence

See [the report](../../local-ring-prototype.md) and
[reproduction commands](../../../tools/embedding_research/ring-prototype.md).

- `angle-manifest.json`: first experiment, fixed cohort and graph-selected targets,
  parameter grid, source hashes and all three executable build manifests.
- `angle-screen.jsonl.gz`: 1,024 paired baseline/angle-only outcomes and unmodified
  bounds-control parity. Full events are retained for the three targets.
- `angle-sweep.jsonl.gz`: 405 target runs (135 each for baseline, angle-only and
  existing global MMFF bounds), including failure snapshots and post-MMFF checks.
- `angle-summary.json.gz`: first-experiment aggregates and diagnostic source audits.
- `ablation-manifest.json`: source and build provenance plus hashes linking the
  first experiment to the targeted planarity-only/combined ablation.
- `ablation-screen.jsonl.gz`: 1,024 paired planarity-only/combined screen records.
- `ablation-sweep.jsonl.gz`: 270 additional target runs (135 per ablation).
- `ablation-summary.json`: paired screen and target-sweep counts.
- `analysis.json`: aggregate and per-preset counts; includes the six combined
  losses versus planarity-only, as well as 38 gains.
- `tests.log`: focused verification, including compiled prototype tests.
- `sha256.json`: integrity hashes for the other evidence files.

Trial counts refer to repeated runs of **three development molecules**. They are
not independent molecular samples. The broader 1,024 screen uses one fixed seed
and initial-coordinate validation only; its `post_mmff_valid_converged: 0`
counters mean post-MMFF was not evaluated there. Post-MMFF is evaluated only in
the target sweeps. An existing MMFF embedding-bounds run is distinct from this
post-embedding MMFF94 optimization.

Source coordinates enter neither executable. Source audits in the first
experiment use the immutable mappings in
`../embedding-constraints-v1/references.json.gz` only after generation completes.
The ablation experiment never reads source coordinates. These files do not
replace the historical QM8 geometry/prediction evidence.

The angle summary is compressed for review size. The ablation manifest records
the SHA256 of its original uncompressed JSON; decompression preserves those exact
bytes. `sha256.json` records the distributed compressed file.
