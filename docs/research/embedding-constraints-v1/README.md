# Frozen embedding constraint diagnosis

Recorded 2026-11-05 UTC. See the
[interpretation and limitations](../../embedding-algorithm-research.md) and
[reproduction commands](../../../tools/embedding_research/README.md).

- `uff-manifest.json`, `mmff-manifest.json`: frozen inputs/settings, tool hashes,
  pinned upstream source commit, compiler command and probe binary identity.
- `uff-records.jsonl.gz`, `mmff-records.jsonl.gz`: 486 runs each. Each record
  contains constraints, rejected snapshots, predicate audits, installed-library
  parity, and the diagnostic source-geometry audit.
- `uff-summary.json`, `mmff-summary.json`: counts by group, preset and mode;
  source-predicate failures and optimizer exceptions.
- `references.json.gz`: 81 source geometries in SMILES/explicit-H atom order,
  full query-to-source mappings and independent coordinate/stereo checks.
- `environment.json`: RDKit conda package identities and linked SDK library
  hashes. The probe compiles the embedding translation unit and links the SDK;
  it is not a full rebuild of all upstream code.
- `prior-art.json`: targeted publication metadata and upstream issue status
  checked for this investigation; not a systematic novelty search.
- `sha256.json`: integrity hashes of the other evidence files.

There are 972 successful probe/library parity checks and 2,628 verified snapshot
replays. Source coordinates were loaded only after generation. UFF and MMFF here
refer to the existing embedding-bounds parameterization, not a post-generation
optimization. These are initial-coordinate experiments, not corrected QM8
prediction results or end-to-end `PrepareConformers` recovery rates.
