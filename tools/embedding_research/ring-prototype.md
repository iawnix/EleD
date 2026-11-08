# Reproduce the local-ring source experiment

The subsequent [nonregression study](../../docs/embedding-nonregression.md)
rejects direct replacement under its quality gate. Build instructions below
reproduce the experimental candidates, not a production recommendation.

Start with the pinned source checkout, SDK and base `embedding_probe` built
using [the diagnostic instructions](README.md). These commands create isolated
executables and generated source files; they do not modify the RDKit checkout or
installed library. Output directories must be empty.

```bash
python tools/embedding_research/build_bounds_variant.py \
  --base-build build/embedding-probe --source build/rdkit-embedding-source \
  --output build/ring-angle-control

python tools/embedding_research/build_bounds_variant.py \
  --base-build build/embedding-probe --source build/rdkit-embedding-source \
  --patch tools/embedding_research/local-ring-angle.patch \
  --output build/ring-angle-candidate

python tools/embedding_research/build_planarity_ablation.py \
  --bounds-build build/ring-angle-control \
  --patch tools/embedding_research/ring-n-planarity-ablation.patch \
  --output build/ring-planarity-only

python tools/embedding_research/build_planarity_ablation.py \
  --bounds-build build/ring-angle-candidate \
  --patch tools/embedding_research/ring-n-planarity-ablation.patch \
  --output build/ring-angle-and-planarity

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build/install-validation/sdk/bin/python tools/embedding_research/evaluate_ring_angle.py \
  --baseline build/embedding-probe/embedding_probe \
  --control build/ring-angle-control/embedding_probe \
  --candidate build/ring-angle-candidate/embedding_probe \
  --output build/local-ring-angle-v1

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build/install-validation/sdk/bin/python tools/embedding_research/evaluate_planarity_ablation.py \
  --initial build/local-ring-angle-v1 \
  --planarity-only build/ring-planarity-only/embedding_probe \
  --combined build/ring-angle-and-planarity/embedding_probe \
  --output build/ring-planarity-ablation-v1
```

Both patches are standard unified diffs against the pinned upstream source.
The planarity patch is applied to the generated instrumented embedding source;
the expected line-number offset comes from the existing logging hooks. Patches
are applied with zero fuzzy context matching. Build manifests retain base build
identity, patch hashes, generated-source hashes, command lines and binary hashes.

The unchanged bounds control must exactly match the original probe across the
1,024-molecule screen. The first experiment rejects unexpected changes outside
the graph-selected target subset. The second reuses and checks the first
experiment's source identities, and checks non-target reports again. Probe runs
have a 120-second subprocess timeout; unexpected failures stop the experiment.
Results are written incrementally, but these scripts do not implement resumption.

Run the optional native regression tests with:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
ELED_RING_ANGLE_CONTROL="$PWD/build/ring-angle-control/embedding_probe" \
ELED_RING_ANGLE_CANDIDATE="$PWD/build/ring-angle-candidate/embedding_probe" \
ELED_RING_ANGLE_COMBINED="$PWD/build/ring-angle-and-planarity/embedding_probe" \
build/install-validation/sdk/bin/python -m pytest -q \
  tests/test_ring_angle_prototype.py tests/test_embedding_constraint_audit.py \
  tests/test_geometry.py
```

The research-native tests skip when the probe variables are absent. This keeps
ordinary package testing independent of a separate RDKit source checkout.
See [the report](../../docs/local-ring-prototype.md) for the negative angle-only
result, factorial comparison and combined-variant regressions.
