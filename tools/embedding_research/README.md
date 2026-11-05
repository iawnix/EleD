# Embedding algorithm research tools

These tools observe RDKit's embedding internals. They do not implement a new
conformer generator and are not part of the installed `rdkit_eth` API.
The initial findings and proposed experiments are in
[the research report](../../docs/embedding-algorithm-research.md).

The probe builds one executable from the **unmodified, pinned** upstream
`Embedder.cpp`, inserting read-only logging hooks into a generated copy. It links
against an existing RDKit 2026.09.1 SDK. Installed libraries are not replaced.
The other RDKit implementation files, including the bounds builder and force
fields, come from that SDK. This is not a complete source rebuild of RDKit.

## Reproduce

Run from the repository root. Paths below match the project's native installation
validation SDK; substitute a compatible SDK and C++20 compiler if necessary.

```bash
git clone --depth 1 --branch Release_2026_09_1 --filter=blob:none --sparse \
  https://github.com/rdkit/rdkit.git build/rdkit-embedding-source
git -C build/rdkit-embedding-source sparse-checkout set \
  Code/GraphMol/DistGeomHelpers Code/DistGeom Code/ForceField/CrystalFF

python tools/embedding_research/build_probe.py \
  --source build/rdkit-embedding-source \
  --sdk build/install-validation/sdk \
  --compiler build/install-validation/sdk/bin/x86_64-conda-linux-gnu-g++ \
  --output build/embedding-probe

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build/install-validation/sdk/bin/python tools/embedding_research/run_diagnosis.py \
  --probe build/embedding-probe/embedding_probe \
  --archive build/qm8/gdb8.tar.gz \
  --force-field UFF --output build/embedding-diagnosis-uff

OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
build/install-validation/sdk/bin/python tools/embedding_research/run_diagnosis.py \
  --probe build/embedding-probe/embedding_probe \
  --archive build/qm8/gdb8.tar.gz \
  --force-field MMFF --output build/embedding-diagnosis-mmff
```

Output directories must be empty. The QM8 archive is the same checksum-verified
input used by `tools/benchmark_qm8.py`; it is never passed to the probe.
The pinned source commit is `fece8caa860bdf6c9c82bdb9253f22a0cf87150c`.

For each molecule/preset/mode, the runner verifies that probe and installed RDKit
return the same result, failure counters, and successful coordinates (absolute
tolerance `1e-10` Angstrom). It also compares replayed chiral and tetrahedral
predicates with direct calls to upstream C++ checks on captured snapshots.
A known `bad direction in linearSearch` exception is explicitly recorded and
checked against the installed library; other unexpected exceptions stop the run.

After **all** generation completes, source geometries are read, mapped to the
SMILES atom order with a full stereochemical graph match, and audited against
the captured constraints. No source geometry is used as a seed or template.

`records.jsonl.gz` contains actual distance bounds, atom-indexed chiral sets,
matched torsion atom quadruplets, up to three rejected snapshots per failure
type, the accepted snapshot if present, and predicate replays. `INITIAL_COORDS`
failures have counters but no snapshot: failed eigendecompositions need not
contain usable coordinates. Snapshot positions retain a fourth coordinate when
present; geometric audits use the first three, as upstream checks do.

This first probe does **not** record every optimizer step, torsional energies,
gradients, or a minimal unsatisfiable constraint subset. Aggregate failure counts
cover all attempts; saved snapshots are a capped diagnostic sample and must not
be treated as an unbiased frequency sample. Generation checkpoints are written
incrementally for debugging, but this script does not resume from them.
