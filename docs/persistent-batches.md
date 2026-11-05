# Persistent reactivity datasets (0.8)

`run_reactivity_batch` calculates fixed-geometry triplets with explicit record IDs,
per-molecule spins/options/fragments, successful-result caching and resumable
checkpoints. It uses the existing libxtb protocol (300 K and reference copies by
default). It does not optimize structures or search spin states.

```python
from rdkit_eth import (
    ReactivityJob, XTBOptions, XTBDiskCache, ReactivityCheckpoint,
    run_reactivity_batch, export_descriptor_reports,
)

# Supply explicit-H molecules with finite 3D conformers.
jobs = [
    ReactivityJob("water", water),
    ReactivityJob("radical", methyl_radical, ionized_uhf=0, attached_uhf=0),
    ReactivityJob("custom", another_mol, options=XTBOptions(accuracy=0.1),
                  fragments={"first_atom": [0]}),
]
cache = XTBDiskCache("results/cache.sqlite")
checkpoint = ReactivityCheckpoint("results/dataset.sqlite")
rows = run_reactivity_batch(jobs, cache=cache, checkpoint=checkpoint)
export_descriptor_reports(rows, "results/descriptors.csv")
export_descriptor_reports(rows, "results/descriptors.jsonl")
```

Call the same function with the same jobs and checkpoint after interruption.
IDs must be nonempty, unique strings; output order follows input order. A job's
`options` overrides the batch `options` object as a whole, rather than merging
individual fields. `population_scheme`, `interface` and `initialization` are
batch-wide settings with the same meanings as `run_reactivity`. Open-shell jobs
must specify both adjacent spin counts. Every input molecule is snapshotted.

## Cache versus checkpoint

The **cache stores only successful results** and can be shared across datasets.
Record IDs are excluded from cache identity, so the same calculation can serve
different IDs. Only complete reactivity triplets are cached; failed reports and
partial wavefunctions are not cached. `XTBDiskCache` also works with
`run_xtb_batch` and `run_descriptor_batch`, preserving tuple/matrix structure in
normalized single-point results. No pickle, raw native handles or new Python
dependency is used.

The **checkpoint records the exact ordered dataset and every committed attempt**,
including failures. A row commits atomically after the complete triplet succeeds
or reports a calculation failure. Interruption during a triplet restarts that
molecule from its reference. A successful cache write just before interruption
can still be reused even if its checkpoint row was not committed.

Completed checkpoint rows are reused by default, including failed rows. Use
`retry_failed=True` to retry failures explicitly under the same protocol.
`checkpoint.history(record_id)` returns all saved attempts; replacing a failed
current row with a successful retry preserves the earlier diagnostic.
`checkpoint.reports()` returns completed rows in dataset order. Interrupted
attempts with no completed report are not counted as committed attempts.

Changing IDs, order, geometry or any calculation protocol rejects the existing
checkpoint. For additions, corrected inputs or new settings, use a **new checkpoint
path** and retain the shared cache; unchanged successful calculations can still
be reused. Cache and checkpoint must use separate files.

## Protocol identity

Cache keys and checkpoint manifests fingerprint:

- indexed topology, geometry, conformer, molecular charge and reference spin;
- GFN method, iteration limit, accuracy, temperature and basis-retention setting;
- adjacent spins, initialization policy, population scheme and fragment indices;
- libxtb binary and parameter identity, relevant environment/thread settings;
- Python implementation sources, compiled adapter binary when used, Python
  major/minor, RDKit and NumPy versions;
- storage/cache, descriptor/schema and reactivity protocol versions.

The checks are deliberately conservative: changing paths/environment or rebuilding
the compiled module can invalidate reuse even if numerical results would match.
Resume requires a working compatible engine; cache files are not an offline
substitute for runtime validation. Restart Python after changing code/libraries
and keep the environment fixed during a run. Open/read already exported reports
directly for offline analysis.

`provenance.batch_protocol` records the fingerprinted recipe. `execution` contains
`source` (`calculation`, `cache`, `checkpoint` or `input_validation`), `key` and
committed attempt number. CSV includes `execution_json` as well as the existing
provenance/settings/error columns; JSONL preserves all nested state/atom data.

## Error policy and storage guarantees

All job configurations and runtime capabilities are checked before calculations.
Invalid molecular geometry (including missing conformers, unexpanded hydrogen
counts or nonfinite coordinates) and `mol=None` produce input-error rows with
null descriptors. `None` can represent an upstream parser failure; the caller
should retain its original source record separately. Actual xTB single-point
failures retain their three-state diagnostics and null derived descriptors.

Missing open-shell spins, invalid spin parity, fragment definitions, duplicate
IDs, unsupported methods, missing/incompatible libraries or basis capabilities
raise. Unexpected programming/analysis errors also raise. Since 0.8,
`run_reactivity(return_errors=True)` no longer suppresses arbitrary exceptions:
allocation/setup/extraction failures raise `XTBEnvironmentError`, while genuine
single-point failures retain `ChargeStateCalculationError` reporting. C++ results
expose `failureKind` for that distinction; rebuild compiled consumers on upgrade.

SQLite commits protect completed rows against interruption. JSON payloads carry
SHA256 checksums; corruption, incompatible storage schemas and I/O failures abort
instead of silently becoming cache misses or molecular failures. Checksums detect accidental record damage.

Checkpoint runs hold a nonblocking OS file lock: another writer to the same
checkpoint fails clearly. Shared-cache writes use SQLite transactions; concurrent
identical misses may calculate twice. The default is sequential; [isolated workers](isolated-batches.md) support
parallel tasks and wall-clock limits. Cache/checkpoint writes stay in the parent
process on a local filesystem; this is not a distributed scheduler. Disk capacity, retention and eviction
are managed by the caller. The batch materializes inputs and returned reports;
large datasets can be split into separately checkpointed chunks sharing a cache.

## Reproduce examples and validation

```bash
python examples/reactivity_batch.py --output build/my-dataset
python examples/reactivity_batch.py --output build/my-dataset  # resumes
python tools/validate_reactivity_dataset.py --output build/corpus-check
python tools/compare_fukui_cli.py \
  --executable .rdkit-eth/runtime/bin/xtb \
  --library .rdkit-eth/runtime/lib/libxtb.so \
  --cases tests/data/reactivity_cases.json --output build/cli-check
```

The frozen corpus includes 37 molecules spanning common functional groups,
cations/anions, halogens, phosphorus/sulfur compounds and four open-shell species.
The dataset validator adds two cold cases and one iteration-limit case per method,
then checks GFN1/2, both initializations and both interfaces (320 triplets),
including resume without calculation. Original CLI comparison covers the 74
standard molecule/method pairs and is limited by its three-decimal Fukui output.
All failures remain in dataset accounting. Unoptimized geometries and successful
API agreement do not establish experimental accuracy or excited-state validity.
