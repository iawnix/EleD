# Isolated parallel libxtb batches (0.11)

`RunXTBBatch`, `RunDescriptorBatch` and `RunReactivityBatch` accept
`max_workers=1`, `timeout=None` and `isolate=False`. The default keeps the existing
in-process behavior. More than one worker, a timeout, or `isolate=True` selects
fresh spawned Python processes, each directly calling the same libxtb library.
No xTB executable is searched for or invoked.

```python
from eled import RunDescriptorBatch, XTBOptions

def main():
    # molecules must have explicit H atoms and prepared 3D coordinates.
    reports = RunDescriptorBatch(
        molecules, options=XTBOptions(keep_overlap=True),
        max_workers=4, timeout=60.0,
    )

if __name__ == '__main__':
    main()
```

Use an importable script and a `__main__` guard, as required by Python's spawn
mode. Interactive environments that cannot import the main module should run
the script instead. Before launching Python, set `OMP_NUM_THREADS=1` and
`OPENBLAS_NUM_THREADS=1` to avoid multiplying numerical-library threads across
workers. The adapter does not change global thread settings behind the caller's
back. Each worker starts from a fresh interpreter and exits after one molecule
(one complete charge-state triplet for reactivity). This favors isolation over
startup throughput; small molecules may be faster with the default sequential API.

## Execution and failures

The timeout is a finite positive number of seconds. It covers process startup,
calculation and result serialization, beginning when a task starts a worker;
queue time is excluded. It is not a separately enforced timer on each charge
state. Expired workers are terminated, then killed if necessary and joined.
Native crashes or an exit without a result produce `XTBWorkerError`; timeouts
produce `XTBTimeoutError`. Both derive from `XTBRunError`.

`RunXTBBatch(return_exceptions=True)` retains these exceptions in input order.
Otherwise it raises and cleans up other active workers. Descriptor batches
produce explicit error rows with null values. Reactivity batches record worker
failure rows and preserve checkpoint attempt history. A killed triplet has no
recoverable partial-state diagnostics; it restarts from the reference on retry.
Configuration, installation and storage failures in the isolated path abort;
they are not converted to apparent SCF failures.

Input molecules are snapshotted with double-precision coordinates, properties,
atom ordering and actual conformer IDs preserved. No native pointers cross
process boundaries. Workers temporarily serialize owned results into private
`eled-worker-*` directories in Python's temporary directory. Normal exit,
timeout and handled interruption remove them; forcibly killing the entire parent
can leave such directories. These are transport files, not xTB calculation logs.
Persistent storage still uses checksummed JSON, not pickle.

## Cache and checkpoint behavior

The parent owns caches and writes checkpoints as individual workers complete;
the returned list retains input order. Only successful results enter caches.
Simultaneous identical misses may calculate twice. Cache hits start no worker.
Native calculation settings and input hashes determine numerical identity;
worker counts and timeouts do not change that identity. Cached elapsed times
describe the original calculation.

Reactivity checkpoint reuse keeps failed rows unless `retry_failed=True`.
Consequently, increasing a timeout alone does not retry old failures:

```python
reports = RunReactivityBatch(
    jobs, checkpoint='dataset.sqlite', cache=cache,
    max_workers=4, timeout=120, retry_failed=True,
)
```

Failures, retries and completed results remain auditable. Inputs/reports are still
materialized in memory; this is a local batch executor, not a streaming or
distributed scheduler. No speedup guarantee is made.
