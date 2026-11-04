# 0.8 persistent batch and expanded-corpus evidence

Validated locally on 2026-11-05 (Asia/Shanghai). The native engine and physical
protocol remain unchanged: GFN1/2, 300 K and independent reference copies for
adjacent states by default; open-shell adjacent spins are explicit.

| Check | Result |
| --- | --- |
| Full development suite, extended engine | 316 passed, 1 skipped |
| Installed wheel outside source checkout, isolated Python | 316 passed, 1 skipped |
| Installed wheel with official libxtb 6.7.1 | 202 passed, 115 capability-dependent skips |
| Standalone C++ core | Passed |
| Standard installation and fresh-venv compiled/pure wheel checks | Passed |
| Separate-process batch example rerun | All four records resumed, including input failure |
| Expanded corpus: 40 cases × 2 methods × 2 initializations × 2 interfaces | 320 reports: 304 completed, 16 calculation failures retained |
| Native/compiled corpus comparisons | 160/160 matched, including failure-state identity |
| Original CLI `--vfukui`: 37 molecules × GFN1/2 | 74/74 matched three-decimal output precision |

Persistence tests cover successful cache roundtrips (including matrices and a
new Python process), interruption before the next row commits, exact-dataset
resume, failed-row reuse/retry/history, protocol/geometry/spin/engine/parameter/
implementation isolation, duplicate IDs, invalid inputs, runtime faults,
concurrent writers, checksum corruption and a simulated storage failure that
preserves the previous commit. No native pointer or executable Python object is
serialized. Unexpected errors and storage faults are not relabeled as SCF rows.

The frozen corpus is `tests/data/reactivity_cases.json`; `inputs.json` records its
hash. It includes 33 earlier geometries and four open-shell species (methyl,
triplet oxygen, nitric oxide and phenoxy). Two cold cases and a one-iteration
ethanol case are added per method. All 37 standard structures completed under
both methods/initializations/interfaces. Cold GFN1 ethanol failed attachment;
cold GFN2 ammonium failed ionization; iteration-limited ethanol failed at the
reference. No failure was silently retried or removed.

`corpus-summary.json` records numerical and failure comparisons. Successful
native/compiled state energies and Fukui values were identical in this local
run (maximum measured differences zero). Full reports, including all failed
states/settings/provenance, are in `corpus-reports.jsonl.gz`; each line has a
validation protocol label and its report. These are protocol tests on fixed,
unoptimized structures, not evidence of experimental prediction accuracy.

`fukui-comparison.json.gz` retains all 74 original CLI comparisons and both
initialization reports. Exact XYZ files and original stdout/stderr are in
`fukui-cli-logs.tar.gz`. CLI output has three decimals: maximum absolute
difference was 0.000499884, below the declared 0.00051 tolerance. This does not
assert full-precision equality with the CLI or equivalence on arbitrary inputs.

To read compressed JSON evidence with the standard library:

```python
import gzip, json
with gzip.open("docs/validation/0.8.0/corpus-reports.jsonl.gz", "rt") as f:
    rows = [json.loads(line) for line in f]
```

`installation.json`, `installed-*.json`, logs and `artifact-manifest.json` record
the local build and installed package checks. Full CSV/JSONL exports and SQLite
databases remain in `dist/0.8.0/corpus`; source/wheel artifacts also remain under
`dist/0.8.0`. The CI uploads its own corpus artifacts. Runtime-specific paths
are provenance, not portable installation instructions. All evidence files are
fingerprinted by `SHA256SUMS.json`. No package-index publication is implied.
