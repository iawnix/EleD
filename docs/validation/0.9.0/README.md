# 0.9 conformer sensitivity and QM8 pilot evidence

Recorded locally on 2026-11-04 UTC. The extended libxtb source pin and existing
single-point/reactivity physical protocols are unchanged. These additions use
native libxtb; no CLI, ASE or tblite runtime dependency is introduced.

| Check | Result |
| --- | --- |
| Full development suite | 340 passed, 3 skipped |
| Installed compiled wheel, isolated Python outside checkout | 340 passed, 3 skipped |
| Installed wheel with official libxtb 6.7.1 | 224 passed, 119 skipped |
| Optional benchmark protocol tests in isolated sklearn environment | 3 passed |
| Standard installation, native/compiled GFN0/1/2 and populations | Passed |
| Pure/compiled wheels installed into fresh venvs | Passed |
| Standalone C++ interface | Passed |
| Installed conformer example vs source example | Identical statistics |
| QM8 pilot selected/common complete cohort | 1,024/1,024; zero calculation or feature failures |

Two full-suite skips are the sklearn group-split tests, run successfully in the
separate benchmark environment; the remaining skip is an existing capability
case. Official-engine skips additionally reflect unavailable extended-basis and
native GFN0 functionality. No tests are made conditional on prediction gains.

`conformer-summary.json` and `conformer-reports.jsonl.gz` retain 11 conformers
across ethanol, ethylene glycol, butane, phenethylamine and benzene. Exact
coordinates, ETKDG seed, MMFF outcomes, selected conformer IDs, xTB provenance,
individual descriptor reports and field coverage are included. These are
uniform descriptive statistics, not Boltzmann averages. Single-conformer zero
spreads are not evidence of stability across unsampled geometries.

The [QM8 protocol and interpretation](../../qm8-pilot.md) describes the
target, source, filtering, feature lists, fixed model and caveats. Evidence:

- `qm8-manifest.json.gz`: archive/script/source hashes, all exclusions, selected
  IDs, chemical identities, target values, geometry hashes, settings and versions.
- `qm8-splits.json.gz`: exact training/test IDs for all six fixed splits.
- `qm8-features.jsonl.gz`: numerical features, quality flags, per-molecule
  engine/parameter provenance and measured cost for the common cohort.
- `qm8-predictions.json.gz`: every test target and prediction for all 18 fits.
- `qm8-summary.json`: MAE/RMSE, sizes, feature counts, software versions and
  timings. `qm8-failures.json` records calculation/feature failures (empty here).

MAE averaged over the three random splits is 0.3469 → 0.2085 → 0.2028 eV
for RDKit, +basic xTB, +orbital distributions. For scaffold splits it is
0.3473 → 0.2105 → 0.2093 eV. Basic xTB helps this pilot; additional orbital
features provide small, inconsistent incremental gains. No significance,
state-of-the-art or general scientific accuracy claim is made. Source
geometries are provided, and their preparation cost is not included.

JSON evidence can be opened using standard-library `gzip.open(path, 'rt')`
and `json.load`, or one `json.loads(line)` per JSONL line. The full QM8 download
stays outside Git in `build/qm8`; raw source structures are not bundled here.
Scripts refuse to overwrite an existing benchmark output directory and verify
the input archive checksum before reading it.

`installation.json`, `artifact-manifest.json`, installed-wheel reports and test
logs identify the locally tested 0.9 artifacts in `dist/0.9.0`. Evidence is
collected after artifact construction. `SHA256SUMS.json` fingerprints the
evidence files. The prior 0.8 commit's remote CI runs 41878065751, 41878065655
and 41878065315 all completed successfully; these are not 0.9 CI results.
No package-index release or portable Linux/macOS/Windows binary is implied.
