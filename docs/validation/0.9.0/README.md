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

Dataset-derived artifacts (QM8/QM9 cohorts, splits, features, predictions,
geometry-preparation audits) were removed from the repository in 0.17.0;
application studies are maintained outside EleD. They remain available in
Git history. `SHA256SUMS.json` still lists the original evidence set.

`installation.json`, `artifact-manifest.json`, installed-wheel reports and test
logs identify the locally tested 0.9 artifacts in `dist/0.9.0`. Evidence is
collected after artifact construction. `SHA256SUMS.json` fingerprints the
evidence files. The prior 0.8 commit's remote CI runs 41878065751, 41878065655
and 41878065315 all completed successfully; these are not 0.9 CI results.
No package-index release or portable Linux/macOS/Windows binary is implied.
