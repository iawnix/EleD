# 0.9.1 validation and matched-geometry evidence

Recorded 2026-11-05 (Asia/Shanghai). This release fixes native libxtb discovery
from a venv inheriting the standard runtime, and adds research tooling for a
matched source/generated geometry comparison. The native engine, descriptor
formulas and physical calculation settings are unchanged.

| Check | Result |
| --- | --- |
| Development suite with extended libxtb | 359 passed, 4 skipped |
| Installed compiled wheel, isolated Python outside checkout | 359 passed, 4 skipped |
| Installed wheel with official libxtb 6.7.1 | 243 passed, 120 skipped |
| Optional benchmark and native-discovery tests in research venv | 23 passed |
| Standard installer, native/compiled GFN0/1/2 and populations | Passed |
| Final pure and compiled wheels, fresh isolated installation checks | Passed |
| Inherited venv with library/prefix overrides removed | Native and compiled population checks passed |
| Matched QM8 selected/common-complete cohort | 1,024 / 983; 41 embedding failures retained |
| Source/generated filtered split IDs | Identical |
| RDKit-only baseline predictions between geometry protocols | Identical within 1e-12 eV |

Three full-suite skips are optional sklearn tests, run separately above; the
other is an existing skip. Additional official-library skips reflect missing
capabilities. Tests check parent-prefix discovery and override priority,
deterministic geometry generation, MMFF convergence/selection, explicit failures,
invalid/frozen partitions, train/test group isolation, paired ID/target matching,
cluster bootstrap behavior and insufficient-group intervals. A test changes
only held-out targets and verifies that fitted predictions do not change.

Dataset-derived artifacts (QM8/QM9 cohorts, splits, features, predictions,
geometry-preparation audits) were removed from the repository in 0.17.0;
application studies are maintained outside EleD. They remain available in
Git history. `SHA256SUMS.json` still lists the original evidence set.

`inherited-native.json` and `inherited-compiled.json` were produced using
isolated Python in `build/benchmark-env` with `XTB_LIBRARY`, `XTBPATH`,
`XTB_PREFIX`, `CONDA_PREFIX` and `PYTHONPATH` removed. `installation.json` and
`artifact-manifest.json` identify tested local artifacts; `installed-pure.json`
and `installed-compiled.json` record final-wheel checks. `SHA256SUMS.json`
fingerprints the evidence. Source/wheel artifacts remain in `dist/0.9.1`.

The prior 0.9.0 commit `ae061cf` passed all three remote CI workflows:
41881848194, 41881848192 and 41881848177. These are not 0.9.1 CI results.
