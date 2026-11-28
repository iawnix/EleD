# 0.10 checked geometry preparation evidence

Validated on 2026-11-05 (Asia/Shanghai). The new preparation layer uses RDKit;
the extended libxtb build and existing electronic descriptor formulas are
unchanged. [API and scope](../../geometry-preparation.md).

| Check | Result |
| --- | --- |
| Full development suite, extended library | 389 passed, 4 skipped |
| Installed compiled wheel, isolated Python outside checkout | 389 passed, 4 skipped |
| Official libxtb 6.7.1 compatibility | 273 passed, 120 skipped |
| Standard installation, native/compiled GFN0/1/2 populations | Passed |
| Installed preparation → SDF/JSON → native xTB example | Passed |
| Preparation on 41 prior embedding failures | 22 accepted, 19 failed |
| Preparation on 40 deterministic prior-success controls | 39 accepted, 1 failed |
| Post-hoc stereo audit of all 983 old generated geometries | 974 pass, 9 rejected |

The four full-suite skips are the same optional/environment cases as 0.9.1.
New geometry tests run without libxtb except for one existing-engine integration
case; that case passes in both native-library profiles. Tests include actual
tetrahedral/E/Z coordinate reconstruction, a reflected conformer that retains
its input tags but must fail, simulated inversion during optimization, rejected
nonconvergence, force-field parameter gaps, partial pools, seed/initialization
retries, unavailable presets, input immutability, unspecified stereo and native
descriptor geometry-hash agreement. Unexpected exceptions propagate.

## Removed preparation audit

Dataset-derived artifacts (QM8/QM9 cohorts, splits, features, predictions,
geometry-preparation audits) were removed from the repository in 0.17.0;
application studies are maintained outside EleD. They remain available in
Git history. `SHA256SUMS.json` still lists the original evidence set.

The preparation layer's contribution is observable acceptance, provenance and
honest failures around RDKit's algorithms. It does not establish that every
SMILES is valid, every molecular geometry is physically plausible, or a returned
conformer is a global minimum or representative equilibrium population.

`installation.json`, `artifact-manifest.json` and wheel checks identify local
artifacts under `dist/0.10.0`; evidence is collected after artifact creation.
`SHA256SUMS.json` fingerprints the archived evidence. All three preceding 0.9.1
CI workflows passed remotely; these local checks are separate from 0.10 CI.
