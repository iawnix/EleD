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

## Separate preparation audit

`tools/validate_geometry_preparation.py` uses the frozen 41 failure IDs from
0.9.1 plus the first 40 prior-success controls in the original target-blind
selection order. It requests **one** conformer with seed 42 under the new
preparation protocol: ETKDGv3 → ETKDGv4 → srETKDGv3 when available, two seeds
for each initialization mode, 300 internal embedding iterations, MMFF94 up to
500 iterations, and coordinate/stereo validation before and after optimization.
The first attempt with an accepted pool ends retries. This differs from the
old eight-conformer benchmark protocol; it is not a controlled comparison of
individual embedding presets and is not a new prediction benchmark.

Of the 41 old failures, 22 now return accepted optimized geometries: eight first
succeed under ETKDGv3/random initialization, twelve under ETKDGv4/standard
initialization, and two under srETKDGv3/standard initialization. Nineteen remain
failures. Of 40 prior-success controls, 39 pass and ID **2226** fails the new
stereo checks. No chirality constraint is disabled. Across all attempts, 33
candidate conformers are rejected for stereo mismatch; this is a candidate
count, not a unique-molecule count. RDKit internal failure counters likewise
count internal attempts/stages, not distinct failed molecules.

`preparation-manifest.json` records exact IDs, input and implementation hashes.
`preparation-reports.jsonl.gz` records every attempt, raw/named failure counters,
candidate checks, optimization outcomes, final coordinates and atom ordering.
`preparation-summary.json` records counts and IDs. The run took 100.5 s locally;
timing is not a portability or performance guarantee.

## Correction to interpretation of the historical geometry benchmark

The failed prior-success control motivated a post-hoc audit of **all 983**
generated geometries saved in 0.9.1. `tools/audit_geometry_stereo.py` reconstructs
the original explicit-H atom ordering, applies the saved coordinates, removes
stereo tags on a temporary copy and reassigns stereo from the geometry.
Nine coordinate sets do not reproduce a specified input tetrahedral center.
IDs and mismatched atom indices are in `prior-stereo-audit.json` together with
input/script/implementation hashes. No E/Z mismatch was found in this audit.

The original procedure enforced chirality during embedding but did not perform
this independent check after force-field optimization. The saved final
coordinates alone do not establish which step introduced each discrepancy.
Original geometries, predictions, scores and checksums remain unchanged in
0.9.1 for reproducibility. Their published table **includes these nine rows**
and is not a fully stereo-validated benchmark; this limitation is now stated in
the main study document. No corrected MAE or renewed accuracy claim is made here.

The preparation layer's contribution is observable acceptance, provenance and
honest failures around RDKit's algorithms. It does not establish that every
SMILES is valid, every molecular geometry is physically plausible, or a returned
conformer is a global minimum or representative equilibrium population.

`installation.json`, `artifact-manifest.json` and wheel checks identify local
artifacts under `dist/0.10.0`; evidence is collected after artifact creation.
`SHA256SUMS.json` fingerprints the archived evidence. All three preceding 0.9.1
CI workflows passed remotely; these local checks are separate from 0.10 CI.
