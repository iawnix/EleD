# 0.6 descriptor protocol and installed-package evidence

Validated locally on 2026-11-04 (Asia/Shanghai), Linux x86-64 / Python 3.11 /
RDKit 2026.09.1. The native engine is the same reviewed source and basis extension
as 0.5; these changes concern adapter settings and descriptor workflows.

| Check | Result |
| --- | --- |
| Extended libxtb, complete development suite | 274 passed, 1 official-only skip |
| Official 6.7.1 compatibility, compiled adapter | 172 passed, 103 capability-dependent skips |
| Installed compiled wheel, isolated Python outside checkout | 274 passed, 1 official-only skip |
| Standalone C++ core | Passed, including settings and clamping rejection |
| Standard installer incremental upgrade | Passed all six method/interface checks; GFN1/2 also tested vertical reactivity |
| Fresh-venv compiled wheel | Both interfaces, basis populations, reports and reactivity passed |
| Fresh-venv pure wheel with official library | Basic report and charge-based reactivity passed |

`python-*.log`, `cpp-core.log`, `installation.json` and `installed-*.json` retain
these checks. `artifact-manifest.json` identifies the source and both wheels by
SHA256. Artifact paths refer to this local installation. The copied artifacts
also reside in `dist/0.6.0`; neither portable binaries nor package-index releases
are implied. `SHA256SUMS.json` fingerprints the evidence records.

The analytical tests check energy-difference signs, the explicit softness
convention, negative EA, unstable reciprocal hardness, Fukui sums, mixed
geometry/settings/engine rejection and failure records. Real-engine tests check
native/compiled parity, charged/open-shell species, Mulliken-versus-charge
differences, Löwdin sums, thermal occupation flags, cache separation and recovery
after failed SCF. A molecule with two valid but different conformers is rejected
as a vertical charge-state triplet.

`sensitivity-inputs.json` freezes the exact unoptimized test geometries.
`sensitivity-reports.json` contains nine installed-package calculations:
water/GFN2, ethanol/GFN1 and ammonium/GFN2 at requested temperatures 0, 300 and
1000 K. Effective cold temperature is 1e-6 K. Seven completed and two cold
protocols failed; failures are retained with null descriptors. Every successful
row keeps three state energies, charge/spin, engine/parameter identity and
settings. This small study is a protocol regression example, not a curated
experimental benchmark.

`scf-probe.json` separately records both interfaces at 250 and 1000 iterations
for those cold problematic structures, plus explicit warm protocols. Increasing
iterations did not fix the cold calculations. Warm protocols converged, but their
energy differences change; they must not be substituted silently. The first
failed state is a diagnostic, not a universal invariant across numerical builds.

These checks establish software/protocol behavior. They do not establish
experimental IP/EA accuracy, bound anions, lowest-energy spin states or
excited-state predictive performance. Conformer ensembles and calibrated
prediction benchmarks remain outside this release.
