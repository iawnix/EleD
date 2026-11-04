# 0.7 reference-wavefunction protocol validation

Validated locally on 2026-11-04. The extended libxtb is unchanged from 0.5/0.6.
The adapter now uses 300 K and reference-wavefunction reuse by default for
reactivity. It records protocol version 2.0; independent initialization and
explicit cold temperatures remain available. Adjacent spins for open-shell
references remain mandatory.

| Check | Result |
| --- | --- |
| Complete extended-engine development suite | 287 passed, 1 skipped |
| Installed wheel, isolated Python outside checkout | 287 passed, 1 skipped |
| Installed wheel with official xTB 6.7.1 | 173 passed, 115 capability-dependent skips |
| Standalone C++ core, including restart triplets and failure recovery | Passed |
| Standard installer upgrade | All six interface/method checks passed; GFN1/2 also test reactivity |
| Fresh-venv compiled and pure wheels | Passed, including reactivity |
| Matched-settings original CLI Fukui comparison | 8/8 matched printed precision |

Tests cover both initialization protocols through native/compiled interfaces,
Mulliken-charge agreement, electron counts, explicit spins, unchanged geometries,
reference-copy identity, independence from neighboring-state order, early-close
cleanup, allocation failures at all three copy steps and recovery after failures.
The C++ triplet owns its handles with RAII and stops at the first failed state.
The shared library's existing `xtb_copyResults` API suffices; no new native patch
or CLI runtime backend was introduced.

`fukui-comparison.json` contains installed-package reports for reference and
independent initialization on water, ethanol, methyl radical and ammonium,
each with GFN1 and GFN2. Settings: 300 K, accuracy 0.1, 250 iterations, fixed
unoptimized coordinates. Adjacent spins explicitly match the CLI choice for
these examples. Exact XYZ input text is retained in `fukui-inputs.json`, and
original stdout/stderr in `cli-logs/`. Executable/library hashes identify the
engine. The comparison tool is developer validation only.

The CLI prints atom Fukui values to **three decimal places**. All cases agreed
within the stated absolute tolerance of 0.00051; the maximum observed absolute
difference was 0.000497389. This establishes agreement at printed precision,
not full-precision identity. Between our two initialization protocols, the
maximum Fukui difference was 7.97e-6 and maximum state-energy difference was
9.62e-12 Hartree on this corpus. These are measured results for these settings,
not guarantees for other molecules or proof of the global SCF minimum.

`installation.json`, `installed-*.json` and `artifact-manifest.json` retain
installed checks and source/wheel fingerprints. `SHA256SUMS.json` hashes the
evidence files. Local artifact copies reside in `dist/0.7.0`; paths are tied to
the matched SDK/runtime prefixes. No package-index release is implied.

Earlier 0.6 temperature/failure records remain unchanged in their own directory.
To reproduce its default protocol, specify both requested temperature 0 K and
`initialization="independent"`. These tests verify numerical/software behavior,
not calibrated IP/EA accuracy or excited-state prediction.
