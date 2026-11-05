# Frozen success-retention and quality checks

See [the acceptance contract, results and reproduction](../../embedding-nonregression.md).
Both direct-replacement prototypes are rejected under the declared quality gate.

- `main-cohort.json`: 40 matching identities outside the 1,024 development cohort,
  plus 40 nonmatching controls. All are selected before generation.
- `main-manifest.json`: fixed settings, executables, source/archive hashes and
  the success/energy gates for the 2,000 paired settings.
- `main-records.jsonl.gz`: complete reports, accepted/optimized coordinates,
  MMFF94 energies, diagnostic timings and flagged-case probe events.
- `main-summary.json`: counts, counterexamples and explicit verdicts.
- `inputs-inputs.json`: three supplementary input strata on the same 40 matching
  molecular identities; unspecified-stereo inputs are derived stress cases.
- `inputs-manifest.json`, `inputs-records.jsonl.gz`, `inputs-summary.json`: the
  600-setting supplementary experiment, with the same metrics and gates.
- `counterexamples.json`: 44 flagged variant/settings comparisons, with energy
  deltas and reproducible settings. They are not 44 different molecules.
- `analysis.json`: combined counts and rejection decision.
- `tests.log`: focused checks, including compiled prototype and gate tests.
- `sha256.json`: hashes of the other distributed evidence files.

The 2,600 paired settings use three executables each: 7,800 runs in total.
All 443 baseline initial successes and 433 baseline post-optimization successes
are preserved. Optimized energy increases produce 16 planarity-only and 28
combined flags; gains on other cases do not offset these. A tolerance of
0.001 kcal/mol was declared before the runs. The largest delta is 56.1756 kcal/mol
and was also independently rerun with identical energies.

The main experiment has only eight baseline successes within the affected
specified-stereo stratum. The specified/MMFF supplementary stratum has none and
is labeled inconclusive for retention. Unspecified input strata provide more
affected baseline successes but pose different stereo-specification tasks.
Reported outcomes do not certify global nonregression or production runtime.

The reference archive supplies molecular identities/SMILES during selection;
source coordinates are never passed to generation. Unchanged original and
patched build identities are recorded. Earlier evidence remains immutable.
