# Validation and verification

EleD validates software contracts, numerical parity and physical invariances
through automated C++, pure-Python and extended-`libxtb` test suites.
Application studies, benchmark datasets and model training are maintained
outside the EleD repository.

## Verification commands

From a Linux x86-64 checkout using the locked prefixes:

1. **Pure-Python / schema / report tests** (`.runtime-env`):

   ```bash
   conda run -p "$PWD/.runtime-env" python -m pytest -q \
     -m "not integration and not benchmark"
   ```

2. **Full extended-engine suite** (`.eled/runtime`, native + compiled `rdXTBTools`):

   ```bash
   ELED_EXPECT_BASIS=1 .eled/runtime/bin/python -m pytest -q -m "not benchmark"
   ```

3. **Standalone C++ core tests** (`.eled/rdxtb-build`):

   ```bash
   ctest --test-dir .eled/rdxtb-build --output-on-failure
   ```

4. **Backend and ABI diagnostics**:

   ```bash
   .eled/runtime/bin/python -m eled doctor --check native --keep-overlap --reactivity
   .eled/runtime/bin/python -m eled doctor --check compiled --keep-overlap --reactivity
   ```

5. **Optional xTB CLI numerical cross-checks** (developer regression tools, not
   runtime backends):

   - `tools/compare_fukui_cli.py`: checks native/compiled vertical Fukui indices
     against `xtb --vfukui` on `tests/data/reactivity_cases.json`.
   - `tools/benchmark_xtb_optimization.py`: checks `OptimizeXTB` geometries and
     energies against the `xtb` CLI optimizer.
   - `tools/validate_reactivity_dataset.py`: checks checkpoint recovery and
     `native`/`compiled` parity across the frozen regression corpus.

## What the test suite checks

- **Engine and ABI parity**: `libxtb` upstream tests, D4 method-switch state
  isolation (upstream PR #1490; see [upstream/xtb-gfn0-upstream-fix.md](upstream/xtb-gfn0-upstream-fix.md)),
  basis export, SCC diagnostics, ESP ABI, extended properties ABI v1, real-space
  ABI v2, and `native` vs compiled `rdXTBTools` parity.
- **Physical and mathematical invariants**: rigid 3D translation/rotation
  invariance, MO phase and degenerate-subspace rotation invariance, AO
  permutation invariance, charge/electron/sum-rule conservation, original RDKit
  atom/bond index preservation, and explicit `None` + `missing` handling
  (undefined quantities are never encoded as zero).
- **Current report schema**: `SCHEMA_VERSION = 2`, `DESCRIPTOR_VERSION = '1.6'`
  (aufbau frontier orbitals), and `REACTIVITY_PROTOCOL_VERSION = '2.2'`.

## Scope and historical records

These checks establish numerical implementation consistency and protocol
behavior; they do not claim experimental prediction accuracy or portability
across unverified xTB releases or operating systems.

Historical pre-0.13 validation logs, wheel-check snapshots and early
investigation archives (`docs/validation/`, `tools/diagnostics/evidence/`,
`tools/diagnostics/archive/tblite-0.3/`) were removed from the working tree to
keep the repository focused on current EleD sources and remain accessible in Git
history.
