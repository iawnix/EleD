# Changes

## 0.8.0 — persistent reactivity datasets

- Added explicit per-molecule jobs, sequential reactivity batches, successful
  disk caching, exact-dataset checkpoints and auditable failed-row retries.
- JSON/SQLite transactions, checksums and exclusive checkpoint writer locks;
  protocol keys include structures/settings/engine/parameters/implementation.
  Disk single-point caching also integrates with existing batch/report APIs.
- Tightened failure handling: runtime/setup/allocation/extraction and unexpected
  errors propagate even with `run_reactivity(return_errors=True)`. Only genuine
  calculation failures become null descriptor rows; C++ exposes `failureKind`.
- Added frozen 37-molecule corpus, cold/iteration-limit cases, cross-process
  persistence tests and expanded original xTB Fukui comparison. Default physical
  protocol remains 300 K/reference copies, with explicit open-shell spins.

## 0.7.0 — reference wavefunctions for charge-state analysis

- Reactivity now defaults to 300 K, matching single-point options. Explicit 0 K
  remains supported. Protocol version 2.0 records this change; descriptor
  formulas/schema remain unchanged.
- Both ctypes and C++ reuse separate deep copies of the reference wavefunction
  for adjacent states. Added `initialization="independent"` for controlled
  comparisons; no automatic retries or selection between converged solutions.
- Open-shell adjacent spins remain explicit. Native restart handles stay scoped
  to the triplet and are released on success, failures and early termination.
- Added copy/lifetime, order-independence, interface parity and standalone C++
  checks, plus an explicit developer comparison against xTB `--vfukui`.

## 0.6.0 — versioned descriptor reports and vertical charge states

- Unified iteration limit, accuracy and electronic temperature across ctypes and
  C++/Python, with validated ranges, effective temperature and cache identity.
- Added versioned descriptor reports, observable quality flags, ordered batch
  extraction and atomic CSV/JSONL export with explicit missing/error semantics.
- Added fixed-geometry GFN1/2 ΔSCF IP/EA, conceptual-DFT quantities, atomic and
  fragment condensed Fukui descriptors from engine charges or Mulliken/Löwdin
  populations. Explicit spin policy, three-state provenance and SCF failure records.
- Added installed-wheel diagnostics and analytical/real-engine protocol tests.
  Approximate energy differences are not calibrated experimental predictions.
  Existing lower-level descriptor APIs and the 0.5 native engine remain usable.

## 0.5.0 — standard native installation and upstream GFN0 repair

- Pinned upstream PR #1490 at `33fe381`; new engine builds apply only the basis
  export patch. Python/C++ verify the upstream repair identity and binary hash;
  previously verified 0.4 local-patch engines remain usable explicitly.
- Added `tools/install_native.py`: locked RDKit and toolchain prefixes, native
  engine build, clean source-to-wheel build, installation and all-method
  population checks from outside the checkout. The engine shares the runtime
  prefix, so no separate xTB activation is needed for ordinary use.
- Added frozen-geometry comparison across 33 molecules × GFN0/1/2 in separate
  processes, plus upstream identity and installer environment-drift tests.


## 0.4.0 — native basis and descriptors

- Removed the tblite backend, Python extra and current comparison test path;
  historical tools and reports are retained as archived evidence.
- Added an opt-in project-local libxtb C ABI exporting same-calculation overlap
  and AO/shell/atom maps for GFN0/1/2. Build and binary hashes verify capability.
- Added `keep_overlap` / `keepOverlap`, owned matrix/map getters and native
  population/frontier descriptors, including spin-summed radical calculations.
- Compiled Python results now retain input snapshot provenance and support
  bond/annotation/spatial helpers. Official libxtb retains basic GFN1/2 support.
- Migration: remove `backend="tblite"`; build/select the extended libxtb and
  request `keep_overlap=True` for populations. No automatic engine fallback.

## 0.3.0 — electronic descriptors

- Added dipole/charge statistics, neighboring frontier energies, atomic WBO
  sums, provenance-checked bond descriptors and disjoint indexed fragments.
- Common scalar/atomic/fragment extraction accepts native, tblite and compiled
  results; annotation now includes atomic WBO sums.
- Added opt-in same-calculation tblite overlap retention, Löwdin/Mulliken
  populations and degeneracy-aware frontier localization/spatial proxies.
  Native results without overlap fail explicitly; no backend fallback.
- Added analytical invariance/conservation tests, engine-charge reconstruction,
  compiled-result agreement, examples and descriptor installation diagnostics.

## 0.2.0 — in-process calculation backends

Breaking API changes; see [migration notes](docs/migration-0.2.md).

- Removed the xTB CLI backend, executable discovery, output parsers and
  automatic fallback. Native libxtb is the default; tblite stays explicit.
- Removed CLI-only options, process/work-directory result fields and optimized
  conformer mapping. The calculation API supports gas-phase single points.
- Batch execution is sequential, retaining caching and per-molecule errors.
- Native GFN0 still requires a verified patched build; missing libraries and
  unverified GFN0 fail explicitly without launching an external process.
- Updated diagnostics, installation/build tooling, CI and tests for the
  shared-library-only calculation path. Existing 0.1 artifacts remain historical.

## 0.1.0 — local release artifacts, 2026-11-03

- Python CLI/native adapters with explicit units, input/conformer provenance,
  atom/bond annotation, batch execution and in-memory caching.
- Optional RDKit-style C++/Python `rdXTBTools` for gas-phase GFN0/1/2.
  Native GFN0 requires the verified pinned xTB build and D4 patch.
- Explicit optional tblite 0.8.0 Python backend for closed-shell GFN1/2.
- Fixed-conformer numerical comparisons and source-reviewed s/p/d coefficient
  alignment, including orbital signs, subspaces and overlap metrics.
- SHA256-locked Linux/Python 3.11 installation profiles, source-to-wheel build
  tooling, installed-wheel checks, and `rdkit-eth doctor` diagnostics.

The compiled wheel depends on the matching local RDKit/Boost and xTB prefixes;
it is not a portable manylinux wheel. Installing official RDKit alone does not
enable this module. No PyPI/conda publication or upstream RDKit merge has been
performed. Excited-state prediction benchmarks remain future work.
