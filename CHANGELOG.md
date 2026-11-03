# Changes

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
