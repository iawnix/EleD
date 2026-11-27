# Changes

## 0.13.0 — Extended libxtb properties ABI, surface ESP, orbital/LMO analysis, implicit solvation and thermochemistry

- **Tier 1 — ESP surface & single-point orbital/bond descriptors**:
  - Added `ESPSurfaceGrid`, `GenerateESPSurfaceGrid` (exposed scaled RDKit vdW Fibonacci shells with atom ownership and area weights), and `CalcESPSurfaceDescriptors` (Politzer–Murray molecular and per-atom `v_min_au`, `v_max_au`, `v_mean_au`, `v_pos_mean_au`, `v_neg_mean_au`, `v_pos_frac`, `v_neg_frac`, `pi_esp_au`, `sigma2_pos_au2`, `sigma2_neg_au2`, `sigma2_tot_au2`, `nu`, `nu_sigma2_au2`, per-atom `sigma2_au2`, plus optional kcal/mol counterparts).
  - Added `CalcOrbitalWeightedFukui` (Pino-Rios / Martínez-Araya orbital-weighted `f_w_plus`, `f_w_minus`, `f_w_radical`, `dual_w` from one retained calculation), `CalcAtomicOrbitalEnergies` (atomic `occ_n_e`, `vir_n`, `occ_e_mean_ev` = negative condensed ALIE proxy, `vir_e_mean_ev`, frontier-weighted energies, and HOMO/LUMO PDOS windows), and `CalcBondOverlapDescriptors` (bond `coop_e`, `ecohp_ev`, and Mulliken/CM5 polarity descriptors `delta_q_e`, `signed_delta_q_e`, `q_prod_e2`, `delta_cm5_q_e`, `signed_delta_cm5_q_e`, `cm5_q_prod_e2`).
- **Tier 2 — Extended `libxtb` properties C ABI v1 & C++/Python parity**:
  - Added reviewed patch `tools/diagnostics/xtb_properties_export.patch` and ABI verification script `tools/diagnostics/verify_properties.py` exporting `xtb_rdkit_getPropertiesAPIVersion`, `xtb_rdkit_getEnergyDecomposition`, `xtb_rdkit_getPolarizability`, `xtb_rdkit_getMultipoles`, and `xtb_rdkit_getNuclearESP` (exact EPN omitting the singular self-nuclear $Z_A/r_{iA}$ term).
  - Exposed all new properties, NumPy getters, and availability flags symmetrically in C++ `rdXTBTools` (`XTBResults`) and Python `XTBResults` / `native` backend, with full `ToDict`/`FromDict`, `XTBDiskCache`, and spawned-worker serialization parity.
  - Added `CalcEnergyDecomposition`, `CalcMultipoleDescriptors` (rotational-invariant norms, eigenvalues, Frobenius norms, anisotropies), `CalcNuclearESP`, `CalcPolarizabilityDescriptors` (analytical D4 `mol_pol_bohr3`, `mol_c6_ha_bohr6`, `mol_c8_ha_bohr8`, atomic `alpha_bohr3` and `c6_aa_ha_bohr6`), and `CalcFiniteFieldPolarizability` (explicitly raising `DescriptorUnavailableError` documenting `libxtb` 6.7.1 C API scope).
- **Tier 3 — Pipek–Mezey LMOs, ALPB/GBSA solvation and vibrational thermochemistry**:
  - Added `CalcLocalizedOrbitals` (Pipek–Mezey 2×2 Jacobi localization in Löwdin and symmetrized Mulliken schemes with lone-pair/bond/delocalized classification, rotated Fock energies, centroids, and explicit open-shell/fractional withholding).
  - Bound `xtb_setSolvent` via optional `solvent`, `solvationModel` (`"alpb"` / `"gbsa"`), and `solvationGrid` across `XTBOptions`, `XTBCalculator`, `RunXTBBatch`, `XTBCacheKey`, and `rdXTBTools.RunMol` (gas-phase defaults and cache keys preserved when `solvent=None`), and added `CalcSolvationDescriptors`.
  - Added `CalcVibrationalThermochemistry` (central finite-difference Cartesian Hessian and dipole derivatives from analytical `GetGradient()`, harmonic frequencies, double-harmonic IR intensities, RRHO/qRRHO ZPE/U/H/S/G, linear/monoatomic handling, and explicit non-stationary geometry flagging).
  - Updated `DescriptorCatalog()` and `DescriptorStatus()` with all new sections and maturity classifications.
- **API & catalog cleanup**:
  - Removed the hardcoded `DescriptorSensitivity` protocol-comparison helper (`gfn12_delta`, `et300_1000_delta`) and the `diagnostic` pseudo-descriptor section from `DescriptorCatalog()` (now 12 physical/electronic descriptor sections), while retaining `RunConformerAnalysis` (conformer statistics & dispersion) and `AssessDescriptorQuality` (single-point conservation and SCC quality checks).

## 0.12.0 — RDKit interface alignment, single-step wrappers and C++ CM5 parity

Breaking pre-1.0 API alignment:

- Unified all public Python function/method names on RDKit `PascalCase`
  (`XTBCalculator.Calculate`, `GetOptions`, `ToDict`, `FromDict`, `XTBCache.Get`/
  `Set`/`GetOrCompute`/`Clear`, `ReactivityCheckpoint.GetHistory`/`GetReports`,
  `ESPChargeFit.GetCharges`/`ToDict`, `SelectionSet.ToDict`,
  `ConformerGenerationResult.GetMolecule`/`GetReport`,
  `XTBOptimizationResult.GetMolecule`/`GetReport`/`GetResult`,
  `GetInstallationInfo`, `CheckBackend`) and all public parameters, dataclass
  fields and result attributes on RDKit `camelCase` (`confId`,
  `numUnpairedElectrons`, `xtbLibrary`, `keepOverlap`, `espPointsAngstrom`,
  `maxIterations`, `electronicTemperature`, `frontierDepth`,
  `occupationThreshold`, `populationScheme`, `recordId`, `returnErrors`,
  `maxWorkers`, etc.). Serialized report/CSV feature column names remain
  compact unit-suffixed lowercase `snake_case`.
- Renamed the primary Python result model to `XTBResults` (`XTBResult` retained
  as a compatibility alias), matching `rdXTBTools.XTBResults` scalar attribute
  names (`totalEnergyHartree`, `orbitalEnergiesEV`, `orbitalOccupancies`,
  `atomicCharges`, `cm5Charges`, `dipoleAtomicUnits`, `wibergBondOrders`,
  `overlapMatrix`, `aoAtomMap`, `aoShellMap`, `shellAtomMap`,
  `shellAngularMomenta`, `aoPDirections`, `espPointsAngstrom`, `espValuesAU`)
  and `Get*()` owned NumPy copy getters.
- Added `CalcMolXTBDescriptors(mol, confId=-1, ..., missingVal=None, silent=False, cache=None)`
  mirroring `rdkit.Chem.Descriptors.CalcMolDescriptors` and
  `RegisterRDKitProperties(...)` registering numeric xTB scalar descriptors with
  `rdkit.Chem.rdMolDescriptors.Properties` via `Descriptors.PropertyFunctor`
  with shared per-conformer single-point deduplication. Added `privateProps` to
  `AnnotateMol` for RDKit `_`-prefixed private properties.
- Implemented `RDKit::XTBTools::calcCM5Charges` and `XTBResults::cm5Charges`
  directly in C++ and moved numerical settings and ESP coordinate validation into
  C++/NumPy C-API (`Wrap.cpp`), decoupling standalone and Python `rdXTBTools`
  execution from Python helper modules.
- Removed EAS study dataset-splitting utilities from top-level `eled.__all__`
  (available in `eled.eas`), removed the legacy `populationScheme="charges"`
  alias in favor of `"xtb_scc"`, standardized provenance
  `metadata["calculation_settings"]` and `CalcMayerBondDescriptors`
  `missing["mbo"]`, privatized internal helper functions, and bumped
  `DESCRIPTOR_VERSION = "1.4"`, `REACTIVITY_PROTOCOL_VERSION = "2.1"`,
  `STORAGE_VERSION = 2` and `XTBCache` schema `6`.

### Earlier 0.12 development notes

- Expanded the machine-readable catalog to exact shell and reconstructed-charge
  field names and corrected the visualization charge label.
- Kept benchmark protocol and model-training tests out of automatic core CI;
  they now run only through the manually triggered benchmark workflow.
- Propagated π scheme, plane tolerance and explicit reference populations through
  density/report/batch/conformer APIs, with JSON-stable settings and validation
  before engine calls. Fixed benchmark bond extraction to read nested sections
  and keep HOMO/LUMO columns separate while excluding failed conformers.
- Corrected π projections for nonorthogonal bases and stored p-component
  directions with AO rows. π charge uses a π reference population; projected
  bond indices use density products. Undefined planes and inapplicable atoms
  or bonds retain null fields with specific reasons.
- Separated HOMO/LUMO bond statistics and atom/bond π statistics. Added weighted
  local summaries, retained-weight coverage and per-field units. Sensitivity
  comparisons validate provenance and matching valid-conformer subsets.
- Unified frontier grouping, validated standalone bond provenance and basis
  maps, and removed the ambiguous reconstructed-charge alias.
- Corrected patched-header hashes in the pinned xTB source manifest and rebuilt
  the runtime. Added analytical invariance and GFN0/1/2 interface-parity tests.

- Descriptor definitions 1.3 (report schema 2) use compact unit-aware feature names such as
  `xtb_energy_ha`, `q_e`, `wbo`, `eta_ev` and `s_1ev`; frontier and shell
  fields follow the same short naming rules.
- YAeHMOP scalar fields now use the same compact convention (`eht_energy_ev`,
  `eht_n_e`). The descriptor catalog also covers fragment summaries, density
  charge reconstructions and `chi_ev`.
- Conformer analysis keeps equal-weight statistics as the default and can
  additionally report 300 K electronic-energy-weighted means and spreads.
  The weighted block records its weights and identifies them as electronic
  energy weights rather than free-energy populations.
- Renamed the distribution and import package from `rdkit-eth`/`rdkit_eth` to
  `EleD`/`eled`, matching the extended `libxtb` backend that defines
  the project. The C++ extension remains `rdXTBTools`.
- Updated the CLI, environment variables, runtime manifest paths, build
  metadata, examples, tests and bilingual documentation to use the new name.
- Synchronized reviewed xTB patch hashes across the manifest, Python runtime
  and C++ extension, and fixed compiled option validation to use the public
  `CalculationSettings` validator.

## 0.11.0 — isolated execution, SCC diagnostics and validated runtime bundles

- Added spawned libxtb workers for single-point, descriptor and reactivity batches,
  with per-task wall-clock limits, crash isolation, ordered results and parent-owned
  cache/checkpoint writes. Preserve double-precision coordinates across processes.
- Added an optional reviewed SCC diagnostic ABI and C++/Python fields for final
  iterations, energy change, mixer RMS and thresholds; GFN0 remains non-SCF.
- Unified every calculation entry point around the complete reviewed extended
  `libxtb` profile (GFN0 repair, basis export and SCC diagnostics). Upstream
  or partially patched libraries now fail before engine handles are created.
- Added `optimize_xtb()` and `prepare_conformers(..., optimize="xTB")`, using
  native xTB gradients, bounded L-BFGS steps, fixed atoms and auditable
  convergence reports without modifying caller-owned molecules.
- Added an in-process-versus-CLI geometry optimization regression test/tool,
  ran the test in the standard CI profile, and recorded extended xTB
  executable/parameter hashes in runtime bundle manifests.
- Hardened clean source-archive wheel builds so the repository's local `build/`
  directory cannot shadow the Python build frontend.
- Added offline Linux runtime bundling and relocation checks for both interfaces.
  Python compiled calls now discover libxtb in the active runtime by default.
- Refit the historical QM8 study after stereo audit, separately evaluate the checked
  preparation workflow, and compare dipoles on 256 nonoverlapping QM9 molecules.
- Paused embedding-algorithm research; no experimental embedding patch is enabled.

## 0.10.0 — checked geometry preparation

- Added `prepare_conformers` and `ConformerGenerationResult`: copied SMILES/Mol
  input, explicit hydrogens, deterministic preset/seed/initialization retries,
  bounded embedding, per-stage failure counters and explicit partial coverage.
- Check specified tetrahedral and E/Z stereochemistry from coordinates before
  and after optional MMFF94/UFF optimization. Reject unconverged force fields,
  invalid/colliding coordinates and stereo changes. No relaxed-stereo or hidden
  source-coordinate fallback; undefined stereo is flagged or rejected explicitly.
- Preserve the original QM8 prediction protocols; add a separate audit on the
  41 frozen embedding failures and 40 prior-success controls.

## 0.9.1 — inherited native runtime and matched geometry benchmark

- Native library discovery now checks `sys.base_prefix` after `sys.prefix`,
  so a venv inheriting the standard runtime can use its libxtb without manual
  overrides. Explicit library selections retain priority and fail without
  falling back silently.
- Added a matched QM8 source-vs-SMILES geometry study: deterministic ETKDGv3
  embedding, explicit MMFF94 convergence/selection, shared successful cohort,
  frozen train/test IDs and paired test-group bootstrap intervals. No tuning,
  new backend or machine-learning runtime dependency.
- The original pilot's model runner now accepts and validates frozen splits;
  tests cover data leakage, failures, selection and clustered uncertainty.

## 0.9.0 — conformer sensitivity and QM8 pilot

- Added `run_conformer_analysis` and `conformer_descriptor_columns`: immutable
  inputs, actual conformer IDs, per-conformer reports/cache, uniform descriptive
  statistics and explicit field-level missing/failure counts. No Boltzmann
  interpretation or automatic optimization. Charged dipoles and Cartesian
  components are excluded from invariant aggregates.
- Input/protocol/runtime/storage/analysis errors propagate; only calculation
  failures become rows. Missing GFN0 parameters now raise the environment-error
  subclass, consistent with other setup failures.
- Added ETKDG/MMFF sensitivity example and a source-checksummed, ID-joined QM8
  E1-CC2 pilot with fixed random/scaffold splits, three nested feature sets,
  paired predictions, failure accounting and measured cost. sklearn/scipy stay
  in a separate optional benchmark environment, outside package dependencies.

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
  tooling, installed-wheel checks, and `eled doctor` diagnostics.

The compiled wheel depends on the matching local RDKit/Boost and xTB prefixes;
it is not a portable manylinux wheel. Installing official RDKit alone does not
enable this module. No PyPI/conda publication or upstream RDKit merge has been
performed. Excited-state prediction benchmarks remain future work.
