# Changes

## Unreleased (0.17.0)

- **Repository scope**: removed application research from the EleD repository, as
  datasets, model training and benchmark results are maintained externally:
  `eled.eas` and both copies of the frozen EAS `protocol.json`; EAS/QM8/QM9
  dataset and benchmark tools (`prepare_eas_dataset.py`, `run_eas_features.py`,
  `benchmark_qm8*.py`, `validate_qm9_dipoles.py`, `benchmark-requirements.txt`);
  RDKit embedding-algorithm research (`tools/embedding_research`, `docs/research`,
  `validate_geometry_preparation.py`, `audit_geometry_stereo.py`); research-plan
  documents; and dataset-derived validation artifacts. Historical validation
  READMEs note the removal; the files remain in Git history.
- `AGENTS.md` (local workspace rules) is no longer tracked.
- **Frontier orbitals (breaking, `DESCRIPTOR_VERSION = '1.6'`)**: HOMO/LUMO now
  follow the aufbau electron count (`n_alpha = (N_e + N_unpaired)/2`) instead of
  an occupation threshold. Under xTB Fermi smearing the old `occupation > 1e-7`
  rule labelled the LUMO as HOMO whenever the gap was below about 0.9 eV at 300 K
  (about 2.9 eV at 1000 K), corrupting HOMO/LUMO/gap, frontier groups,
  orbital-weighted Fukui, linear response and frontier bond descriptors.
  `occupationThreshold` is removed from `CalcXTBDescriptors`,
  `CalcMolXTBDescriptors`, `RegisterRDKitProperties`, `CalcFrontierDescriptors`,
  `CalcFrontierBondDescriptors`, `CalcOrbitalWeightedFukui`,
  `CalcAtomicOrbitalEnergies`, `CalcBondOverlapDescriptors`,
  `CalcLinearResponseMatrix`, `CalcElectronicCoupling*`,
  `CalcFragmentInteraction`, `AssessDescriptorQuality` and the report APIs; the
  `xtb_occ_cut` column is removed. `AssessDescriptorQuality` replaces
  `frontier_threshold_sensitive` / `occupied_counts_threshold_window` with
  `fractional_frontier_occupation` / `frontier_occupations`.
- COOP/eCOHP (bond and inter-fragment) use the full occupation-weighted density
  over all orbitals rather than a thresholded occupied subset.
- `CalcDescriptorReport` no longer passes the frontier threshold (1e-7) to
  `CalcLocalizedOrbitals`, whose `occupationThreshold` (default 1.5) selects
  doubly occupied orbitals; reports now use that default.
- **Descriptor definitions (breaking names)**:
  - qRRHO now uses Grimme's fixed `B_av = 1e-44 kg m^2` (as in xtb) instead of the
    molecule's mean principal moment; recorded as `settings.qrrho_b_av_kg_m2`.
  - Same-name fields with different meanings are renamed: FOD Fermi levels
    `fermi_*_ev` -> `fod_fermi_*_ev`; solvation charge shifts `delta_q_e` /
    `delta_cm5_q_e` -> `solv_delta_q_e` / `solv_delta_cm5_q_e`; orbital-weighted
    Fukui `mu_ev` -> `mu_orb_ev`; `atomic_softness_e_per_ev` (-chi_AA, not the
    conceptual-DFT local softness) -> `neg_chi_diag_e_per_ev`.
  - `DescriptorCatalog()` now has one entry per concrete field name (grouped keys
    such as `'alie_ev / lea_ev'` are expanded, with per-field units), and
    `DescriptorStatus()` lists concrete names with no field in two tiers
    (`mul_ovlp_e` is exploratory; `mul_s_e`... are recommended).
  - `xtb_energy_ha` is documented as the total (not electronic) energy.
  - Real-space metadata records the live ABI version instead of a fixed `1`.
  - Docs note that SambVca %Vbur literature values use Bondi radii x1.17.
- **Real-space LEA under Fermi smearing (libxtb rebuild required)**: the realspace
  extension weighted the virtual density/LEA only over orbitals with occupation
  `<= occ_thr`, so a partially occupied LUMO was dropped from the LEA. Acceptor
  orbitals are now weighted by their vacancy `2 - n_i` for every orbital with a
  vacancy above `occ_thr`, and the LEA window starts at the lowest orbital with at
  least one vacancy. The ABI stays v2; the patch hash changes, and the v1/old-v2
  patch acceptance paths are removed from `runtime.py`, `rdXTBTools` and the
  ctypes binding.
- **One interface layer (breaking API)**:
  - `_as_result` is the single normalizer for native and compiled results; it
    copies the compiled extended properties directly (the previous
    `_as_property_result` reconstructed `energyDecompositionFlags` from
    non-`None` fields; `rdXTBTools.XTBResults.energyDecompositionFlags` is now
    exposed). Compiled results carry the same engine identity, version,
    HOMO-LUMO gap and solvation provenance as native ones, and CM5 is no longer
    recomputed as a fallback.
  - Native and compiled results hold the same fields: the gradient is always
    returned and stored once (`gradientHartreePerBohr`; the duplicate
    `metadata['gradient_hartree_per_bohr']` is removed), and MO coefficients are
    retained only with `keepOverlap=True`.
  - `interface='native'|'compiled'` is dispatched in one place
    (`calculator._calculate`) instead of four copies of `rdXTBTools.RunMol`.
  - `xtb_backend` / `xtb_version` are removed from the scalar descriptor table
    (they are provenance), so native and compiled tables are identical.
  - Removed: `XTBOptions.backend`, the `XTBResult` alias,
    `XTBCalculator.options` / `__call__` (use `GetOptions()` / `Calculate()`),
    `XTBResults.fermiEnergyEV` (it returned the HOMO), and the loose
    single-point keyword arguments of `CalcSolvationDescriptors`,
    `CalcVibrationalThermochemistry` and `CalcLocalBondForceConstants`, which
    were silently ignored when `options` was given; pass `options=XTBOptions(...)`.
- **Validate once, compute once**:
  - Single-point settings are validated in one place, `XTBOptions`; the native
    backend builds its options through it instead of repeating the checks, the
    thin `run_xtb_native` wrapper is removed, and `_run_xtb` no longer re-resolves
    charge/spin and reloads the library before the backend does it again.
  - Same-basis analyses share one cached `_orbital_basis(result)` (positive
    definite S, C^T S C check, S^(1/2), P, PS and S^(1/2) P S^(1/2), read-only).
    A full descriptor report previously re-diagonalized S about ten times.
    All analyses use one `ORTHONORMALITY_TOLERANCE = 1e-5` (populations used
    1e-6, the rest 1e-5); `CalcOrbitalPopulations` loses its
    `orthonormalityTolerance` argument.
  - Physical constants live in `eled.constants` (previously redefined in about
    eight modules). Length conversions now use one Bohr radius everywhere: the
    CODATA 2014 value that the ctypes backend and rdXTBTools already used for
    coordinates (`properties` had used CODATA 2018 for Bohr^3 -> Angstrom^3,
    a 4e-10 relative inconsistency).
  - `XTBOptions.ToDict()` always returns every field (the "preserve cache key"
    branch was dead: cache keys already hash all package sources).
  - Removed assertion-only `json.dumps(...)` calls from analysis functions (kept
    at report/export boundaries), discarded validation-only calls, `is None`
    checks on tuple fields, and checks in `AnnotateMol` that repeated
    `_validate_mapping`.

## 0.16.0 — DIPRO electronic coupling (J_AB), full 0.13–0.15 pipeline integration, and catalysis steric complements

- **DIPRO dimer-projection electronic coupling / transfer integral (`eled.dipro`)**:
  - Added `CalcElectronicCouplingFromResults(mol, complexResult, fragmentAResult, fragmentBResult, fragmentA, fragmentB, ...)` (0 SCF calls) and `CalcElectronicCoupling(mol, fragmentA, fragmentB, ...)` (3 matched `keepOverlap=True` single points, `native`/`compiled` parity, zero C ABI changes).
  - Reconstructs the exact complex AO Hamiltonian $H_{\text{AO}} = (S_{AB}C_{AB})\operatorname{diag}(\varepsilon_{AB})(S_{AB}C_{AB})^T$, verifies fragment-to-dimer AO overlap alignment, and performs symmetric Löwdin orthogonalization $H_{\text{eff}} = S_{\text{sub}}^{-1/2} H_{\text{sub}} S_{\text{sub}}^{-1/2}$ for hole (`j_hole_ev`/`mev`), electron (`j_electron_ev`/`mev`), and cross charge-transfer (`j_ct_ab_ev`/`mev`, `j_ct_ba_ev`/`mev`) channels, invariant to rigid 3D motion, MO phase flips, degenerate frontier subspace rotations (e.g. benzene dimers), and monomer swaps.
- **Full `0.13`–`0.15` pipeline integration (`reports`, `conformers`, `visualization`, `diagnostics`)**:
  - Upgraded `CalcDescriptorReport` (`DESCRIPTOR_VERSION = '1.5'`, optional `includeExtended=True`) to populate zero-extra-SCF sections (`energy_decomposition`, `polarizability`, `multipoles`, `nuclear_esp`, `orbital_weighted_fukui`, `atomic_orbital_energies`, `bond_overlap`, `localized_orbitals`, `linear_response`, and sampled `realspace`) identically across `native` and `compiled` `XTBResults`.
  - Added 10 rigid-motion-invariant extended scalar columns (`e_elec_ha`, `e_rep_ha`, `e_es_ha`, `e_disp_ha`, `mol_pol_bohr3`, `mol_c6_ha_bohr6`, `mol_c8_ha_bohr8`, `epn_min_au`, `epn_max_au`, `epn_mean_au`) to `DescriptorColumns()` and `ConformerDescriptorColumns()`, and extended `RunConformerAnalysis` local statistics to track per-atom/per-bond extended sections.
  - Updated `DrawDescriptorMap` (`eled.visualization`) with semantic labels, diverging/sequential palettes, unit rules, and section resolvers for all `0.13`–`0.16` atom and bond descriptors.
  - Extended `CheckBackend` (`eled doctor`) to report `properties_export_verified`, `esp_api_verified`, `realspace_export_verified`, and `realspace_api_version`.
- **Catalysis steric & geometric complements (`eled.local_geometry`)**:
  - Added optional molecularly aligned quadrant and hemisphere breakdown (`zAxisAtom`, `xzPlaneAtom` → `vbur_q1_percent`..`vbur_q4_percent`, `vbur_z_pos_percent`, `vbur_z_neg_percent`) to `CalcBuriedVolume`.
  - Added `CalcPyramidalization(mol, atomIdx, *, neighborIndices=None, confId=-1)` reporting Haddon POAV1 `pyramidalization_p`, `pyramidalization_angle_deg`, `poav_angle_deg`, `angle_sum_deg`, and `pyramid_height_a` (with explicit `None` + `missing` for degenerate geometries).
  - Added `CalcConeAngle(mol, atomIndices, *, center, axisAtom=None, ...)` computing the exact Bilbrey–Allen (2013) minimum enclosing spherical cone angle (or Tolman axial cone angle when `axisAtom` is given), with explicit `apex_inside_atom_sphere` missing handling.
  - Wired optional `pyramidalization` and `coneAngle` requests into `AnalyzeLocalGeometry` and `RunLocalGeometryBatch`.
- **Catalog & exports**: Exported `CalcElectronicCoupling`, `CalcElectronicCouplingFromResults`, `CalcPyramidalization`, and `CalcConeAngle` at package top level and registered all new fields in `DescriptorCatalog()` and `DescriptorStatus()`.

## 0.15.0 — Real-space ABI v2 (ELF/LOL/bond ellipticity/NCI), directed σ/π-hole & H-bond extrema, fragment EDA, and Cremer–Kraka local bond force constants

- **Real-space wavefunction ABI v2 & valence density topology (`eled.realspace`)**:
  - Upgraded `tools/diagnostics/xtb_realspace_export.patch` (`api_version = 2`, backward-compatible with v1) and `tools/diagnostics/verify_realspace.py` with `xtb_rdkit_evalRealSpaceDerivatives`: exact analytical Cartesian GTO first/second derivatives yielding $\rho(\mathbf{r})$, $\nabla\rho(\mathbf{r})$, packed Hessian $\nabla\nabla\rho(\mathbf{r})$, Laplacian $\nabla^2\rho$, positive-definite kinetic energy density $\tau(\mathbf{r})$, Becke–Edgecombe ELF, Schmider–Becke LOL, reduced density gradient RDG, and $\mathrm{sign}(\lambda_2)\rho(\mathbf{r})$.
  - Added `CalcRealSpaceTopology(mol, pointsAngstrom, ...)` (point-level analytical density derivatives, ELF, LOL, RDG, and $\mathrm{sign}(\lambda_2)\rho$), `CalcBondDensityTopology(mol, *, refineAlongBond=True, ...)` (per-RDKit-bond saddle/midpoint valence density topology, ELF, LOL, and bond ellipticity $\varepsilon = \lambda_1/\lambda_2 - 1$), and `CalcNCIDescriptors(mol, ...)` (volumetric NCI/RDG attractive, van der Waals, and steric repulsion descriptors in one SCF, with optional intermolecular `fragmentA`/`fragmentB` interface focus).
- **Directed surface extrema (`eled.surface_extrema`)**:
  - Added `CalcDirectedSurfaceDescriptors(mol, *, isovalueAU=0.001, ...)` locating exact `0.001 a.u.` isosurface extrema in one native SCF along halogen $\sigma$-hole axes (and 15° cones), polar `X–H` hydrogen-bond donor axes (`v_hbd_au`/`kcalmol`), heteroatom acceptor Voronoi surfaces (`v_hba_min_au`/`kcalmol`), and aromatic ring $\pi$-hole normals (`v_pi_max_au`/`kcalmol`).
- **Fragment interaction EDA & charge transfer (`eled.fragment_eda`)**:
  - Added `CalcFragmentInteraction(mol, fragmentA, fragmentB, ...)` running 3 matched fixed-geometry single points (`AB`, `A`, `B`, with `keepOverlap=True` and `native`/`compiled` parity) to report energy decomposition differences (`delta_e_total`, `delta_e_elec`, `delta_e_rep`, `delta_e_es`, `delta_e_aes`, `delta_e_axc`, `delta_e_disp`, `delta_e_xb`, `delta_g_solv`), charge transfer (`ct_scc_a_to_b_e`, `ct_cm5_a_to_b_e`, `ct_lowdin_a_to_b_e`, `ct_mulliken_a_to_b_e`), inter-fragment `wbo_inter`, `coop_inter_e`, `ecohp_inter_ev`, closed-shell `chi_inter_e_per_ev`, and `delta_dipole_debye`.
- **Cremer–Kraka local bond force constants & compliance constants (`eled.local_modes`)**:
  - Added `CalcLocalBondForceConstantsFromHessian(mol, hessianHartreePerBohr2, ...)` (0 SCF calls) and `CalcLocalBondForceConstants(mol, *, extraPairs=(), ...)` computing coordinate-invariant Grunenberg compliance constants $C_{ii} = \mathbf{b}_i^T H_{\mathrm{int}}^+ \mathbf{b}_i$, Cremer–Kraka local stretching force constants $k^a = 1/C_{ii}$ (`mDyn/Angstrom`, `Hartree/Bohr^2`), and adiabatic local frequencies $\omega^a$ (`cm^-1`) for every RDKit bond and optional noncovalent `extraPairs`.
  - Added optional `includeHessian=False` to `CalcVibrationalThermochemistry` so its Hessian can be reused directly without extra gradient evaluations.
- **Catalog**: Registered all new fields across `bond`, `fragment`, and `realspace` in `DescriptorCatalog()` and `DescriptorStatus()`.

## 0.14.0 — Real-space wavefunction ABI (density/ALIE/LEA & 0.001 a.u. isosurface), Grimme FOD, higher-order CDFT, linear response matrix, and bond electric fields

- **Real-space wavefunction ABI v1 & 0.001 a.u. isosurface (`eled.realspace`)**:
  - Added reviewed patch `tools/diagnostics/xtb_realspace_export.patch` and ABI verifier `tools/diagnostics/verify_realspace.py` exporting `xtb_rdkit_getRealSpaceAPIVersion` and `xtb_rdkit_evalRealSpace` (exact BLAS-accelerated Cartesian GTO evaluation of valence density $\rho(\mathbf{r})$ in `e/Bohr^3`, Average Local Ionization Energy $\mathrm{ALIE}(\mathbf{r})$ in `eV`, virtual density in `1/Bohr^3`, and Local Electron Affinity $\mathrm{LEA}(\mathbf{r})$ in `eV` at arbitrary 3D points without extra SCF calls).
  - Exposed symmetric `realSpacePointsAngstrom`, `realSpaceOccupationThreshold`, `realSpaceLEAWindowEV` options and `realSpaceSampled` / `GetRealSpace*()` NumPy getters across `XTBOptions`, `XTBResults`, `native` backend, and C++ `rdXTBTools` (verified by CTest and Python parity tests).
  - Added `CalcRealSpaceProperties(result)` (pure analysis; points below `densityCutoffAU` return `None` rather than zero), `CalcRealSpaceSurfaceDescriptors(mol, result, ...)` (area-weighted molecular and per-atom surface ALIE/LEA statistics), and `CalcDensityIsosurface(mol, *, isovalueAU=0.001, ...)` (Voronoi-partitioned batched radial bisection in **one** native `libxtb` SCF returning molecular volume, surface area, sphericity, `electron_density_isosurface` grid, and optional isosurface ESP and ALIE/LEA surface descriptors).
- **Grimme FOD static-correlation descriptors (`eled.fod`)**:
  - Added `CalcFODFromResult(result, *, scheme="lowdin", electronicTemperature=None)` and `RunFOD(mol, *, electronicTemperature=5000.0, scheme="lowdin", ...)` implementing Grimme & Hansen's spin-resolved Fermi–Dirac smearing and `fodenmak` FOD weights (`n_fod_e`, `fermi_alpha_ev`, `fermi_beta_ev`, `orbital_fod_e`, per-atom `fod_pop_e` and `s/p/d` shell breakdown), matching `xtb --fod --etemp 5000` across both `native` and `compiled` interfaces.
- **Higher-order Conceptual DFT & atom-condensed linear response**:
  - Extended `CalcReactivityDescriptors` and `RunReactivity` (`REACTIVITY_PROTOCOL_VERSION = "2.2"`) with Gázquez electrodonating/electroaccepting powers (`omega_minus_ev`, `omega_plus_ev`), local `omega_minus_loc_ev`, `omega_plus_loc_ev`, multiphilic `delta_omega_loc_ev`, and local hyper-softness `s2_loc_1ev2` ($S^2 \Delta f_A$), explicitly returning `None` + `missing` when hardness $\eta \le 10^{-8}$.
  - Added `CalcLinearResponseMatrix(result, *, scheme="lowdin", ...)` computing the closed-shell independent-particle atom-condensed linear response matrix $\chi_{AB}^{\text{IP}}$ (`chi_e_per_ev`, `chi_e_per_ha`, per-atom `chi_diag_e_per_ev`, `atomic_softness_e_per_ev`, `nonlocal_coupling_e_per_ev`) with exact zero-sum rule verification and explicit withholding for open-shell, unknown-spin, fractional-occupation, or degenerate-gap states.
- **Bond-aligned local electric fields (`eled.electric_field`)**:
  - Added `CalcBondElectricField(mol, *, stepAngstrom=1e-3, ...)` and pure-analysis `CalcBondElectricFieldFromESP(mol, result, ...)` evaluating $\mathbf{E}(\mathbf{r}_{\text{mid}}) = -\nabla V(\mathbf{r}_{\text{mid}})$ via 6-point central finite differences of `libxtb` in-memory ESP at all RDKit bond midpoints in one batched single-point call (`e_parallel_au`/`mv_cm`, `e_perp_au`/`mv_cm`, `e_norm_au`/`mv_cm`, `e_vec_au`; `1 a.u. = 5142.2067476326 MV/cm`).
- **Catalog**: Updated `DescriptorCatalog()` (14 sections: adding `fod` and `realspace`, plus new `atom`, `bond`, and `reactivity` fields) and `DescriptorStatus()`.

## 0.13.0 — Extended libxtb properties ABI, surface ESP, orbital/LMO analysis, implicit solvation and thermochemistry

- **Tier 1 — ESP surface & single-point orbital/bond descriptors**:
  - Added `ESPSurfaceGrid`, `GenerateESPSurfaceGrid` (exposed scaled RDKit vdW Fibonacci shells with atom ownership and area weights), and `CalcESPSurfaceDescriptors` (Politzer–Murray molecular and per-atom `v_min_au`, `v_max_au`, `v_mean_au`, `v_pos_mean_au`, `v_neg_mean_au`, `v_pos_frac`, `v_neg_frac`, `pi_esp_au`, `sigma2_pos_au2`, `sigma2_neg_au2`, `sigma2_tot_au2`, `nu`, `nu_sigma2_au2`, per-atom `sigma2_au2`, plus optional kcal/mol counterparts).
  - Added `CalcOrbitalWeightedFukui` (Pino-Rios / Martínez-Araya orbital-weighted `f_w_plus`, `f_w_minus`, `f_w_radical`, `dual_w` from one retained calculation), `CalcAtomicOrbitalEnergies` (atomic `occ_n_e`, `vir_n`, `occ_e_mean_ev` = negative condensed ALIE proxy, `vir_e_mean_ev`, frontier-weighted energies, and HOMO/LUMO PDOS windows), and `CalcBondOverlapDescriptors` (bond `coop_e`, `ecohp_ev`, and Mulliken/CM5 polarity descriptors `delta_q_e`, `signed_delta_q_e`, `q_prod_e2`, `delta_cm5_q_e`, `signed_delta_cm5_q_e`, `cm5_q_prod_e2`).
- **Tier 2 — Extended `libxtb` properties C ABI v1 & C++/Python parity**:
  - Added reviewed patch `tools/diagnostics/xtb_properties_export.patch` and ABI verification script `tools/diagnostics/verify_properties.py` exporting `xtb_rdkit_getPropertiesAPIVersion`, `xtb_rdkit_getEnergyDecomposition`, `xtb_rdkit_getPolarizability`, `xtb_rdkit_getMultipoles`, and `xtb_rdkit_getNuclearESP` (exact EPN omitting the singular self-nuclear $Z_A/r_{iA}$ term).
  - Exposed all new properties, NumPy getters, and availability flags symmetrically in C++ `rdXTBTools` (`XTBResults`) and Python `XTBResults` / `native` backend, with full `ToDict`/`FromDict`, `XTBDiskCache`, and spawned-worker serialization parity.
  - Added `CalcEnergyDecomposition`, `CalcMultipoleDescriptors` (rotational-invariant norms, eigenvalues, Frobenius norms, anisotropies), `CalcNuclearESP`, and `CalcPolarizabilityDescriptors` (analytical D4 `mol_pol_bohr3`, `mol_c6_ha_bohr6`, `mol_c8_ha_bohr8`, atomic `alpha_bohr3` and `c6_aa_ha_bohr6`; finite-field polarizability remains out of scope for `libxtb` 6.7.1 without `tblite`/PTB).
- **Tier 3 — Pipek–Mezey LMOs, ALPB/GBSA solvation and vibrational thermochemistry**:
  - Added `CalcLocalizedOrbitals` (Pipek–Mezey 2×2 Jacobi localization in Löwdin and symmetrized Mulliken schemes with lone-pair/bond/delocalized classification, rotated Fock energies, centroids, and explicit open-shell/fractional withholding).
  - Bound `xtb_setSolvent` via optional `solvent`, `solvationModel` (`"alpb"` / `"gbsa"`), and `solvationGrid` across `XTBOptions`, `XTBCalculator`, `RunXTBBatch`, `XTBCacheKey`, and `rdXTBTools.RunMol` (gas-phase defaults and cache keys preserved when `solvent=None`), and added `CalcSolvationDescriptors`.
  - Added `CalcVibrationalThermochemistry` (central finite-difference Cartesian Hessian and dipole derivatives from analytical `GetGradient()`, harmonic frequencies, double-harmonic IR intensities, RRHO/qRRHO ZPE/U/H/S/G, linear/monoatomic handling, and explicit non-stationary geometry flagging).
  - Updated `DescriptorCatalog()` and `DescriptorStatus()` with all new sections and maturity classifications.
- **API & catalog cleanup and `CalcMCI` fix**:
  - Fixed `CalcMCI` to use the exact closed-shell Giambiagi–Bultinck AO-block matrix multiplication and trace (`mci` and `i_ring`, strictly invariant to rigid 3D motion, AO orthogonal/sign/permutation transforms, and ring start/direction) and explicitly withhold open-shell/unknown-spin systems (`None` + `missing`).
  - Removed the hardcoded `DescriptorSensitivity` protocol-comparison helper (`gfn12_delta`, `et300_1000_delta`) and the `diagnostic` pseudo-descriptor section from `DescriptorCatalog()` (now 12 physical/electronic descriptor sections), while retaining `RunConformerAnalysis` (conformer statistics & dispersion) and `AssessDescriptorQuality` (single-point conservation and SCC quality checks).
  - Removed `CalcMayerBondDescriptors` (redundant with closed-shell xTB `wbo` and unavailable for open shells), the legacy YAeHMOP/EHT adapter module (`eled.descriptors`: `CalcEHTDescriptors`, `CalcOrbitalDescriptors`), and the always-raising `CalcFiniteFieldPolarizability` stub.

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
