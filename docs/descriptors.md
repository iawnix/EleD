# Electronic descriptors and density analysis (1.5)

For standardized reports, quality flags, CSV/JSONL export and vertical charge-state
analysis, see [descriptor reports](descriptor-reports.md).

Result-based extraction functions consume an existing `XTBResults` and perform no additional
quantum calculation. Single-call `CalcMolXTBDescriptors` and `RegisterRDKitProperties`
bridge directly to `rdkit.Chem.Descriptors` and `rdMolDescriptors.Properties`.
Units are encoded in field names. Results are JSON-compatible
Python dictionaries/lists, with `None` for unavailable scalar values rather than
fabricated zeros. All calculation paths use libxtb; no tblite or ASE Python dependency is required.

## Interface overview

| Function | Output | Supported input |
| --- | --- | --- |
| `CalcMolXTBDescriptors(mol, confId=-1, ..., missingVal=None, silent=False, cache=None)` | Flat molecular features in one call | Explicit-H RDKit `Mol` with a 3D conformer |
| `RegisterRDKitProperties(*, method="GFN2-xTB", confId=-1, ..., prefix="")` | Tuple of registered numeric property names | Integrates with `rdMolDescriptors.Properties` |
| `CalcXTBDescriptors(result, frontierDepth=2)` | Flat molecular features | Native `XTBResults` or compiled `rdXTBTools.XTBResults` |
| `CalcXTBAtomDescriptors(result)` | One row per atom | Both result paths |
| `CalcXTBBondDescriptors(mol, result)` | One row per existing RDKit bond | Either result type with matching input provenance |
| `CalcXTBFragmentDescriptors(result, fragments)` | Named fragment records, fragment-pair records, unassigned atoms | Both result paths; explicit original atom indices |
| `CalcOrbitalPopulations(result, scheme="lowdin", orbIdx=None)` | Atom × selected-MO populations and electron populations | Same-calculation coefficients, overlap and AO map; verified native basis extension |
| `CalcFrontierDescriptors(mol, result, fragments=None)` | HOMO/LUMO subspace populations, localization and atom-centroid proxies | Same as population analysis, plus matching molecule/conformer |
| `CalcDensityDescriptors(mol, result, mciRings=None)` | Density conservation, shell populations, Mayer/Wiberg, frontier-bond, π and optional MCI descriptors | Same-calculation C/S/maps; `keepOverlap=True` |
| `CalcOrbitalWeightedFukui(result, *, deltaEV=0.5, scheme="lowdin")` | Single-point orbital-weighted Fukui (`f_w_plus`, `f_w_minus`, `f_w_radical`, `dual_w`) | Same-calculation basis (`keepOverlap=True`) |
| `CalcAtomicOrbitalEnergies(result, *, scheme="lowdin", windowEV=2.0)` | Atomic `occ_n_e`, `vir_n`, `occ_e_mean_ev` (negative condensed ALIE proxy), `vir_e_mean_ev`, frontier-weighted energies and PDOS windows | Same-calculation basis (`keepOverlap=True`) |
| `CalcBondOverlapDescriptors(mol, result)` | Per-bond `coop_e`, `ecohp_ev` and Mulliken/CM5 bond polarity (`delta_q_e`, `signed_delta_q_e`, `q_prod_e2`, `*_cm5_*`) | Same-calculation basis (`keepOverlap=True`) + matching molecule |
| `CalcLocalizedOrbitals(mol, result, *, scheme="lowdin", occupationThreshold=1.5, maxIterations=200, tolerance=1e-8)` | Pipek–Mezey LMOs (`lone_pair`, `bond`, `delocalized`), rotated Fock energies, centroids and per-atom lone-pair/bond counts | Closed-shell integer-occupied retained basis (`keepOverlap=True`) |
| `CalcEnergyDecomposition(result)` | Energy terms (`e_elec_ha`…`e_xb_ha`, `g_solv_ha` + sub-terms) and Fermi levels (`fermi_alpha_ha`/`ev`, `fermi_beta_ha`/`ev`) | Native or compiled `XTBResults` from extended properties ABI v1 |
| `CalcMultipoleDescriptors(result, *, mol=None)` | CAMM atomic dipoles/traceless quadrupoles and molecular traceless quadrupole with rotational invariants | GFN2/GFN1 `XTBResults` (GFN0 returns `available=False`) |
| `CalcPolarizabilityDescriptors(result, *, mol=None)` | D4 molecular `mol_pol_bohr3`, `mol_c6_ha_bohr6`, `mol_c8_ha_bohr8` and atomic `alpha_bohr3`, `c6_aa_ha_bohr6` | GFN2/GFN0 `XTBResults` (GFN1 DFT-D3 returns `available=False`) |
| `CalcNuclearESP(result, *, mol=None)` | Exact valence electronic + other-nuclei electrostatic potential at nuclei (`epn_au`, `epn_kcalmol`) | Native or compiled `XTBResults` from extended properties ABI v1 |
| `CalcSolvationDescriptors(mol, solvent, *, solvationModel="alpb", solvationGrid=None, ...)` | Matched gas vs ALPB/GBSA `delta_g_solv_*`, direct `g_solv`, polarization energy, dipole/HOMO/LUMO/charge shifts | Explicit-H 3D RDKit `Mol`; runs 2 single points |
| `CalcVibrationalThermochemistry(mol, *, temperatureK=298.15, pressureAtm=1.0, rotorCutoffCm1=100.0, symmetryNumber=1, ...)` | Harmonic frequencies (cm⁻¹), IR intensities (km/mol), rotational constants and RRHO/qRRHO ZPE/U/H/S/G | Explicit-H 3D RDKit `Mol`; central-difference Hessian of `GetGradient()` |
| `CalcFODFromResult(result, *, scheme="lowdin", electronicTemperature=None)` / `RunFOD(mol, *, electronicTemperature=5000.0, ...)` | Grimme FOD static-correlation index (`n_fod_e`), spin-resolved Fermi levels, orbital FOD weights and atomic/shell `fod_pop_e` (`s/p/d`) | Retained result (`keepOverlap=True` for atomic/shell FOD) or explicit-H 3D `Mol` (1 SCF) |
| `CalcLinearResponseMatrix(result, *, scheme="lowdin", ...)` | Closed-shell independent-particle atom-condensed linear response $\chi_{AB}^{\mathrm{IP}}$ (`chi_e_per_ev`, `chi_e_per_ha`, per-atom softness & nonlocal coupling) | Closed-shell integer-occupied retained basis (`keepOverlap=True`) |
| `CalcBondElectricField(mol, *, stepAngstrom=1e-3, ...)` / `CalcBondElectricFieldFromESP(mol, result, ...)` | Bond-midpoint electric field $\mathbf{E}=-\nabla V$ projected along every RDKit bond (`e_parallel_*`, `e_perp_*`, `e_norm_*`, `e_vec_au`) | Explicit-H 3D `Mol` (1 batched ESP SCF) or matching stencil result |
| `CalcRealSpaceProperties(result)` / `CalcRealSpaceSurfaceDescriptors(mol, result, ...)` / `CalcDensityIsosurface(mol, *, isovalueAU=0.001, ...)` | Exact Cartesian GTO valence density $\rho(\mathbf{r})$, real-space $\mathrm{ALIE}(\mathbf{r})/\mathrm{LEA}(\mathbf{r})$, and `0.001 a.u.` isosurface volume, area, sphericity, ESP & ALIE/LEA | `realspace_extension` ABI v2 (0 SCF for analyzers; 1 SCF for `CalcDensityIsosurface`) |
| `CalcRealSpaceTopology(mol, pts)` / `CalcBondDensityTopology(mol)` / `CalcNCIDescriptors(mol)` | Analytical valence density derivatives, ELF, LOL, RDG, $\mathrm{sign}(\lambda_2)\rho$, bond ellipticity $\varepsilon$, and volumetric NCI/RDG regimes | `realspace_extension` ABI v2 (1 native SCF) |
| `CalcDirectedSurfaceDescriptors(mol, *, isovalueAU=0.001, ...)` | Directed `0.001 a.u.` isosurface extrema for halogen $\sigma$-holes, polar X–H HBD, heteroatom HBA, and aromatic $\pi$-holes | `realspace_extension` + `esp_extension` (1 native SCF) |
| `CalcFragmentInteraction(mol, fragmentA, fragmentB, ...)` | Fixed-geometry fragment EDA ($\Delta E_{\text{elec/rep/es/disp/solv}}$), charge transfer (`ct_*_a_to_b_e`), and inter-fragment WBO/COOP/ECOHP/$\chi_{AB}$ | 3 matched single points (`AB`, `A`, `B`, `keepOverlap=True`) |
| `CalcElectronicCoupling(mol, fragmentA, fragmentB, ...)` / `CalcElectronicCouplingFromResults(...)` | DIPRO/Valeev electronic coupling transfer integrals (`j_hole_ev`/`mev`, `j_electron_ev`/`mev`, `j_ct_ab_ev`/`mev`, `j_ct_ba_ev`/`mev`) invariant to degenerate frontier rotations | 3 matched `keepOverlap=True` single points (0 SCF for `CalcElectronicCouplingFromResults`) |
| `CalcLocalBondForceConstants(mol)` / `CalcLocalBondForceConstantsFromHessian(mol, hess)` | Cremer–Kraka local bond stretching force constants $k^a$ (`mDyn/Å`), Grunenberg compliance constants $C_{ii}$, and adiabatic local frequencies $\omega^a$ | Cartesian Hessian (0 SCF from `CalcVibrationalThermochemistry(..., includeHessian=True)`, or `1 + 6N` gradient evaluations) |

The Python `rdXTBTools.RunMol` wrapper records input provenance from its owned
molecule snapshot. Both Python result types support bond mapping, annotation and
frontier spatial descriptors. Standalone C++ results do not carry the Python
provenance record; molecule-bound helpers reject missing provenance.

## Molecular, atomic and fragment descriptors

```python
from eled import (
    XTBCalculator, CalcMolXTBDescriptors, CalcXTBDescriptors, CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors, CalcXTBFragmentDescriptors, AnnotateMol,
)

# mol must already have explicit hydrogen atoms and a 3D conformer.
result = XTBCalculator(method="GFN2-xTB").Calculate(mol)
scalars = CalcXTBDescriptors(result, frontierDepth=2)
atoms = CalcXTBAtomDescriptors(result)
bonds = CalcXTBBondDescriptors(mol, result)
fragments = CalcXTBFragmentDescriptors(result, {"first": [0], "second": [1, 2]})
annotated = AnnotateMol(mol, result)
# Or in a single step: CalcMolXTBDescriptors(mol, frontierDepth=2)
```

To expose numeric xTB scalar descriptors through RDKit's property registry (sharing one single-point evaluation per conformer across all requested properties):

```python
from rdkit.Chem import rdMolDescriptors
from eled import RegisterRDKitProperties

names = RegisterRDKitProperties(method="GFN2-xTB")
props = rdMolDescriptors.Properties(["xtb_energy_ha", "xtb_gap_ev", "xtb_dip_debye"])
energy_ha, gap_ev, dipole_debye = props.ComputeProperties(mol)
```

`CalcXTBDescriptors` retains the original method/backend/version, total energy
(Hartree), HOMO/LUMO/gap (eV) and electron/orbital counts. It adds:

- `xtb_homo_idx`, `xtb_lumo_idx`: zero-based indices; `None` if absent.
- `xtb_homo_m1_ev`, `xtb_lumo_p1_ev`, etc., up to `frontierDepth`.
  Values beyond the occupied/virtual lists are `None`; depth zero omits them.
- `xtb_dip_x_ebohr`, `xtb_dip_y_ebohr`, `xtb_dip_z_ebohr`, `xtb_dip_ebohr` and
  `xtb_dip_debye`. One e·Bohr = 2.541746473 Debye.
- `xtb_q_sum_e`, `_min_e`, `_max_e`, `_mean_e`, `_sd_e`, `_abs_e`.
  Standard deviation uses all atoms and `ddof=0`, including explicit hydrogens.
- `xtb_wbo_sum`: sum over all `i < j`, including pairs outside RDKit bonds.

HOMO and LUMO follow the aufbau electron count, not an occupation threshold.
With `n_alpha = (N_e + N_unpaired) / 2`, the HOMO is MO `n_alpha - 1` (the SOMO
for open shells) and the LUMO is MO `n_alpha`. xTB smears occupations with a
Fermi distribution, so at 300 K a LUMO already holds more than 1e-7 electrons
when the gap is below about 0.9 eV (about 2.9 eV at 1000 K); a threshold rule
would then mislabel the frontier. `AssessDescriptorQuality` reports the actual
HOMO/LUMO occupations as `frontier_occupations` and raises
`fractional_frontier_occupation` when smearing is significant. COOP/eCOHP use
the full occupation-weighted density over all orbitals. The gap is not an optical
excitation energy; the input unpaired-electron count is not a predicted spin state.
For charged molecules, the dipole and its norm depend on the supplied coordinate
origin. Total energies have method-specific references.

Atom rows contain `atom_index`, `q_e` and `wbo_sum`; `q_e` is the engine charge
set and its exact `charge_model` is recorded in report provenance. For GFN1/2,
`cm5_q_e` is the separately derived CM5 value. `wbo_sum` is the sum of all
off-diagonal WBOs involving that atom and represents a bond-order aggregate.
Bond rows contain `bond_index`, `begin_atom_index`, `end_atom_index` and
`wbo`. They do not create bonds or change RDKit bond orders.

CM5 adds an antisymmetric, distance-dependent pair correction to SCC Mulliken
charges, in elementary-charge units. H/C/N/O pairs use the xTB pair parameters;
other element pairs use the ordered difference of their atomic parameters.
The correction conserves total charge and follows the original atom order under
renumbering. Coincident atom coordinates raise `ValueError`. The corrected
implementation fixes an extra sign reversal for non-tabulated element pairs;
previously archived CM5 values involving such pairs should be regenerated from
their stored geometry and SCC charges. SCC charges and Wiberg bond orders are
unaffected by this correction.

For absent CM5 fields, reports record `cm5_unsupported_method` for GFN0 and
`cm5_not_reported` for GFN1/2. Each populated CM5 array must contain one charge per atom.

Charge arrays follow the calculation atom order. Direct extraction, Fukui analysis and reports check the charge model and atom count. Missing overlap, orbital coefficients or AO-to-atom mapping leaves orbital sections unavailable while the report retains scalar results. Invalid retained matrix dimensions or atom indices raise an input error.

Fragment groups are nonempty **disjoint** collections of unique integer indices
into the calculation's original atom order. Indices include hydrogen atoms;
hydrogens are not added to a fragment implicitly. A complete partition is optional.
Duplicate, negative, Boolean, out-of-range or overlapping indices are rejected.
Each named fragment reports:

- `q_e`, `q_min_e`, `q_max_e`, `q_sd_e`, `num_atoms`;
- `cm5_q_e` and CM5 summary statistics when GFN1/2 CM5 charges are available;
- `wbo_intra`: each pair inside the fragment counted once;
- `wbo_inter`: all pairs to atoms outside, including unassigned atoms.

`pairs` reports each fragment pair's WBO sum once. `unassigned_atom_indices`
identifies omitted atoms. For a complete partition, fragment charges sum to the
molecular charge, and internal plus interfragment WBO sums recover the total
pair sum. No spatial criterion or bond-order cutoff is silently applied.

`AnnotateMol` writes `xtb_wbo_sum` on atoms, in addition to partial charges,
bond WBOs and molecular scalar descriptors. `atomWibergProp` can override its
name, and `privateProps=True` prefixes written property names with `_` following
RDKit's private property convention. Atomic charge and WBO property names must differ.
All validation precedes in-place writes; geometry/topology changes are rejected.

## Canonical selections and multicenter index

`NormalizeSelections(mol, selections)` validates a shared mapping with named
`atoms`, `reactionSites` (or `"reaction_sites"`), `bonds` and ordered `rings`. Original RDKit atom
indices are retained. The same selection object can label electronic reaction
sites and drive local geometry analysis.

`CalcMCI(mol, result, ringAtomIndices, scheme='mulliken')` calculates the
closed-shell Giambiagi–Bultinck multicenter index (`mci`) and ring index
(`i_ring`) for an ordered bonded ring (3–9 atoms) by AO-block matrix
multiplication and trace on `P*S` (Mulliken) or `S^(1/2)*P*S^(1/2)` (Löwdin).
It is strictly invariant to rigid 3D motion, AO orthogonal/sign/permutation
transformations, MO phases, and ring start/direction. Open-shell or
unknown-spin results return `None` with explicit `missing` reasons.
`CalcMCIDescriptors` accepts a named mapping of rings. MCI requires retained
same-calculation basis data. It is an exploratory density descriptor and is not
obtained by averaging `pi_bo`.

## Obtain same-calculation basis data

```python
result = XTBCalculator(method="GFN2-xTB").Calculate(mol, keepOverlap=True)
print(result.aoAtomMap)  # one zero-based atom index per AO
print(result.overlapMatrix)  # AO × AO, same basis/order as coefficients
```

Select the [extended libxtb installation](xtb-patched-build.md) before running
this example. `keepOverlap` is also supported in `XTBOptions` and
`RunXTBBatch`, and is part of cache identity. It defaults to false, accepts
only a Boolean and implies retaining coefficients. Extra storage costs O(NAO²).

Both `XTBResults` and compiled `rdXTBTools.XTBResults` provide
`GetOverlapMatrix()`, `GetAOAtomMap()`, `GetAOShellMap()`, `GetShellAtomMap()`
and `GetShellAngularMomenta()` returning independent NumPy copies. Native tuple
fields are `overlapMatrix`, `aoAtomMap`, `aoShellMap`, `shellAtomMap` and
`shellAngularMomenta`. Maps use zero-based indices in native spherical AO
order; angular momenta use s=0, p=1, d=2. S and C are dimensionless.

The extension captures S and maps in the successful GFN0/1/2 calculation.
Official libxtb lacks this capability and rejects overlap retention before
calculation. No cross-engine overlap substitution or coefficient-square
approximation is used. Radicals are supported with spin-summed occupations;
separate alpha/beta densities and spin populations are not supplied.

## Löwdin and Mulliken populations

```python
from eled import CalcOrbitalPopulations

lowdin = CalcOrbitalPopulations(result)
mulliken = CalcOrbitalPopulations(result, scheme="mulliken")
selected = CalcOrbitalPopulations(result, orbIdx=[3, 4])
```

For real AO coefficients `C`, overlap `S`, and atom A:

- Löwdin amplitudes are `S^(1/2) C`; sum their squares over A's AOs.
- Mulliken populations sum `C_mu,i * (S C)_mu,i` over A's AOs.

The analysis requires finite square matrices, a symmetric positive-definite
well-conditioned overlap, a valid AO map covering every atom, and
`max(abs(C.T S C - I)) <= 1e-6` by default. The measured orthonormality error
is reported. After this check, each MO column is divided by its measured
norm to remove roundoff; large norm defects are rejected rather than hidden.

`atom_pop` has axes `[atom, selected orbital]`, and each column sums
to one. `orb_idx` records its column identities. `mo_n_e`
weights those columns by their actual occupations; selecting only virtual MOs
gives zero electrons despite nonzero orbital populations. Mulliken values can
be negative or exceed one and are never clipped. Löwdin values are nonnegative.
These are basis-dependent population definitions, not unique observables.
The electron populations are not atomic charges without the corresponding
method's reference valence populations. In GFN1/2, Mulliken reconstruction
can be checked against engine charges. GFN0 reports EEQ charges, which are
not defined by these MO Mulliken populations and must not be equated with them.

Missing basis data raises `DescriptorUnavailableError`; malformed or inconsistent
data raises `ValueError`. Descriptor extraction never invokes another engine.

## Frontier subspaces and spatial proxies

```python
from eled import CalcFrontierDescriptors

frontier = CalcFrontierDescriptors(
    mol, result,
    fragments={"donor": donor_atom_indices, "acceptor": acceptor_atom_indices},
    degeneracyEV=1e-3,
)
print(frontier["homo"]["frag_pop"])
print(frontier["lumo"]["pr"])
print(frontier["hl_dist_a"])
```

HOMO/LUMO groups collect contiguous orbitals in the same occupied/virtual class,
within `degeneracyEV` of the frontier energy and `occupationTolerance=1e-6`
of its occupation. Each group's Löwdin atom populations are equally averaged
and normalized to one. Rotating orbitals within that group leaves its descriptors
unchanged; individual orbital populations need not be invariant. Near-degenerate
grouping is a user-controlled descriptor convention, recorded in the output.

Each group returns orbital indices, group size, energy range, atom/fragment
populations, inverse participation ratio `sum(p_A²)`, participation ratio
`1 / sum(p_A²)`, entropy `-sum(p_A ln p_A)` in nats, and population-weighted
atomic coordinate centroid in Angstrom. A fully localized group has participation
ratio 1; equal weights on N atoms give N. Missing frontier groups are `None`.

The two-group summary includes:

- `hl_ovlp = sum_A sqrt(p_H,A * p_L,A)`; a dimensionless
  atom-distribution similarity, not an orbital overlap integral.
- Distance between the two atom-weighted centroids. This is a ground-state
  spatial proxy, **not an excited-state electron–hole separation** and not the
  continuous real-space MO density centroid.

Spatial values require matching original geometry and atom order. Their distance
is invariant to a common translation/rotation of the coordinates. No excitation
energy, oscillator strength or transition dipole is included in this interface.
The separate [report and reactivity APIs](descriptor-reports.md) provide
versioned export and fixed-geometry GFN1/2 ΔSCF/Fukui workflows.

## Validation

Analytical tests cover nonorthogonal metrics, unit conversion, charge and pair
accounting, missing/invalid data, signed Mulliken values, AO permutation, orbital
phase changes, degenerate rotations, coordinate transforms and provenance checks.
Real GFN1/2 tests check 33 frozen conformers per method (66 cases), including
halogens and d shells: CᵀSC, shell/map consistency and population conservation.
Mulliken charge reconstruction is checked where reference valences are known.
GFN0/1/2 on ethanol, ammonium and methyl radical compare compiled and ctypes
matrices/populations. Benzene verifies degenerate frontier subspace invariance.
The C ABI diagnostic also checks copies, disabled retention, invalid buffers,
failed SCF, Hessian invalidation and recovery. No tblite or ASE is imported.

## Shell, frontier-bond and π descriptors

```python
from eled import CalcDensityDescriptors

result = XTBCalculator(method="GFN2-xTB").Calculate(mol, keepOverlap=True)
density = CalcDensityDescriptors(mol, result)
```

`CalcDensityMatrix` exposes `P`, `tr_ps_e` and the electron-conservation
error. `CalcShellPopulations` reports Mulliken and Löwdin s/p/d electron
counts, p/s ratios and shell rows. `CalcFrontierBondDescriptors` uses
near-degenerate HOMO/LUMO projectors and signed Mulliken overlap populations.
`CalcPiDescriptors` reports p-shell projection for aromatic/sp2 atoms and
π-bond order where a local plane exists. Units, invariance and missing-value
semantics are included in each section. See [density-descriptors](density-descriptors.md)
for formulas, units and applicability.

## Orbital-weighted reactivity, ALIE/PDOS, COOP/ECOHP and Pipek–Mezey LMOs

From one retained calculation (`keepOverlap=True`), EleD computes orbital-weighted
reactivity, atomic orbital energy / PDOS statistics, bond overlap populations,
and Pipek–Mezey localized molecular orbitals without additional engine calls:

```python
from eled import (
    XTBCalculator, CalcOrbitalWeightedFukui, CalcAtomicOrbitalEnergies,
    CalcBondOverlapDescriptors, CalcLocalizedOrbitals,
)

result = XTBCalculator(method="GFN2-xTB").Calculate(mol, keepOverlap=True)
owf = CalcOrbitalWeightedFukui(result, deltaEV=0.5, scheme="lowdin")
aoe = CalcAtomicOrbitalEnergies(result, scheme="lowdin", windowEV=2.0)
bond_ovlp = CalcBondOverlapDescriptors(mol, result)
lmo = CalcLocalizedOrbitals(mol, result, scheme="lowdin")
```

- `CalcOrbitalWeightedFukui` evaluates Pino-Rios / Martínez-Araya Gaussian-weighted
  `f_w_plus`, `f_w_minus`, `f_w_radical` and `dual_w` centered at the orbital midpoint `mu_orb_ev` $= (\mathrm{HOMO}+\mathrm{LUMO})/2$ (not the ΔSCF `mu_ev`)
  with exact sum rules ($\sum f_w^+ = \sum f_w^- = 1$, $\sum \Delta f_w = 0$).
- `CalcAtomicOrbitalEnergies` reports per-atom occupied/virtual populations
  (`occ_n_e`, `vir_n`), mean occupied orbital energy `occ_e_mean_ev` (the negative
  of the condensed ALIE proxy: $\mathrm{ALIE}_A = -\bar{\varepsilon}_{\mathrm{occ},A}$),
  `vir_e_mean_ev`, Gaussian frontier-weighted energies, and PDOS windows
  (`pdos_homo_win_e` in e, `pdos_lumo_win` dimensionless).
- `CalcBondOverlapDescriptors` evaluates per-RDKit-bond `coop_e` (positive bonding),
  `ecohp_ev` (negative bonding, eigenvalue-weighted $\mathrm{Tr}(W_{AB}S_{BA})$),
  and engine/CM5 bond polarity (`delta_q_e`, `signed_delta_q_e`, `q_prod_e2`,
  `delta_cm5_q_e`, `signed_delta_cm5_q_e`, `cm5_q_prod_e2`; CM5 entries are `None`
  for GFN0).
- `CalcLocalizedOrbitals` maximizes the Pipek–Mezey functional via 2×2 Jacobi
  rotations (`"lowdin"` or symmetrized `"mulliken"`), classifying orbitals into
  `"lone_pair"`, `"bond"`, and `"delocalized"` with rotated Fock energies
  (`lmo_energy_ev`), centroids (`centroid_a`), and per-atom `n_lone_pairs`,
  `n_localized_bonds`, and `max_lone_pair_energy_ev`. Open-shell, unknown-spin,
  or fractional-occupation results are explicitly withheld (`None` + `missing`).

## Extended `libxtb` properties (energy decomposition, multipoles, polarizability, EPN)

Every successful GFN0/1/2 single point under the extended properties ABI v1
populates energy decomposition, Fermi levels, nuclear ESP, D4 polarizabilities
(GFN2/GFN0), and CAMM multipoles (GFN2/GFN1) on both native `XTBResults` and
compiled `rdXTBTools.XTBResults`:

```python
from eled import (
    XTBCalculator, CalcEnergyDecomposition, CalcMultipoleDescriptors,
    CalcPolarizabilityDescriptors, CalcNuclearESP,
)

result = XTBCalculator(method="GFN2-xTB").Calculate(mol)
decomp = CalcEnergyDecomposition(result)
multipoles = CalcMultipoleDescriptors(result, mol=mol)
pol = CalcPolarizabilityDescriptors(result, mol=mol)
epn = CalcNuclearESP(result, mol=mol)
```

Method-inapplicable quantities are never encoded as zero: `e_aes_ha`/`e_axc_ha`
are `None` for GFN0/1, `e_xb_ha` is `None` for GFN2, `g_solv_ha` is `None` in gas
phase, Fermi levels are `None` when `electronicTemperature <= 0.1 K`, multipoles
return `available=False` for GFN0, and D4 polarizability/C6/C8 returns
`available=False` for GFN1 (which uses DFT-D3).
Finite-field electronic polarizability is outside scope for `libxtb` 6.7.1
without `tblite`/PTB (external electric fields are not supported for GFN0/1/2
and point-charge embedding couples only to monopole shell charges); use
analytical D4 `CalcPolarizabilityDescriptors(result)` instead.

## Implicit solvation (ALPB/GBSA) and vibrational thermochemistry

`XTBOptions`, `XTBCalculator`, `RunXTBBatch`, and `rdXTBTools.RunMol` accept
optional `solvent` (default `None` for gas phase, preserving gas-phase cache
keys), `solvationModel` (`"alpb"` or `"gbsa"`), and `solvationGrid` (Lebedev
grid size). `CalcSolvationDescriptors` runs matched gas-phase and solvated single
points at fixed geometry; `CalcVibrationalThermochemistry` computes the Cartesian
Hessian and dipole derivatives by central finite differences of analytical
`GetGradient()`:

```python
from eled import CalcSolvationDescriptors, CalcVibrationalThermochemistry

solv = CalcSolvationDescriptors(mol, "water", method="GFN2-xTB", solvationModel="alpb")
thermo = CalcVibrationalThermochemistry(
    mol, method="GFN2-xTB", temperatureK=298.15, pressureAtm=1.0, symmetryNumber=1
)
```

`CalcSolvationDescriptors` reports `delta_g_solv_ha`/`kcalmol`, direct `g_solv_direct_ha`,
`solute_polarization_energy_ha`, Born/SASA/HB/shift sub-terms (GFN1/2), dipole
shifts, HOMO/LUMO/gap shifts, and per-atom charge shifts (`solv_delta_q_e`, `solv_delta_cm5_q_e`).
`CalcVibrationalThermochemistry` reports harmonic `frequencies_cm1`, `ir_intensities_km_mol`,
`rotational_constants_cm1`, ZPE, thermal corrections, entropy, and Gibbs free energy
using RRHO (`rotorCutoffCm1=0`) or Grimme qRRHO (`rotorCutoffCm1=100.0` default;
free-rotor reduced moment uses Grimme's fixed `B_av = 1e-44 kg m^2`, as in xtb).
Non-stationary geometries (large gradients or imaginary frequencies) are explicitly
flagged via `is_stationary_point=False` and `warnings`.

## Grimme FOD, condensed linear response, bond electric fields and real-space `0.001 a.u.` isosurfaces

```python
from eled import (
    RunFOD, CalcFODFromResult, CalcLinearResponseMatrix,
    CalcBondElectricField, CalcDensityIsosurface,
    CalcRealSpaceProperties, CalcRealSpaceSurfaceDescriptors,
)

fod = RunFOD(mol, electronicTemperature=5000.0)         # matches xtb --fod --etemp 5000
chi = CalcLinearResponseMatrix(result, scheme="lowdin") # closed-shell chi_AB (e/eV, e/Ha)
ef = CalcBondElectricField(mol, stepAngstrom=1e-3)      # 1 batched ESP SCF for all bonds
iso = CalcDensityIsosurface(mol, isovalueAU=0.001)      # 1 SCF exact 0.001 a.u. isosurface
```

- `RunFOD` and `CalcFODFromResult` implement Grimme & Hansen's spin-resolved Fermi–Dirac
  fractional occupation number weighted electron density (FOD) diagnostic, returning
  `n_fod_e`, `fod_fermi_alpha_ev`/`fod_fermi_beta_ev` (re-solved at the FOD
  temperature; distinct from the engine `fermi_alpha_ev` in `CalcEnergyDecomposition`), `orbital_fod_e`, and per-atom `fod_pop_e`
  with `fod_s_e`/`fod_p_e`/`fod_d_e` shell breakdown when `keepOverlap=True`.
- `CalcReactivityDescriptors` / `RunReactivity` (`protocol_version = "2.2"`) also report
  Gázquez electrodonating/electroaccepting powers (`omega_minus_ev`, `omega_plus_ev`),
  local `omega_minus_loc_ev`/`omega_plus_loc_ev`, multiphilic `delta_omega_loc_ev`, and
  local hyper-softness `s2_loc_1ev2` ($S^2 \Delta f_A$), explicitly returning `None` +
  `missing` whenever hardness $\eta \le 10^{-8}$.
- `CalcLinearResponseMatrix` evaluates the closed-shell independent-particle atom-condensed
  linear response matrix $\chi_{AB}^{\mathrm{IP}} = -4\sum_{i\in\mathrm{occ},a\in\mathrm{vir}} Q^A_{ia} Q^B_{ia}/(\varepsilon_a - \varepsilon_i)$
  with exact electron-conservation row/column sum verification ($\sum_B \chi_{AB} = 0$) and
  explicit withholding (`available=False`, `missing`) for open-shell, unknown-spin,
  fractional-occupation, or degenerate-gap systems.
- `CalcBondElectricField` / `CalcBondElectricFieldFromESP` evaluate $\mathbf{E}(\mathbf{r}_{\mathrm{mid}}) = -\nabla V(\mathbf{r}_{\mathrm{mid}})$
  by 6-point central finite differences of `libxtb` in-memory ESP at every RDKit bond midpoint
  in one batched single point (`e_parallel_au`/`mv_cm`, `e_perp_au`/`mv_cm`, `e_norm_au`/`mv_cm`,
  `e_vec_au`; `1 a.u. = 5142.2067476326 MV/cm`).
- `CalcRealSpaceProperties`, `CalcRealSpaceSurfaceDescriptors`, and `CalcDensityIsosurface`
  use the reviewed `realspace_extension` ABI v2 in `libxtb` to evaluate exact Cartesian GTO
  valence electron density $\rho(\mathbf{r})$, real-space $\mathrm{ALIE}(\mathbf{r})$,
  virtual density, and $\mathrm{LEA}(\mathbf{r})$ without additional SCF calls. Under Fermi
  smearing orbitals are weighted rather than thresholded: ALIE uses occupations $n_i$ and
  the virtual density/LEA use vacancies $2-n_i$; the LEA window starts at the lowest orbital
  with at least one vacancy (aufbau LUMO, or SOMO for open shells). Points below
  `densityCutoffAU` return `None` rather than zero. `CalcDensityIsosurface` locates the exact
  `0.001 a.u.` isosurface in one SCF via Voronoi-partitioned batched radial bisection and
  returns molecular volume, surface area, sphericity, surface grid, and optional isosurface
  Politzer–Murray ESP and ALIE/LEA surface statistics.

## Real-space valence topology (ELF/LOL/ellipticity/NCI), directed σ/π-hole extrema, fragment EDA and Cremer–Kraka local bond force constants

```python
from eled import (
    CalcBondDensityTopology, CalcNCIDescriptors, CalcDirectedSurfaceDescriptors,
    CalcFragmentInteraction, CalcLocalBondForceConstants,
)

topo = CalcBondDensityTopology(mol)                                  # bond rho, ELF, LOL, ellipticity
nci = CalcNCIDescriptors(dimer, fragmentA=[0, 1, 2], fragmentB=[3, 4, 5])
extrema = CalcDirectedSurfaceDescriptors(mol, isovalueAU=0.001)      # sigma/pi-holes & HBD/HBA
eda = CalcFragmentInteraction(dimer, [0, 1, 2], [3, 4, 5])           # 3-SCF EDA + charge transfer
k_local = CalcLocalBondForceConstants(mol, extraPairs=[(1, 3)])      # Cremer-Kraka k^a (mDyn/A)
```

- `CalcRealSpaceTopology`, `CalcBondDensityTopology`, and `CalcNCIDescriptors` use
  `realspace_extension` ABI v2 analytical Cartesian GTO derivatives ($\nabla\psi, \nabla^2\psi$)
  to evaluate $\rho$, $\nabla\rho$, Hessian eigenvalues, $\nabla^2\rho$, kinetic energy density
  $\tau$, Becke–Edgecombe ELF, Schmider–Becke LOL, RDG, $\mathrm{sign}(\lambda_2)\rho$, bond
  ellipticity $\varepsilon = \lambda_1/\lambda_2 - 1$, and volumetric NCI/RDG attractive/vdW/repulsive
  integrals in one SCF. Outputs explicitly note that GFN-xTB uses a valence-only minimal basis
  (no core electrons).
- `CalcDirectedSurfaceDescriptors` locates exact `0.001 a.u.` isosurface extrema along halogen
  $\sigma$-hole axes (and 15° cones), polar `X–H` hydrogen-bond donor axes, heteroatom acceptor
  Voronoi patches, and aromatic ring $\pi$-hole normals in one SCF.
- `CalcFragmentInteraction` decomposes fixed-geometry fragment interactions (`AB - A - B`) across
  electronic, repulsion, electrostatic, dispersion, halogen-bond, and solvation terms alongside
  SCC/CM5/Löwdin/Mulliken charge transfer, inter-fragment WBO/COOP/ECOHP, closed-shell $\chi_{AB}$,
  and interaction dipole shifts.
- `CalcLocalBondForceConstants` and `CalcLocalBondForceConstantsFromHessian` evaluate coordinate-invariant
  Grunenberg compliance constants $C_{ii} = \mathbf{b}_i^T H_{\mathrm{int}}^+ \mathbf{b}_i$,
  Cremer–Kraka local stretching force constants $k^a = 1/C_{ii}$ (`mDyn/Å`), and adiabatic local
  bond frequencies $\omega^a$ (`cm^-1`) for RDKit bonds and optional noncovalent `extraPairs`.
