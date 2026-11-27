# Descriptor catalog

The public descriptor API has five layers. The complete machine-readable
catalog is returned by `eled.DescriptorCatalog()`; this page states which
layers are appropriate for validation.

Public descriptor functions and class methods follow RDKit's `PascalCase`
(`Calc...`, `Run...`, `Calculate`, `Get*`), and Python parameters and result
attributes use `camelCase` (`confId`, `keepOverlap`, `totalEnergyHartree`).
Serialized report/CSV feature column names use compact unit-suffixed lowercase
`snake_case`, for example `xtb_gap_ev`, `xtb_energy_ha`, `xtb_dip_ebohr` and
`s_1ev`. There is one public spelling for each field. `DescriptorCatalog()` expands grouped
definitions into exact fields such as `mul_s_e`, `low_d_e` and `low_qrec_e`.

The naming rules keep feature columns short while retaining their meaning:

- `ha`, `ev`, `e`, `ebohr`, `debye` and `a` identify Hartree, eV, electron,
  e·Bohr, Debye and Ångström values;
- `q`, `wbo`, `mbo`, `mul`, `low`, `eta`, `s` and `omega` identify charge,
  Wiberg, Mayer, Mulliken, Löwdin, hardness, softness and electrophilicity;
- `idx`, `n`, `sd`, `ipr`, `pr` and `ent` mean index, count, standard deviation,
  inverse participation ratio, participation ratio and entropy;
- `xtb_` marks molecular fields directly reported by xTB. Atom and bond rows
  use their shorter local names (`q_e`, `cm5_q_e`, `wbo`, `mul_p_e`, `low_ps`).
  `q_e` is the engine charge set, with its exact model recorded in
  `charge_model`; `cm5_q_e` is the geometry-corrected CM5 field.

| Layer | Main API | What it means | Engine calls during extraction |
| --- | --- | --- | --- |
| Molecular | `CalcXTBDescriptors`, `CalcMolXTBDescriptors`, `RegisterRDKitProperties` | Energy, frontier energies, dipole, charge and Wiberg summaries | 0 for `CalcXTBDescriptors`; 1 single point for single-step wrappers |
| Atom/bond topology | `CalcXTBAtomDescriptors`, `CalcXTBBondDescriptors` | Engine charges, incident Wiberg sums and RDKit-bond Wiberg orders | 0 |
| Orbital density | `CalcOrbitalPopulations`, `CalcShellPopulations`, `CalcDensityMatrix` | Populations derived from the same `C`, `S`, AO and shell maps | 0; requires `keepOverlap=True` |
| Frontier/bond | `CalcFrontierDescriptors`, `CalcFrontierBondDescriptors`, `CalcMayerBondDescriptors`, `CalcBondOverlapDescriptors` | Degeneracy-invariant HOMO/LUMO localization, Mayer orders, COOP/ECOHP and bond polarity | 0; requires `keepOverlap=True` |
| Orbital & LMO analysis | `CalcOrbitalWeightedFukui`, `CalcAtomicOrbitalEnergies`, `CalcLocalizedOrbitals` | Single-point orbital-weighted Fukui, condensed ALIE/PDOS and Pipek–Mezey LMOs | 0; requires `keepOverlap=True` |
| ESP surface & EPN | `GenerateESPSurfaceGrid`, `CalcESPSurfaceDescriptors`, `CalcNuclearESP` | Politzer–Murray vdW-shell ESP statistics and exact nuclear ESP (EPN) | 0 (`CalcESPSurfaceDescriptors` requires sampled ESP) |
| Extended properties | `CalcEnergyDecomposition`, `CalcMultipoleDescriptors`, `CalcPolarizabilityDescriptors` | Energy decomposition, Fermi levels, CAMM multipoles and D4 polarizability/C6/C8 | 0; properties ABI v1 |
| Solvation & thermochemistry | `CalcSolvationDescriptors`, `CalcVibrationalThermochemistry` | Matched ALPB/GBSA solvation shifts and finite-difference RRHO/qRRHO thermochemistry | 2 single points (`CalcSolvationDescriptors`); `1 + 6N` gradient single points (`CalcVibrationalThermochemistry`) |
| Reactivity | `CalcReactivityDescriptors`, `RunReactivity` | Vertical ΔSCF IP/EA, Fukui and local conceptual-DFT indices | 0 for `CalcReactivityDescriptors`; 3 single points for `RunReactivity` |

Every value has four pieces of meaning:

1. **Unit**: energy is Hartree or eV, charge/population is electrons, dipole
   is e·Bohr or Debye, bond order and ratios are dimensionless, and atom/MO
   indices are integer indices.
2. **Indexing**: atom and bond rows use the original RDKit order. Reordering a
   molecule requires remapping the result and its provenance.
3. **Invariance**: scalar energies and populations do not change under rigid
   translation/rotation; density descriptors also survive AO permutation and
   orbital phase changes when `S`, `C`, maps and p-component directions are transformed together.
4. **Missing semantics**: `None` means unavailable or undefined under the
   stated protocol. It never means zero. The surrounding report stores a
   reason such as `basis_not_retained`, `spin_density_unavailable_for_open_shell`,
   or `local_plane_nonplanar`.

The stable feature groups for the EAS study are:

- **Core features** (the default starting point for a model): xTB total energy,
  HOMO/LUMO/gap, dipole norm, atomic partial charge, and RDKit-bond Wiberg
  order. These are available from one ordinary single-point calculation.
- **Density features** (the candidate extension for the EAS study): Mulliken
  and Löwdin s/p/d populations and p/s ratios. They require the extended
  libxtb result with `keepOverlap=True`, but add no engine call.
- **Reactivity features**: vertical IP/EA, global hardness/softness and
  atom-level Fukui/local indices. They require three matched GFN1/2 charge
  states at the same geometry, so they are a separate analysis rather than
  part of the ordinary single-point feature vector.
- **Diagnostics**: conformer mean/std/min/max, GFN1-vs-GFN2 differences and
  electronic-temperature differences. These describe reliability and are not
  model inputs by default.

`DescriptorStatus()` exposes the maturity labels used by the project:
`recommended` is suitable for a first validation model (including core electronic
fields, Politzer–Murray ESP surface statistics, EPN, D4 polarizability/C6/C8,
energy decomposition, ALPB/GBSA solvation shifts, and harmonic thermochemistry),
`diagnostic` records calculation or protocol sensitivity (conformer/method/temperature
deltas, Fermi levels, and direct solvation polarization decomposition), and
`exploratory` is retained for testing new hypotheses (Mayer bond order, signed
frontier-bond overlap, π projections, local nucleophilicity proxy, orbital-weighted
Fukui, condensed ALIE/PDOS, COOP/ECOHP/bond polarity, Pipek–Mezey LMOs, and CAMM
atomic/molecular multipoles). `DescriptorCatalog()` organizes all definitions under
`molecular`, `atom`, `fragment`, `bond`, `density`, `reactivity`, `local_geometry`,
`esp_surface`, `localized_orbitals`, `properties`, `solvation`, `thermochemistry`,
and `diagnostic`. Their values are still fully specified and missing values are
explicit; the label only says whether the repository has established stable
predictive utility for them.

Frontier centroid separation remains a method-defined ground-state proxy and
should not be described as an excited-state observable. GFN0 shell populations
are density populations; GFN0 EEQ charges are a separate model quantity and
must not be used as a shell-population charge reconstruction.

For the full formulas and literature positioning, see
[density descriptors](density-descriptors.md) and the
[literature audit](literature-density-ml-audit.md).
