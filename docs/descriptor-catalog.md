# Descriptor catalog

The public descriptor API has five layers. The complete machine-readable
catalog is returned by `eled.DescriptorCatalog()`; this page states which
layers are appropriate for validation.

Public descriptor functions follow RDKit's `Calc...` / `Run...` naming.
Report fields and Python analysis arguments use `snake_case`. Unit suffixes are lower-case and explicit, for example
`xtb_gap_ev`, `xtb_energy_ha`, `xtb_dip_ebohr` and
`s_1ev`. There is one public spelling for each field.

The naming rules keep feature columns short while retaining their meaning:

- `ha`, `ev`, `e`, `ebohr`, `debye` and `a` identify Hartree, eV, electron,
  e·Bohr, Debye and Ångström values;
- `q`, `wbo`, `mbo`, `mul`, `low`, `eta`, `s` and `omega` identify charge,
  Wiberg, Mayer, Mulliken, Löwdin, hardness, softness and electrophilicity;
- `idx`, `n`, `sd`, `ipr`, `pr` and `ent` mean index, count, standard deviation,
  inverse participation ratio, participation ratio and entropy;
- `xtb_` marks molecular fields directly reported by xTB. Atom and bond rows
  use their shorter local names (`q_e`, `wbo`, `mul_p_e`, `low_ps`).

| Layer | Main API | What it means | Engine calls during extraction |
| --- | --- | --- | --- |
| Molecular | `CalcXTBDescriptors` | Energy, frontier energies, dipole, charge and Wiberg summaries | 0 |
| Atom/bond topology | `CalcXTBAtomDescriptors`, `CalcXTBBondDescriptors` | Engine charges, incident Wiberg sums and RDKit-bond Wiberg orders | 0 |
| Orbital density | `CalcOrbitalPopulations`, `CalcShellPopulations`, `CalcDensityMatrix` | Populations derived from the same `C`, `S`, AO and shell maps | 0; requires `keep_overlap=True` |
| Frontier/bond | `CalcFrontierDescriptors`, `CalcFrontierBondDescriptors`, `CalcMayerBondDescriptors` | Degeneracy-invariant HOMO/LUMO localization and signed bond decompositions | 0; closed-shell Mayer only |
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
  libxtb result with `keep_overlap=True`, but add no engine call.
- **Reactivity features**: vertical IP/EA, global hardness/softness and
  atom-level Fukui/local indices. They require three matched GFN1/2 charge
  states at the same geometry, so they are a separate analysis rather than
  part of the ordinary single-point feature vector.
- **Diagnostics**: conformer mean/std/min/max, GFN1-vs-GFN2 differences and
  electronic-temperature differences. These describe reliability and are not
  model inputs by default.

`DescriptorStatus()` exposes the maturity labels used by the project:
`recommended` is suitable for a first validation model, `diagnostic` records
calculation or protocol sensitivity, and `exploratory` is retained for testing
new hypotheses. Mayer bond order, signed frontier-bond overlap, π projections,
and the local nucleophilicity proxy are exploratory. Their values are still
fully specified and missing values are explicit; the label only says that the
repository has not yet shown stable predictive utility for them.

Frontier centroid separation remains a method-defined ground-state proxy and
should not be described as an excited-state observable. GFN0 shell populations
are density populations; GFN0 EEQ charges are a separate model quantity and
must not be used as a shell-population charge reconstruction.

For the full formulas and literature positioning, see
[density descriptors](density-descriptors.md) and the
[literature audit](literature-density-ml-audit.md).
