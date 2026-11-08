# Descriptor catalog

The public descriptor API has five layers. The complete machine-readable
catalog is returned by `eled.descriptor_catalog()`; this page states which
layers are appropriate for validation.

Python functions, report fields and analysis arguments use RDKit/Python's
`snake_case`. Unit suffixes are lower-case and explicit, for example
`xtb_gap_ev`, `xtb_total_energy_hartree`, `xtb_dipole_norm_e_bohr` and
`softness_per_ev`. There is one public spelling for each field.

| Layer | Main API | What it means | Engine calls during extraction |
| --- | --- | --- | --- |
| Molecular | `calc_xtb_descriptors` | Energy, frontier energies, dipole, charge and Wiberg summaries | 0 |
| Atom/bond topology | `calc_xtb_atom_descriptors`, `calc_xtb_bond_descriptors` | Engine charges, incident Wiberg sums and RDKit-bond Wiberg orders | 0 |
| Orbital density | `calc_orbital_populations`, `calc_shell_populations`, `calc_density_matrix` | Populations derived from the same `C`, `S`, AO and shell maps | 0; requires `keep_overlap=True` |
| Frontier/bond | `calc_frontier_descriptors`, `calc_frontier_bond_descriptors`, `calc_mayer_bond_descriptors` | Degeneracy-invariant HOMO/LUMO localization and signed bond decompositions | 0; closed-shell Mayer only |
| Reactivity | `calc_reactivity_descriptors`, `run_reactivity` | Vertical ΔSCF IP/EA, Fukui and local conceptual-DFT indices | 0 for `calc_*`; 3 single points for `run_*` |

Every value has four pieces of meaning:

1. **Unit**: energy is Hartree or eV, charge/population is electrons, dipole
   is e·Bohr or Debye, bond order and ratios are dimensionless, and atom/MO
   indices are integer indices.
2. **Indexing**: atom and bond rows use the original RDKit order. Reordering a
   molecule requires remapping the result and its provenance.
3. **Invariance**: scalar energies and populations do not change under rigid
   translation/rotation; density descriptors also survive AO permutation and
   orbital phase changes when `S`, `C` and maps are transformed together.
4. **Missing semantics**: `None` means unavailable or undefined under the
   stated protocol. It never means zero. The surrounding report stores a
   reason such as `basis_not_retained`, `spin_density_unavailable_for_open_shell`,
   or `local_plane_normal_undefined`.

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

`descriptor_status()` exposes the maturity labels used by the project:
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
