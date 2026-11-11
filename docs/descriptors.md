# Electronic descriptors and density analysis (1.3)

For standardized reports, quality flags, CSV/JSONL export and vertical charge-state
analysis, see [descriptor reports](descriptor-reports.md).

All extraction functions below use an existing result and perform no additional
quantum calculation. Units are encoded in field names. Results are JSON-compatible
Python dictionaries/lists, with `None` for unavailable scalar values rather than
fabricated zeros. All calculation paths use libxtb; no tblite or ASE Python dependency is required.

## Interface overview

| Function | Output | Supported input |
| --- | --- | --- |
| `CalcXTBDescriptors(result, frontier_depth=2, occupation_threshold=1e-7)` | Flat molecular features | Native `XTBResult` or compiled `rdXTBTools.XTBResults` |
| `CalcXTBAtomDescriptors(result)` | One row per atom | Both result paths |
| `CalcXTBBondDescriptors(mol, result)` | One row per existing RDKit bond | Either result type with matching input provenance |
| `CalcXTBFragmentDescriptors(result, fragments)` | Named fragment records, fragment-pair records, unassigned atoms | Both result paths; explicit original atom indices |
| `CalcOrbitalPopulations(result, scheme="lowdin", orb_idx=None)` | Atom × selected-MO populations and electron populations | Same-calculation coefficients, overlap and AO map; verified native basis extension |
| `CalcFrontierDescriptors(mol, result, fragments=None)` | HOMO/LUMO subspace populations, localization and atom-centroid proxies | Same as population analysis, plus matching molecule/conformer |
| `CalcDensityDescriptors(mol, result)` | Density conservation, shell populations, Mayer/Wiberg, frontier-bond and π descriptors | Same-calculation C/S/maps; `keep_overlap=True` |

The Python `rdXTBTools.RunMol` wrapper records input provenance from its owned
molecule snapshot. Both Python result types support bond mapping, annotation and
frontier spatial descriptors. Standalone C++ results do not carry the Python
provenance record; molecule-bound helpers reject missing provenance.

## Molecular, atomic and fragment descriptors

```python
from eled import (
    XTBCalculator, CalcXTBDescriptors, CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors, CalcXTBFragmentDescriptors, AnnotateMol,
)

# mol must already have explicit hydrogen atoms and a 3D conformer.
result = XTBCalculator(method="GFN2-xTB").calculate(mol)
scalars = CalcXTBDescriptors(result, frontier_depth=2)
atoms = CalcXTBAtomDescriptors(result)
bonds = CalcXTBBondDescriptors(mol, result)
fragments = CalcXTBFragmentDescriptors(result, {"first": [0], "second": [1, 2]})
annotated = AnnotateMol(mol, result)
```

`CalcXTBDescriptors` retains the original method/backend/version, total energy
(Hartree), HOMO/LUMO/gap (eV) and electron/orbital counts. It adds:

- `xtb_homo_idx`, `xtb_lumo_idx`: zero-based indices; `None` if absent.
- `xtb_homo_m1_ev`, `xtb_lumo_p1_ev`, etc., up to `frontier_depth`.
  Values beyond the occupied/virtual lists are `None`; depth zero omits them.
- `xtb_dip_x_ebohr`, `xtb_dip_y_ebohr`, `xtb_dip_z_ebohr`, `xtb_dip_ebohr` and
  `xtb_dip_debye`. One e·Bohr = 2.541746473 Debye.
- `xtb_q_sum_e`, `_min_e`, `_max_e`, `_mean_e`, `_sd_e`, `_abs_e`.
  Standard deviation uses all atoms and `ddof=0`, including explicit hydrogens.
- `xtb_wbo_sum`: sum over all `i < j`, including pairs outside RDKit bonds.
- `xtb_occ_cut`: records the frontier convention used.

An orbital is occupied when its occupation exceeds the specified threshold.
HOMO is the last occupied MO and LUMO the first unoccupied MO in energy order.
Fractional occupations can change that classification. The gap is not an optical
excitation energy; the input unpaired-electron count is not a predicted spin state.
For charged molecules, the dipole and its norm depend on the supplied coordinate
origin. Total energies have method-specific references.

Atom rows contain `atom_index`, `q_e` and `wbo_sum`. The latter
is the sum of all off-diagonal WBOs involving that atom, not formal valence.
Bond rows contain `bond_index`, `begin_atom_index`, `end_atom_index` and
`wbo`. They do not create bonds or change RDKit bond orders.

Fragment groups are nonempty **disjoint** collections of unique integer indices
into the calculation's original atom order. Indices include hydrogen atoms;
hydrogens are not added to a fragment implicitly. A complete partition is optional.
Duplicate, negative, Boolean, out-of-range or overlapping indices are rejected.
Each named fragment reports:

- `q_e`, `q_min_e`, `q_max_e`, `q_sd_e`, `num_atoms`;
- `wbo_intra`: each pair inside the fragment counted once;
- `wbo_inter`: all pairs to atoms outside, including unassigned atoms.

`pairs` reports each fragment pair's WBO sum once. `unassigned_atom_indices`
identifies omitted atoms. For a complete partition, fragment charges sum to the
molecular charge, and internal plus interfragment WBO sums recover the total
pair sum. No spatial criterion or bond-order cutoff is silently applied.

`AnnotateMol` now writes `xtb_wbo_sum` on atoms, in addition to partial charges,
bond WBOs and molecular scalar descriptors. `atom_wiberg_prop` can override its
name. Atomic charge and WBO property names must differ. All validation precedes
in-place writes; geometry/topology changes are rejected.

The optional YAeHMOP adapter follows the same compact convention. Its scalar
fields are `eht_homo_ev`, `eht_lumo_ev`, `eht_gap_ev`, `eht_fermi_ev`,
`eht_energy_ev`, `eht_n_orb` and `eht_n_e`; the `eht_` prefix distinguishes
the RDKit EHT calculation from xTB fields.

## Obtain same-calculation basis data

```python
result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
print(result.ao_atom_indices)  # one zero-based atom index per AO
print(result.overlap_matrix)  # AO × AO, same basis/order as coefficients
```

Select the [extended libxtb installation](xtb-patched-build.md) before running
this example. `keep_overlap` is also supported in `XTBOptions` and
`RunXTBBatch`, and is part of cache identity. It defaults to false, accepts
only a Boolean and implies retaining coefficients. Extra storage costs O(NAO²).

Compiled callers use `rdXTBTools.RunMol(mol, keepOverlap=True)` and retrieve
`GetOverlapMatrix()`, `GetAOAtomMap()`, `GetAOShellMap()`, `GetShellAtomMap()`
and `GetShellAngularMomenta()`. Native result fields are `overlap_matrix`,
`ao_atom_indices`, `ao_shell_indices`, `shell_atom_indices` and
`shell_angular_momenta`. Maps use zero-based indices in native spherical AO
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
selected = CalcOrbitalPopulations(result, orb_idx=[3, 4])
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
    degeneracy_ev=1e-3,
)
print(frontier["homo"]["frag_pop"])
print(frontier["lumo"]["pr"])
print(frontier["hl_dist_a"])
```

HOMO/LUMO groups collect contiguous orbitals in the same occupied/virtual class,
within `degeneracy_ev` of the frontier energy and `occupation_tolerance=1e-6`
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

## Shell, Mayer, frontier-bond and π descriptors

```python
from eled import CalcDensityDescriptors

result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
density = CalcDensityDescriptors(mol, result)
```

`CalcDensityMatrix` exposes `P`, `tr_ps_e` and the electron-conservation
error. `CalcShellPopulations` reports Mulliken and Löwdin s/p/d electron
counts, p/s ratios and shell rows. `CalcMayerBondDescriptors` places Mayer
and Wiberg values beside each RDKit bond; open-shell Mayer entries are `None`
with an explicit spin-density reason. `CalcFrontierBondDescriptors` uses
near-degenerate HOMO/LUMO projectors and signed Mulliken overlap populations.
`CalcPiDescriptors` reports p-shell projection for aromatic/sp2 atoms and
π-bond order where a local plane exists. Units, invariance and missing-value
semantics are included in each section. See [density-descriptors](density-descriptors.md)
for formulas, units and applicability. The scope and literature
positioning audit is in [literature-density-ml-audit](literature-density-ml-audit.md).
