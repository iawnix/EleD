# Electronic descriptors and population analysis (0.6)

For standardized reports, quality flags, CSV/JSONL export and vertical charge-state
analysis, see [descriptor reports](descriptor-reports.md).

All extraction functions below use an existing result and perform no additional
quantum calculation. Units are encoded in field names. Results are JSON-compatible
Python dictionaries/lists, with `None` for unavailable scalar values rather than
fabricated zeros. All calculation paths use libxtb; no tblite or ASE Python dependency is required.

## Interface overview

| Function | Output | Supported input |
| --- | --- | --- |
| `calc_xtb_descriptors(result, frontier_depth=2, occupation_threshold=1e-7)` | Flat molecular features | Native `XTBResult` or compiled `rdXTBTools.XTBResults` |
| `calc_xtb_atom_descriptors(result)` | One row per atom | Both result paths |
| `calc_xtb_bond_descriptors(mol, result)` | One row per existing RDKit bond | Either result type with matching input provenance |
| `calc_xtb_fragment_descriptors(result, fragments)` | Named fragment records, fragment-pair records, unassigned atoms | Both result paths; explicit original atom indices |
| `calc_orbital_populations(result, scheme="lowdin", orbital_indices=None)` | Atom × selected-MO populations and electron populations | Same-calculation coefficients, overlap and AO map; verified native basis extension |
| `calc_frontier_descriptors(mol, result, fragments=None)` | HOMO/LUMO subspace populations, localization and atom-centroid proxies | Same as population analysis, plus matching molecule/conformer |

The Python `rdXTBTools.RunMol` wrapper records input provenance from its owned
molecule snapshot. Both Python result types support bond mapping, annotation and
frontier spatial descriptors. Standalone C++ results do not carry the Python
provenance record; molecule-bound helpers reject missing provenance.

## Molecular, atomic and fragment descriptors

```python
from rdkit_eth import (
    run_xtb, calc_xtb_descriptors, calc_xtb_atom_descriptors,
    calc_xtb_bond_descriptors, calc_xtb_fragment_descriptors, annotate_mol,
)

# mol must already have explicit hydrogen atoms and a 3D conformer.
result = run_xtb(mol, method="GFN2-xTB")
scalars = calc_xtb_descriptors(result, frontier_depth=2)
atoms = calc_xtb_atom_descriptors(result)
bonds = calc_xtb_bond_descriptors(mol, result)
fragments = calc_xtb_fragment_descriptors(result, {"first": [0], "second": [1, 2]})
annotated = annotate_mol(mol, result)
```

`calc_xtb_descriptors` retains the original method/backend/version, total energy
(Hartree), HOMO/LUMO/gap (eV) and electron/orbital counts. It adds:

- `xtb_homo_index`, `xtb_lumo_index`: zero-based indices; `None` if absent.
- `xtb_homo_minus_1_eV`, `xtb_lumo_plus_1_eV`, etc., up to `frontier_depth`.
  Values beyond the occupied/virtual lists are `None`; depth zero omits them.
- `xtb_dipole_x_au`, `_y_au`, `_z_au`, `xtb_dipole_norm_au` and
  `xtb_dipole_norm_debye`. One e·Bohr = 2.541746473 Debye.
- `xtb_charge_sum_e`, `_min_e`, `_max_e`, `_mean_e`, `_std_e`, `_abs_sum_e`.
  Standard deviation uses all atoms and `ddof=0`, including explicit hydrogens.
- `xtb_wiberg_pair_sum`: sum over all `i < j`, including pairs outside RDKit bonds.
- `xtb_occupation_threshold`: records the frontier convention used.

An orbital is occupied when its occupation exceeds the specified threshold.
HOMO is the last occupied MO and LUMO the first unoccupied MO in energy order.
Fractional occupations can change that classification. The gap is not an optical
excitation energy; the input unpaired-electron count is not a predicted spin state.
For charged molecules, the dipole and its norm depend on the supplied coordinate
origin. Total energies have method-specific references.

Atom rows contain `atom_index`, `partial_charge_e` and `wiberg_sum`. The latter
is the sum of all off-diagonal WBOs involving that atom, not formal valence.
Bond rows contain `bond_index`, `begin_atom_index`, `end_atom_index` and
`wiberg_order`. They do not create bonds or change RDKit bond orders.

Fragment groups are nonempty **disjoint** collections of unique integer indices
into the calculation's original atom order. Indices include hydrogen atoms;
hydrogens are not added to a fragment implicitly. A complete partition is optional.
Duplicate, negative, Boolean, out-of-range or overlapping indices are rejected.
Each named fragment reports:

- `charge_e`, `charge_min_e`, `charge_max_e`, `charge_std_e`, `num_atoms`;
- `internal_wiberg_sum`: each pair inside the fragment counted once;
- `external_wiberg_sum`: all pairs to atoms outside, including unassigned atoms.

`pairs` reports each fragment pair's WBO sum once. `unassigned_atom_indices`
identifies omitted atoms. For a complete partition, fragment charges sum to the
molecular charge, and internal plus interfragment WBO sums recover the total
pair sum. No spatial criterion or bond-order cutoff is silently applied.

`annotate_mol` now writes `xtb_wiberg_sum` on atoms, in addition to partial charges,
bond WBOs and molecular scalar descriptors. `atom_wiberg_prop` can override its
name. Atomic charge and WBO property names must differ. All validation precedes
in-place writes; geometry/topology changes are rejected.

## Obtain same-calculation basis data

```python
result = run_xtb(mol, method="GFN2-xTB", keep_overlap=True)
print(result.ao_atom_indices)  # one zero-based atom index per AO
print(result.overlap_matrix)  # AO × AO, same basis/order as coefficients
```

Select the [extended libxtb installation](xtb-patched-build.md) before running
this example. `keep_overlap` is also supported in `XTBOptions` and
`run_xtb_batch`, and is part of cache identity. It defaults to false, accepts
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
from rdkit_eth import calc_orbital_populations

lowdin = calc_orbital_populations(result)
mulliken = calc_orbital_populations(result, scheme="mulliken")
selected = calc_orbital_populations(result, orbital_indices=[3, 4])
```

For real AO coefficients `C`, overlap `S`, and atom A:

- Löwdin amplitudes are `S^(1/2) C`; sum their squares over A's AOs.
- Mulliken populations sum `C_mu,i * (S C)_mu,i` over A's AOs.

The analysis requires finite square matrices, a symmetric positive-definite
well-conditioned overlap, a valid AO map covering every atom, and
`max(abs(C.T S C - I)) <= 1e-6` by default. The measured orthonormality error
is reported. After this check, each MO column is divided by its measured
norm to remove roundoff; large norm defects are rejected rather than hidden.

`atom_populations` has axes `[atom, selected orbital]`, and each column sums
to one. `orbital_indices` records its column identities. `electron_populations_e`
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
from rdkit_eth import calc_frontier_descriptors

frontier = calc_frontier_descriptors(
    mol, result,
    fragments={"donor": donor_atom_indices, "acceptor": acceptor_atom_indices},
    degeneracy_ev=1e-3,
)
print(frontier["homo"]["fragment_populations"])
print(frontier["lumo"]["participation_ratio"])
print(frontier["homo_lumo_atom_centroid_separation_angstrom"])
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

- `homo_lumo_atom_overlap = sum_A sqrt(p_H,A * p_L,A)`; a dimensionless
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
