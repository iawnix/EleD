# Density descriptors

EleD extracts shell populations, bond indices and local π features from a
retained xTB result. Extraction adds no engine calls. `keepOverlap=True`
retains the overlap matrix `S`, coefficients `C`, occupations `n`, AO/shell
maps and p-component directions from the same calculation.

```python
from eled import XTBCalculator, CalcDensityDescriptors, CalcPiDescriptors

result = XTBCalculator(method="GFN2-xTB").Calculate(mol, keepOverlap=True)
density = CalcDensityDescriptors(mol, result)
pi = CalcPiDescriptors(mol, result, scheme="lowdin", planeTolAngstrom=0.1)
```

Atom and bond indices follow the original RDKit molecule. Functions taking
`mol` verify its topology, atom order and coordinates against the calculation.
`CalcDescriptorReport` includes these sections by default when overlap is retained.

Analyses of unchanged basis data share the overlap decomposition and density
matrices. Changes to input arrays or the atom count trigger revalidation and
recomputation. Cached arrays are owned by EleD. AO/shell maps must be
one-dimensional integer arrays.

## Density and shell populations

For spin-summed occupations, define

```text
P = C diag(n) Cᵀ
M = P S
D = S½ P S½
N = tr(P S)
```

`Cᵀ S C = I` is checked before analysis. `P` is an AO-basis matrix: its entries
change with the basis, while `tr(P S)` is invariant and measured in electrons.
`n_err_e = tr(P S) - n_e`, where `n_e` is the engine's modeled electron count
(or the occupation sum when that count is absent).

A shell's Mulliken population is the sum of its diagonal entries in `M`;
its Löwdin population is the corresponding sum in `D`. Atom rows collect
shells with the same angular momentum.

| Fields | Meaning | Unit |
| --- | --- | --- |
| `mul_s_e`, `mul_p_e`, `mul_d_e` | Mulliken s/p/d populations | e |
| `low_s_e`, `low_p_e`, `low_d_e` | Löwdin s/p/d populations | e |
| `mul_ps`, `low_ps` | p population / s population | dimensionless |
| `mul_qrec_e`, `low_qrec_e` | Reference valence minus the corresponding total population | e |
| `mul_n_e`, `low_n_e` (shell rows) | Population of one shell | e |
| `tr_ps_e`, `n_e`, `n_err_e` | Density electron count, reference count, difference | e |
| `shell_err_e` | Maximum shell-to-atom sum discrepancy for each scheme | e |

Absent s/p/d channels contribute zero; a zero s denominator gives `None` for
the ratio. Shells above d are outside this interface. Charge reconstruction
uses known reference valences for H, B, C, N, O, F, Si, P, S, Cl, Ge, As, Se,
Br and I; other elements return `None` with `reference_valence_unavailable`.
For GFN1/2 the Mulliken reconstruction can be checked against engine charges.
Löwdin charges are a different partition. GFN0's EEQ charges are a separate
model quantity, so both reconstructed charge fields are `None` for GFN0.

## Frontier-bond overlap descriptors

For each HOMO/LUMO group G, the frontier overlap descriptor uses

```text
PG = C_G C_Gᵀ / |G|
mul_ovlp_e(A,B) = 2 Σ[μ∈A,ν∈B] PGμν Sμν
```

This is a **one-electron-normalized group average**, in e, rather than the
occupied group's total electron population. Positive values indicate bonding
overlap; negative values indicate antibonding overlap. `frontier_group` and
`orb_idx` identify the group. Energy and occupation tolerances are shared with
`CalcFrontierDescriptors`; the group projector is invariant to orbital phases
and rotations within the selected subspace. Missing groups retain bond rows
with null overlap values and a `frontier_group_unavailable` reason.

## Local π projections

`CalcPiDescriptors` applies to sp²/aromatic atoms with exactly one complete
p shell and a defined plane. For a ring atom the fitted plane uses all atoms
in rings containing it; for other atoms it uses the atom and its neighbors.
At least three noncollinear points are required. The maximum distance from
the fitted plane must not exceed `planeTolAngstrom` (default 0.1 Å).

The result records `aoPDirections` with each AO row. These directions turn
the local unit normal into a vector `vA` supported on atom A's p shell.
Thus the projection remains well-defined after an AO permutation or sign change.
With `scheme="lowdin"` (default) or `"mulliken"`:

| Field | Löwdin | Mulliken | Unit |
| --- | --- | --- | --- |
| `pi_n_e` | `vAᵀ D vA` | `vAᵀ sym(M) vA` | e |
| `pi_q_e` | `nπ,ref − pi_n_e` | `nπ,ref − pi_n_e` | e |
| `pi_bo` | `(vAᵀ D vB)²` | `(vAᵀ M vB)(vBᵀ M vA)` | dimensionless |

The default π reference is **1 e for neutral, nonradical sp²/aromatic carbon**.
For other atoms `pi_n_e` can still be available, while `pi_q_e` requires an
explicit reference, for example `piRef={atom_index: 2.0}`. References must
lie between 0 and 2 e. Total atomic valence is not a π reference.
The returned `pi_ref` records explicit overrides with JSON-compatible string keys.
`CalcDensityDescriptors`, reports, batches and conformer analysis accept
`piScheme`, `piPlaneTolAngstrom` and `piRef`; see [report settings](descriptor-reports.md).

`pi_bo` is a projected bond index, not an absolute overlap population or a
complete decomposition of the total bond order. It requires a closed shell
and two valid endpoint projections. Reversing either normal leaves it unchanged.
The Löwdin expression is nonnegative; the Mulliken product can be signed.

Every atom and RDKit bond retains a row. Unavailable fields contain `None`
and a field-specific `missing` reason, such as `local_plane_degenerate`,
`local_plane_nonplanar`, `p_component_directions_unavailable`,
`pi_reference_unspecified` or `spin_density_unavailable_for_open_shell`.
The returned `normal` is coordinate metadata, not a scalar descriptor.

## Invariance and validation

Scalar density descriptors are invariant to rigid motion, MO phase changes
and consistent AO permutations/sign changes of `C`, `S`, maps and component
directions. `P` itself transforms with the AO basis. Rotating equally occupied
MOs preserves the density; rotating a complete selected frontier group preserves
its group-average overlap. Local normal signs do not affect π values.

Tests cover analytical nonorthogonal densities, electron and shell conservation,
undefined planes, atom/geometry provenance, degenerate rotations and actual
GFN0/1/2 calculations through ctypes and compiled interfaces. Benzene, pyridine,
ethylene and furan also exercise π invariance under rigid motion.

The density analyses follow Mulliken (1955), Löwdin (1950) and Mayer (1983).
Definitions and units are available programmatically from `DescriptorCatalog()`.
