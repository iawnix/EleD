# Local steric and ring descriptors

EleD analyzes selected atoms of an explicit-hydrogen RDKit molecule with a 3D
conformer. All indices follow the original molecule. Every function accepts
`conf_id=-1`, which selects RDKit's default conformer, and records its actual ID.
These functions require RDKit and NumPy and analyze existing coordinates.

```python
from eled import PrepareConformers, CalcHOMA, CalcRingPuckering

prepared = PrepareConformers('c1ccc2ccccc2c1')
if not prepared.success:
    raise RuntimeError(prepared.report)
mol = prepared.molecule
rings = mol.GetRingInfo().AtomRings()
reports = [CalcHOMA(mol, ring) for ring in rings]
print([r['values']['homa'] for r in reports])
```

Separate calls can analyze overlapping fused rings. A ring selection is a closed
bonded path of distinct atoms, listed consecutively, with the last bonded to the
first. The first atom appears once. The caller chooses which cycles to analyze.

## Outputs and provenance

Outputs use `schema_version=1`, a descriptor-specific `report_type` and
`protocol_version='1.0'`. They contain `values`, `units`, `missing`, `selection`,
`analysis_settings`, and `provenance`. Provenance includes the conformer ID,
indexed structure and coordinate hashes, atomic numbers and coordinates.
Undefined quantities are `None` with a reason. Invalid input raises `ValueError`.
Reports can be serialized with `json.dumps(report, allow_nan=False)`.
`ExportDescriptorReports` handles the electronic report schema; save local
geometry reports as JSON or JSONL directly.

| Function | Selection | Scalar results |
| --- | --- | --- |
| `CalcSterimol(mol, atom_indices, axis=(base, attached))` | Substituent atoms and directed bonded axis | `l_a`, `b1_a`, `b5_a`, Angstrom |
| `CalcBuriedVolume(mol, atom_indices, center=i)` | Occupying atoms and excluded center atom | `vbur_percent`, percent |
| `CalcPlaneDescriptors(mol, atom_indices, other_atom_indices=None)` | One or two atom sets, at least three atoms each | `plane_rms_a`, `plane_max_a`, optional `centroid_distance_a`, Angstrom; `plane_angle_deg`, degree |
| `CalcHOMA(mol, ring_atom_indices)` | Ordered ring | `homa`, dimensionless; individual `bonds` with `length_a` |
| `CalcRingPuckering(mol, ring_atom_indices)` | Ordered ring of at least four atoms | `q_total_a`, optional `q_alternating_a`, Angstrom; six-membered `theta_deg`; per-mode `q_a` and `phi_deg` |

Selections explicitly include the hydrogen atoms to be counted. Steric functions
use RDKit van der Waals radii with `radius_scale=1.0`. Pass `radii_angstrom` as a
full-length vector in original atom order for another radius set. Effective
radii and the RDKit version are saved. Comparisons require the same radius set,
hydrogen selection, geometric preparation and descriptor parameters.

## Sterimol

The axis runs from the base to the attached atom. The base must be outside the
selection and the attached atom inside it. Positions are measured relative to
the base. Each selected atom is a sphere with its specified radius.

- `L = max(axial_position + radius) + l_correction_angstrom`, default correction
  0.4 Angstrom.
- `B1 = min_direction max_atom(transverse_position · direction + radius)`.
- `B5 = max_atom(distance_to_axis + radius)`.

B1 samples a full transverse circle with `angular_samples=7200`. Its upward
sampling error is bounded by `max_transverse_distance * pi / angular_samples`,
recorded as `b1_error_bound_angstrom`. Increasing the sample count tightens this
bound. A single on-axis sphere gives B1=B5=radius. The geometry and chosen radii
fully specify this sphere-envelope implementation. Radius sets and L conventions
must match when comparing published Sterimol values.

## Buried volume

`CalcBuriedVolume` integrates the union of selected atom spheres within a sphere
of `sphere_radius_angstrom=3.5` centered on the requested atom. The center atom
must be excluded explicitly. Hydrogens participate if selected.

A centered Cartesian midpoint grid uses requested `spacing_angstrom=0.1`.
The diameter is divided into an integer number of cells, and the effective
spacing and occupied/total point counts are recorded. The percentage is
`100 * occupied_grid_points / grid_points`. Overlaps count once. Rotation
invariance is approximate at finite resolution. Repeat with smaller spacing to
assess numerical sensitivity. The implementation supports 4–256 cells per
sphere diameter and processes one grid plane at a time. Radius scaling and
hydrogen treatment are explicit, including for comparisons to ligand %Vbur
literature.

## Planes and ring puckering

Plane descriptors use unweighted least squares. `plane_rms_a` and `plane_max_a`
measure absolute deviations from that plane. The angle between two unoriented
planes lies in [0,90] degrees. Collinear selections and a degenerate smallest
eigenvalue give `plane_normal_not_unique`.

Ring puckering follows Cremer–Pople (JACS 1975, 97, 1354–1358). For centered
coordinates r_j, j=0..N−1, the normal is parallel to
`sum(r_j sin(2*pi*j/N)) × sum(r_j cos(2*pi*j/N))`. Projected heights are z_j.
For modes m=2..floor((N−1)/2), cosine and sine components are
`sqrt(2/N)*sum(z_j*cos(2*pi*m*j/N))` and
`-sqrt(2/N)*sum(z_j*sin(2*pi*m*j/N))`.
Even rings additionally have the signed alternating mode
`sum((-1)^j*z_j)/sqrt(N)`. Q is `sqrt(sum(z_j**2))`.
Six-membered rings report `theta=atan2(q2,q3)` in [0,180] degrees.
Mode phases use [0,360) degrees and depend on the ring start and direction.
Amplitudes below 1e-10 Angstrom have undefined phases. Planar six-membered rings
have undefined theta. Degenerate normals give explicit missing values.

## HOMA

`HOMA = 1 - mean(alpha_AB * (R_AB - Ropt_AB)**2)`.
The built-in C–C parameters are Ropt=1.388 Angstrom and alpha=257.7 Angstrom^-2
from Krygowski, J. Chem. Inf. Comput. Sci. 1993, 33, 70–78.
Negative results are retained. Every ring edge contributes, including the
closing edge. Aromaticity labels are not needed to evaluate the geometric index.

For heteroatom rings, supply a complete parameter table such as
`parameters={(6, 6): (Ropt_CC, alpha_CC), (6, 7): (Ropt_CN, alpha_CN)}` together
with a nonempty `parameter_source`. Custom tables replace the default table.
An unsupported pair makes the whole HOMA value missing with reason
`unsupported_element_pair`; all measured bond lengths remain in the report.
HOMA depends on the geometry preparation and parameter set. It provides a
geometric measure of aromaticity, to be interpreted with the chemical system.

## Multiple conformers

Call the same function and selections with each actual conformer ID. Retain
those individual reports. For example, compute the equal-weight mean and
population standard deviation of Sterimol B5 with:

```python
import numpy as np
from eled import CalcSterimol

rows = [CalcSterimol(mol, substituent, axis=axis, conf_id=c.GetId())
        for c in mol.GetConformers()]
b5 = np.array([row['values']['b5_a'] for row in rows])
mean, std = float(b5.mean()), float(b5.std())
# For supplied normalized conformer weights w:
# mean = float(np.average(b5, weights=w))
# std = float(np.sqrt(np.average((b5-mean)**2, weights=w)))
```

Keep missing values and coverage when aggregating other descriptors. Circular
phase angles need circular statistics. Boltzmann weights require appropriate
conformer energies and a specified temperature.

See [the runnable example](../examples/local_geometry.py) for anisole Sterimol,
local volume, the two rings of naphthalene and cyclohexane puckering. Analytical
tests cover a single-sphere Sterimol reference, intersecting-sphere volumes,
ideal HOMA lengths, Fourier ring modes, and rigid-motion/index invariance.
