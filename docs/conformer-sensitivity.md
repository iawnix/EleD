# Conformer sensitivity (1.3)

`RunConformerAnalysis` evaluates **explicitly supplied conformers of one RDKit
molecule** using native libxtb. It does not generate or optimize geometries.
Use [`OptimizeXTB`](geometry-preparation.md) before this
analysis when relaxed structures are required. The default requests the
extended library's overlap/basis data, as do descriptor reports.

```python
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import RunConformerAnalysis, XTBDiskCache

mol = Chem.AddHs(Chem.MolFromSmiles("OCCO"))
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.numThreads = 1
ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=8, params=params))
report = RunConformerAnalysis(
    mol, conf_ids=ids, record_id="ethylene-glycol",
    cache=XTBDiskCache("conformers.sqlite"),
)
print(report["statistics"]["xtb_gap_ev"])
with open("sensitivity.json", "w") as f:
    json.dump(report, f, indent=2, allow_nan=False)
```

Here the geometries are unoptimized ETKDG samples. A separate
[`examples/conformer_sensitivity.py`](../examples/conformer_sensitivity.py)
demonstrates explicit MMFF94 optimization, recording convergence, rejected
geometries, coordinates and the generation seed before selecting conformers.
MMFF is RDKit's force field, not an xTB geometry optimizer.

## Contract

- `conf_ids=None` uses all actual conformer IDs in RDKit order. IDs need not be
  contiguous; duplicate, missing and negative IDs are rejected. Select with
  `conf_ids`, leaving `options.conf_id=-1`. Molecule inputs are copied; all
  selected geometries are validated before calculation starts.
- All conformers use the same method, charge, spin and numerical settings.
  Results must match their input hashes; mixed engine/parameter/settings
  provenance raises. This high-level workflow uses the native Python adapter;
  individual `rdXTBTools.RunMol` calculations retain their existing API.
- `statistics[name]` contains `mean`, population `std` (`ddof=0`), `min`, `max`,
  `range`, `n_total`, `n_valid`, `n_missing`, `n_failed`. Units are in `units`.
  Missing values are not zero. No valid values yields null statistics. One
  conformer yields zero spread; this is not evidence of conformational stability.
- A calculation failure or error-level quality check excludes that conformer
  from all statistics and increments `n_failed`. Other absent fields increment
  `n_missing`. Their sum with `n_valid` always equals `n_total`. The full report
  and error/quality flags remain in `reports` in the requested order.
- `status` describes selected-column coverage: `complete`, `partial`, or
  `failed` (all conformers failed). `complete` does not override quality warnings
  such as fractional occupations. Inspect the individual reports.
- Invalid inputs/settings, unavailable engine capabilities, environment errors,
  corrupted storage and unexpected analysis errors **raise**. With
  `return_errors=False`, single-point calculation failures also raise.
- `XTBCache` and `XTBDiskCache` reuse successful individual single points with
  the existing geometry/settings/engine/implementation keys. Analysis choices
  do not change the underlying calculation; reports are recomputed from results.

## Which descriptors are summarized?

`ConformerDescriptorColumns(charge=0, frontier_depth=2)` returns the default
ordered scalar schema. Energies, frontier gaps, charge/WBO summaries, frontier
participation/entropy/overlap and centroid **separation** are included. Counts
are retained to expose changes in frontier classification. Indices, categorical
metadata, numerical settings, Cartesian dipole components and absolute orbital
centroids are excluded. Charged molecules also exclude dipole norms, since
those depend on the coordinate origin. No alignment is needed for these scalar
aggregates; rotation/translation invariance is tested on neutral and charged
molecules. `columns=[...]` may select a nonempty subset of this schema.

Atom and bond scalar fields also receive indexed statistics (see below).
Fragment and individual orbital sections remain per-conformer data; this API
does not equate a specific MO index across geometries.
Near-degenerate frontier groups follow existing report tolerances, and can
change composition with conformation; inspect counts and quality flags.

## Interpretation

These are equally weighted statistics over the **supplied geometries**, including
any duplicates. This is the default and remains the main sensitivity report.
For an optional ensemble-style summary, pass `boltzmann_weighted=True`:

```python
report = RunConformerAnalysis(
    mol, conf_ids=ids, boltzmann_weighted=True,
    boltzmann_temperature=300.0,
)
weighted = report["boltzmann"]["statistics"]["xtb_gap_ev"]
weights = report["boltzmann"]["weights"]
```

The optional weights use xTB electronic energies,
`exp(-(E_i-E_min)/(k_B T))`. They describe an energy-weighted descriptor
aggregate; they are not a free-energy conformer population because vibrational
entropy, zero-point energy and solvent effects are not included. The temperature
is a conformer-aggregation parameter and is separate from xTB's electronic
occupation temperature. Both the equal-weight statistics and the weighted
statistics are retained in the report.

The recorded example retains 11 MMFF-converged conformers across five molecules.
The GFN2 HOMO–LUMO gap spans about 0.250 eV for ethylene glycol, 0.333 eV for
butane and 0.151 eV for phenethylamine. Ethanol and benzene each retain only one
geometry; their zero spread is uninformative about unsampled conformations.
These orbital gaps are not vertical optical excitation energies.

π extraction accepts `pi_scheme`, `pi_plane_tol_a` and `pi_ref` as described
in [report settings](descriptor-reports.md). The settings apply to every supplied
conformer and are retained in both the ensemble and individual reports.
Equal-weight and electronic-energy-weighted π statistics use that definition.

## Atom and bond statistics

`atom_statistics[atom_index][field]` and `bond_statistics[bond_index][field]`
report `mean`, `std`, `min`, `max` and coverage in original RDKit order. JSON
index keys are strings. `atom_units` and `bond_units` map field names to units.
Each conformer contributes at most one value per indexed field:

```python
charge = report["atom_statistics"]["0"]["atoms.q_e"]
homo = report["bond_statistics"]["0"]["frontier_bonds.homo.mul_ovlp_e"]
lumo = report["bond_statistics"]["0"]["frontier_bonds.lumo.mul_ovlp_e"]
pi_bond = report["bond_statistics"]["0"]["pi.pi_bo"]
```

Unavailable values retain their slots: `n_total = n_valid + n_missing + n_failed`.
An entirely missing field has null statistics, and one valid observation has
zero population standard deviation (`ddof=0`). Atom π populations never enter
bond statistics; HOMO and LUMO overlap populations have separate fields.
`frontier_group_composition_changed` records changes in the selected orbital
indices. It does not assert that orbitals with the same indices have been tracked
between conformers; the per-conformer reports retain each group's definition.

With Boltzmann weighting, `boltzmann.atom_statistics` and
`boltzmann.bond_statistics` contain weighted means and standard deviations.
`weight_sum` records the available probability mass before renormalization for
each field. A finite observation with an underflowed zero weight still counts
as valid; if all available values have zero mass, the weighted result is null
with `no_positive_weight_for_available_values`. Weighted and equal-weight local
statistics use the same unit maps.

`DescriptorSensitivity` compares molecular and local statistics across matched
method/temperature reports. It checks identical ordered geometries and protocol,
then requires matching valid-conformer subsets for each compared field.
See [density descriptors](density-descriptors.md#sensitivity-diagnostics).
