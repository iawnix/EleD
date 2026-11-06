# Conformer sensitivity (0.9)

`run_conformer_analysis` evaluates **explicitly supplied conformers of one RDKit
molecule** using native libxtb. It does not generate or optimize geometries.
Use [`optimize_xtb`](../README.md#xtb-geometry-optimization) before this
analysis when relaxed structures are required. The default requests the
extended library's overlap/basis data, as do descriptor reports.

```python
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_conformer_analysis, XTBDiskCache

mol = Chem.AddHs(Chem.MolFromSmiles("OCCO"))
params = AllChem.ETKDGv3()
params.randomSeed = 42
params.numThreads = 1
ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=8, params=params))
report = run_conformer_analysis(
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

`conformer_descriptor_columns(charge=0, frontier_depth=2)` returns the default
ordered scalar schema. Energies, frontier gaps, charge/WBO summaries, frontier
participation/entropy/overlap and centroid **separation** are included. Counts
are retained to expose changes in frontier classification. Indices, categorical
metadata, numerical settings, Cartesian dipole components and absolute orbital
centroids are excluded. Charged molecules also exclude dipole norms, since
those depend on the coordinate origin. No alignment is needed for these scalar
aggregates; rotation/translation invariance is tested on neutral and charged
molecules. `columns=[...]` may select a nonempty subset of this schema.

Atom/fragment/bond and individual orbital sections remain per-conformer data;
this API does not average them or equate a specific MO index across geometries.
Near-degenerate frontier groups follow existing report tolerances, and can
change composition with conformation; inspect counts and quality flags.

## Interpretation

These are equally weighted statistics over the **supplied geometries**, including
any duplicates. They are not Boltzmann populations, free energies or predictions
at a conformational temperature. The xTB electronic temperature describes
orbital occupations, not conformer weights. No geometries are deduplicated or
chosen by their calculated energy inside this API.

The recorded example retains 11 MMFF-converged conformers across five molecules.
The GFN2 HOMO–LUMO gap spans about 0.250 eV for ethylene glycol, 0.333 eV for
butane and 0.151 eV for phenethylamine. Ethanol and benzene each retain only one
geometry; their zero spread is uninformative about unsampled conformations.
These orbital gaps are not vertical optical excitation energies.
