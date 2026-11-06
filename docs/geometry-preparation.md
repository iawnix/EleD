# Checked SMILES-to-3D preparation (0.10)

`prepare_conformers` builds an explicit-H RDKit molecule plus an auditable
preparation report. It uses RDKit embedding and MMFF94/UFF by default. With
`optimize="xTB"`, it can run the separate gradient optimizer from
`rdkit_eth.optimize_xtb` after embedding. This is a workflow around RDKit,
not a new embedding algorithm.

```python
from rdkit_eth import prepare_conformers, run_xtb

prepared = prepare_conformers(
    'C[C@H](O)C(=O)O', num_conformers=8, random_seed=42,
    optimize='MMFF94',
)
if prepared.success:
    molecule = prepared.molecule
    result = run_xtb(molecule, conf_id=prepared.report['recommended_conformer_id'])
else:
    print(prepared.report['error'], prepared.report['attempts'])
```

The molecule can also be passed to `run_conformer_analysis`, including all its
accepted conformers. Requested and accepted counts can differ. The
[`prepare_geometry.py` example](../examples/prepare_geometry.py) writes the
report and an SDF, with an optional `--xtb` calculation.

## Retry and acceptance policy

Defaults request one conformer, ETKDGv3, two deterministic seeds per strategy,
and up to 300 internal embedding iterations per requested conformer. Try standard
distance-geometry initialization, then random coordinates. Repeat these strategies
with ETKDGv4 and srETKDGv3 if the installed RDKit exposes them. Missing default
presets are recorded; explicitly requesting an unavailable or unknown preset
raises. `fallback_methods=()` disables other presets. `random_seed=42` is fixed
by default; retries use `(seed + 104729 * attempt_index) % 2147483647`.

Every attempt enforces specified stereochemistry and uses one thread. Within
an embedding call, heavy-atom RMS pruning is 0.2 Å. Each candidate must have
finite 3D coordinates and no atom pair closer than 0.3 Å. Specified tetrahedral
atom tags and E/Z bonds are checked **from the actual coordinates**, not just
from the tags on the input molecule. Those tags are removed on a temporary
copy before 3D stereochemistry is reassigned and compared in the same atom order.

Default MMFF94 optimization must converge within 500 iterations. `optimize='UFF'`
explicitly chooses UFF; missing force-field parameters produce a failed report
before embedding, with no force-field substitution. `optimize=None` returns
checked, unoptimized embeddings and labels them accordingly. After optimization,
coordinate and stereo checks run again. Nonconverged or stereochemically changed
coordinates never enter the returned molecule. `optimize='xTB'` uses GFN2-xTB by
default (override with `xtb_method`), records the full gradient optimization
report under each candidate, and requires the reviewed extended libxtb. The
standalone `optimize_xtb()` API exposes convergence tolerances, bounded step
size, fixed atoms and a separate `scf_max_iterations` limit. `prepare_conformers`
uses `optimization_iterations` for the geometry-step limit and keeps the other
optimizer controls at their documented defaults. Unexpected exceptions
propagate; they are not hidden as ordinary embedding failures.

Stop at the **first strategy with any accepted conformers**. Return only its
accepted pool, with status `complete` or `partial`. A partial pool does not
trigger attempts to pad it with duplicates. There is no cross-attempt pooling
or post-optimization RMS deduplication. Therefore requesting N conformers does
not guarantee N distinct minima. The lowest force-field energy in this accepted
pool is recommended, with ID as a tie-breaker; with no optimization, recommend
the first accepted conformer. This is not a global minimum search.

If all configured attempts fail, return `status='failed'`, a zero-conformer
molecule and the recorded reasons. `success` means at least one accepted
conformer; it is true for a partial pool. The API never disables chirality or
uses hidden source-coordinate/template fallbacks. Experiments that disable
chirality belong in separate diagnostics, not this production API.

## Input and stereochemistry scope

Inputs may be a SMILES string or RDKit Mol. Molecules are copied, existing
conformers removed from the copy, and hydrogens expanded to vertices. Input
atom order is retained with added H atoms appended. The source object is not
modified; its existing coordinates are not used. Record kind and number of
discarded conformers are in the report.

Invalid SMILES, empty molecules, dummy atoms, disconnected fragments, enhanced
stereo groups and unsupported non-tetrahedral/non-E/Z stereo are rejected.
Relative geometries for salts/complexes require explicit user preparation;
already prepared structures can still be passed directly to `run_xtb`.

Unspecified stereochemistry is flagged in `unspecified_stereo` and `warnings`.
Set `require_assigned_stereo=True` to reject it. Otherwise the output can be
one of several allowed stereoisomers; it is not an enumeration or a unique
chemical identity. Checks protect only **specified and supported** stereo.
Embedding/optimization success cannot establish that an input SMILES inferred
from an external structure was chemically or stereochemically correct.

## Report and limits

`ConformerGenerationResult` holds `molecule` and a JSON-compatible `report`.
These are caller-owned mutable objects. The report records RDKit version,
canonical input, atom/bond stereo constraints, unresolved stereo, the retry
policy, actual seeds, embedding outcomes, optimization status/energy, checks
before and after optimization, accepted IDs and structure/geometry SHA256.
Energies are kcal/mol; positions and collision/pruning thresholds are Å.
Use `json.dump(report, ..., allow_nan=False)` for export.

`failure_counts` retains both RDKit's raw counter array and enum-indexed names,
including `UNKNOWN_SLOT_N` for counters unnamed by the installed API. These
count internal failed attempts/stages, not failed molecules, and stage names
alone do not prove a unique physical root cause. Reports distinguish an empty
embedding result from coordinate rejection, stereo rejection, missing force-field
parameters and failure to converge. With `optimize='xTB'`, each candidate also
records the native SCF/gradient optimization report; other modes perform no SCF
calculation in this layer.

Reliability here means explicit scope, bounded attempts, validation and recorded
failure. It does not guarantee coverage of every molecule, a correct protonation
state, chemical plausibility of every strained structure, thermodynamic conformer
populations or numerical invariance across RDKit versions. This protocol is
separate from the frozen 0.9.1 QM8 benchmark; improving preparation does not
retroactively change that benchmark's cohort or prediction results.

The separate [0.10 audit](validation/0.10.0/README.md) accepts 22 of the 41
previous embedding failures and 39 of 40 deterministic prior-success controls.
The remaining control fails coordinate stereo checks. Rechecking all 983
generated geometries from the old study finds nine stereo mismatches. These
results motivate validation of coordinates; they do not constitute a new
excitation-energy benchmark or prove monotonic improvement in all datasets.
