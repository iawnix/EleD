# Descriptor reports and vertical reactivity (1.1)

These APIs organize libxtb results into reproducible features. They do not add
a new electronic Hamiltonian, excited-state solver or fitted IP/EA model.
Both ctypes and compiled results work with the extraction functions. There
are no additional Python dependencies.

## Numerical settings

| Python / `XTBOptions` | `rdXTBTools.RunMol` | Default / constraint |
| --- | --- | --- |
| `max_iterations` | `maxIterations` | 250; positive 32-bit integer |
| `accuracy` | `accuracy` | 1.0; finite, 1e-4 to 1e3; smaller is tighter |
| `electronic_temperature` | `electronicTemperature` | 300 K; finite, nonnegative |

These are xTB numerical controls, not accuracy estimates. Booleans are rejected.
Values outside the supported accuracy range are rejected instead of silently
clamped. xTB floors temperature at 1e-6 K; metadata records both requested and
effective temperature. Settings participate in cache identity. GFN0 is non-SCF;
retaining its settings does not imply an SCC iteration or finite-T SCF protocol.

## Report API and export

```python
from eled import (
    XTBCalculator, CalcDescriptorReport, RunDescriptorBatch,
    ExportDescriptorReports, XTBOptions, XTBCache,
)

result = XTBCalculator(accuracy=0.1).calculate(mol, keep_overlap=True)
report = CalcDescriptorReport(mol, result, record_id="molecule-1")
reports = RunDescriptorBatch(molecules, record_ids=identifiers, cache=XTBCache())
ExportDescriptorReports(reports, "features.csv")
ExportDescriptorReports(reports, "features.jsonl")
```

`CalcDescriptorReport` performs no engine calls. It checks original atom order,
topology and geometry before attaching atom/bond/spatial descriptors. It combines
the [existing descriptors](descriptors.md), both population definitions and
degeneracy-aware frontier distributions. `RunDescriptorBatch` retains overlap by default and supports
[isolated parallel execution](isolated-batches.md). The same extended library
is used for basic and overlap-aware features; set `keep_overlap=False` only
when the additional basis snapshot is not needed.

The report schema version is **2** and descriptor-definition version is **1.1**;
neither is the package version. `DescriptorColumns(frontier_depth=2)` returns
the ordered scalar schema and units. Each report includes:

- `values`, `units`: fixed molecular/frontier scalar fields. Numeric units are
  explicit; nonnumeric identity fields have no physical unit.
- `atoms`, `bonds`, `fragments`, `populations`, `frontier`: nested indexed data.
  Frontier sections retain the actual MO indices used for each group.
- `density`, `shell_populations`, `mayer_bonds`, `frontier_bonds`, `pi`: zero-cost
  sections derived from retained `P`, `S`, `C` and shell maps. Each section
  carries units, invariance and explicit missing-value reasons.
- `analysis_settings`: thresholds, degeneracy tolerance and extraction controls.
- `provenance`: input hashes, charge/spin, numerical settings, engine and parameter
  identity. Failed calculations retain requested options where no result exists.
- `quality`, `missing`, `error`: observed flags, missing reasons and errors.

Statuses are `ok`, `partial`, `warning`, `error`; flags and `missing` remain the
authoritative detail when issues coexist. Unretained basis gives null population
sections, not zero populations. Invalid matrices or mismatched provenance raise
on extraction; batching converts a row failure into null scalar fields and an
error record. IDs and input order are preserved. Global analysis-argument errors
are rejected before calculations. Fragment indices are shared across batch rows;
for different fragment definitions, extract each report separately.

CSV contains scalar features and JSON columns for units/settings/provenance/
quality/missing/error. Null numeric features become blank CSV cells. JSONL keeps
all nested atom, bond, fragment and orbital arrays. Export rejects incompatible
schema/descriptor versions, report types or columns and nonfinite JSON values,
then atomically replaces the destination. Export single-point and reactivity
reports into separate files.

Quality checks include fractional occupations, frontier-threshold sensitivity,
charge/electron conservation, missing frontier groups, spin-summed open-shell
interpretation and charged-dipole origin dependence. The [SCC diagnostics
extension](scf-diagnostics.md) reports the actual iteration count, final energy
change and mixed-variable RMS with their thresholds. The runtime rejects
libraries without this extension; GFN0 is non-SCF.
Convergence criteria are not an experimental accuracy guarantee.

## Fixed-geometry charge-state analysis

```python
from eled import RunReactivity, XTBOptions, CalcReactivityDescriptors

# Default protocol: GFN2, 300 K, reference wavefunction copies, fixed geometry.
reaction = RunReactivity(mol, population_scheme="lowdin", return_errors=True)
print(reaction["values"]["vertical_ip_ev"])
print(reaction["atoms"])  # f_plus, f_minus, f_radical, dual

# Explicit independent initialization for comparison; no automatic selection.
reaction = RunReactivity(mol, options=XTBOptions(electronic_temperature=300.),
                          interface="compiled", initialization="independent", return_errors=True)
# Existing compatible results can instead be analyzed without new calculations:
# CalcReactivityDescriptors(mol, reference, ionized, attached)
```

Each atom row also contains `local_softness_per_ev = s f+`,
`local_electrophilicity_ev = ω f+`, and the explicitly named
`local_nucleophilicity_per_ev = s f−` proxy. Their units and missing-value
semantics are in `atom_units`, `atom_invariance`, and
`atom_missing_semantics`; no absolute nucleophilicity reference is assumed.

Only GFN1/2 support this ΔSCF workflow. It evaluates reference charge Q,
ionized Q+1 and attached Q−1 at the **same input geometry**, method, numerical
settings, engine and parameter identity. Inputs require explicit H atoms and
a 3D conformer, as for `XTBCalculator`. No optimization or spin-state search occurs.
For a closed-shell reference, adjacent states default to one unpaired electron.
For an open-shell reference, explicitly supply both `ionized_uhf` and
`attached_uhf`; spin parity is validated. The selected states are not asserted
to be ground states.

Since 0.7, `RunReactivity()` and `XTBOptions()` both default to 300 K.
Explicit `XTBOptions(electronic_temperature=0.)` still requests the 1e-6 K
engine floor. Temperature is an electronic occupation setting, not a nuclear
thermal ensemble or a solution free energy. Finite-T results retain their flag.

The default `initialization="reference"` computes the reference independently,
then supplies each adjacent state with its **own deep copy** via the official
`xtb_copyResults` API. GFN1/2 resets electron count and spin from the target
molecule during SCF. The ionized state is never used as the attached-state
initial guess. Ctypes and compiled paths use the same protocol; no engine patch
or additional dependency is needed. All restart handles are released within the
triplet workflow, including failures. Returned results contain owned arrays.

`initialization="independent"` uses fresh guesses for every state, reproducing
0.6 initialization. To reproduce the complete 0.6 default protocol, also set
`electronic_temperature=0.` explicitly. Reports record `protocol_version="2.0"`
and requested initialization in `analysis_settings`; successful state provenance
records the actual initialization. Descriptor version 1.1 and schema 2 identify
the added nested density sections; existing scalar columns remain unchanged.

Reference reuse can improve convergence, but neither initialization guarantees
the lowest SCF solution. Compare protocols explicitly if this matters; the API
does not automatically retry or choose the lower result. The CLI's automatic
open-shell spin choice remains distinct from our explicit spin policy.

With N modeled electrons and energies converted from Hartree to eV:

| Field | Definition |
| --- | --- |
| `vertical_ip_ev` | E(N−1) − E(N) |
| `vertical_ea_ev` | E(N) − E(N+1) |
| `chemical_potential_ev` | −(IP+EA)/2 |
| `electronegativity_ev` | (IP+EA)/2 |
| `hardness_ev` | η = (IP−EA)/2 |
| `softness_per_ev` | 1/(IP−EA) = 1/(2η) |
| `electrophilicity_ev` | μ²/(2η) |

The softness convention is explicit because literature conventions differ.
When η ≤ 1e-8 eV, softness/electrophilicity are null with a reason. Negative IP
or EA is retained and flagged; a converged attached state does not establish a
bound anion. These values use the libxtb total-energy API and are approximate
GFN descriptors, **not calibrated experimental IP/EA or optical excitations**.

For atomic electron populations nA, f+ = nA(N+1)−nA(N), f− = nA(N)−nA(N−1),
f0 = (f++f−)/2, dual = f+−f−. `population_scheme="charges"` uses equivalent
differences of engine charges: f+ = qA(N)−qA(N+1), f− = qA(N−1)−qA(N).
`"mulliken"` and `"lowdin"` require the extended library and retain basis data
automatically. Signed values are preserved; sums are checked against 1, 1 and
0. Explicit disjoint `fragments` aggregate these quantities. These finite
charge differences are distinct from frontier-MO population approximations.

Any charge-state calculation failure stops subsequent states. Runtime/setup,
allocation/extraction and unexpected errors propagate even with `return_errors=True`
(since 0.8); only genuine single-point failures receive error reports. By default this
raises `ChargeStateCalculationError`, with `.state` and `.states` diagnostics.
`return_errors=True` returns null derived features and retains completed states,
failed-state diagnostics, requested settings/spins and `not_run` states. Invalid
inputs, unsupported options and engine resolution errors before the state loop,
or incompatible results during analysis, still raise; this flag is not a blanket
exception suppressor. Dataset callers should handle those errors explicitly.

There is no automatic change of temperature, spin or method to force convergence.
In the archived 0.6 independent-initialization tests, cold GFN1 ethanol and GFN2 ammonium charge states
failed even at 1000 iterations. Explicit 300 K protocols converged, with different
energies. The first failing state can depend on numerical details; all failure
rows must be retained in dataset accounting.

## Reproduce checks

```bash
python -m eled doctor --check compiled --keep-overlap --reactivity
python examples/descriptor_reports.py --output build/reports --sensitivity
```

The example saves exact geometries, settings, all three state records and failure
rows for a small temperature study. Its unoptimized structures test protocol
behavior, not predictive accuracy. Conformer sensitivity and a QM8 descriptor pilot are now available separately.
Solvent support and calibrated IP/EA validation remain outside this protocol.


The `charges` and `mulliken` options for gas-phase GFN1/2 are two extraction/
validation paths for the same population definition, rather than independent
new descriptors. Löwdin is a distinct, established partition. Original xTB
already supports atomic charge-difference Fukui via `--vfukui`; our contribution
is the RDKit/libxtb workflow, analysis options, fragment mapping and provenance.
For an explicit development comparison, run `tools/compare_fukui_cli.py` with
`--executable`, `--library` and a new `--output` directory. This diagnostic keeps
CLI logs and frozen coordinates and compares both initialization policies.
It is not a production calculation backend or a package dependency. Original
CLI Fukui values print to three decimals, so matching those values only proves
agreement within the reported tolerance, not full-precision numerical identity.

For resumable datasets, per-molecule options/spins and persistent caching, see
[persistent batches](persistent-batches.md).
