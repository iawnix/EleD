# rdkit-eth

An RDKit adapter for in-process xTB electronic-structure calculations. It
validates indexed molecular input, returns quantities with explicit units,
and maps results back to the original atoms and bonds. RDKit's YAeHMOP
`rdEHTTools` interface remains available separately.

Version 0.10 supports **gas-phase single points through libxtb only**.
The tblite backend and Python extra have been removed. Our reviewed native
extension exports overlap and basis maps from the same calculation. The external xTB CLI backend and `auto`
fallback were removed in 0.2. See the [0.2 migration notes](docs/migration-0.2.md).

New in 0.9: [conformer sensitivity](docs/conformer-sensitivity.md) with explicit
coverage/failures and invariant statistics, plus a reproducible
[QM8 excitation-energy pilot](docs/qm8-pilot.md). The pilot finds a clear benefit
from basic xTB features on its fixed cohort; additional orbital distribution
features show only small, split-dependent gains. Machine-learning tools remain
optional developer dependencies in a separate environment.

0.9.1 fixes automatic libxtb discovery in venvs inheriting the standard runtime.
A [matched geometry study](docs/qm8-geometry-study.md) records exploratory
SMILES-generated geometry results, failures and paired uncertainty intervals.
The subsequent coordinate audit found nine stereo mismatches in its historical
cohort; that limitation is recorded alongside the original results.

## Installation

0.10 adds [`prepare_conformers`](docs/geometry-preparation.md): deterministic
RDKit retries with coordinate-derived tetrahedral/E/Z checks before and after
optional force-field optimization. It returns an explicit-H molecule ready for
the existing xTB APIs, plus complete provenance and failure/coverage records.

The standard profile is Linux x86-64 / Python 3.11. With conda available, run:

```bash
python tools/install_native.py --jobs 6
conda activate "$PWD/.rdkit-eth/runtime"
python -m rdkit_eth doctor --check compiled --keep-overlap
```

This installs matched RDKit dependencies, builds our extended libxtb and the
C++/Python adapter, then checks GFN0/1/2 and population descriptors through
both interfaces. The shared library and parameters live in the Python runtime
prefix, so ordinary use needs no separate xTB activation or executable lookup.
Keep `.rdkit-eth/runtime` and `.rdkit-eth/toolchain` in place; the build is tied
to these local native-library prefixes.

The pinned source includes the official [GFN0 repair, PR #1490](https://github.com/grimme-lab/xtb/pull/1490).
Only our native basis-export patch is applied. Previously verified 0.4 engines
remain usable explicitly; official libraries retain basic GFN1/2 compatibility.

See the [installation guide](docs/install.md) for upgrades, manual builds and
artifact checks. This is a local source installer; plain pip installation does
not install libxtb, and no PyPI/conda-forge full-runtime package is published.

## RDKit-style C++/Python interface

The optional **`rdXTBTools` module** provides a direct `ROMol` → libxtb
interface in RDKit's `RunMol()` style:

```python
from rdkit_eth import rdXTBTools
ok, result = rdXTBTools.RunMol(mol, method="GFN2-xTB")
if not ok:
    raise RuntimeError(result.errorMessage)
print(result.totalEnergyHartree)
print(result.GetOrbitalEnergies())  # eV
print(result.GetWibergBondOrders())  # dimensionless; not EHT overlap populations
```

Build it with `RDKIT_ETH_BUILD_RDXTB=1` in the matched RDKit/Boost development
environment. It supports gas-phase GFN0/1/2, requires a verified patched
library for GFN0, and also provides an installable C++ library. See the
[build guide, API contract and unit comparison](docs/rdxtb-tools.md).

## Versioned features and vertical reactivity

Version 0.6 adds reports with units, provenance, quality flags and CSV/JSONL
export, plus fixed-geometry GFN1/2 ΔSCF IP/EA, conceptual-DFT descriptors and
charge/Mulliken/Löwdin condensed Fukui analysis. Numerical settings are shared
by the ctypes and compiled APIs. No additional Python backend is required. Reactivity now defaults to 300 K and
separate copies of the reference wavefunction for both adjacent charge states;
`initialization="independent"` retains the fresh-guess protocol. Open-shell
references still require explicit adjacent spins.

```python
from rdkit_eth import run_descriptor_batch, export_descriptor_reports, run_reactivity

reports = run_descriptor_batch(molecules, record_ids=identifiers)
export_descriptor_reports(reports, "features.csv")
reactivity = run_reactivity(mol, population_scheme="lowdin", return_errors=True)
print(reactivity["values"])  # missing/failed values are None, never fabricated zeros
```

Charge-state calculations can fail to converge; errors retain settings and state
diagnostics. These are approximate GFN energy differences, not calibrated IP/EA
predictions or optical excitation energies. See the [protocol and API guide](docs/descriptor-reports.md)
and `examples/descriptor_reports.py` for complete examples and spin conventions.

## Persistent datasets

Version 0.8 adds `ReactivityJob`, `run_reactivity_batch`, `XTBDiskCache` and
`ReactivityCheckpoint`. Completed molecules resume across process restarts;
failed rows retain attempt history and require explicit retry. Cache identity
includes geometry, charge/spin, numerical settings, initialization, analysis,
engine/parameter hashes and adapter implementation. Configuration, installation
and storage errors abort instead of being swallowed as molecular failures.

See [persistent batch usage](docs/persistent-batches.md) and
`examples/reactivity_batch.py`. The frozen regression corpus now includes 37
molecules plus difficult calculation protocols, with all failures retained.

## Single molecules and results

The pure Python adapter uses ctypes to call the same native library:

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_xtb, calc_xtb_descriptors, annotate_mol

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
result = run_xtb(mol, method="GFN2-xTB")  # backend="native" by default
print(calc_xtb_descriptors(result))
annotated = annotate_mol(mol, result)
print(annotated.GetAtomWithIdx(0).GetDoubleProp("xtb_partial_charge"))
print(annotated.GetBondWithIdx(0).GetDoubleProp("xtb_wiberg_order"))
```

`run_xtb` never changes the molecule. A missing/incompatible library or failed
calculation raises an error. Geometry optimization, solvation and GFN-FF are
outside the current API. Calculations do not create temporary work directories.

Native GFN0 requires the build's capability record, the reviewed upstream fix identity,
and a matching SHA-256 of the selected library. An xTB/API version alone cannot
identify the fix. Missing or mismatched verification raises an error before
GFN0 initialization. See the [root-cause analysis](docs/gfn0-native-state.md).
After activating the standard installation:

```bash
conda activate "$PWD/.rdkit-eth/runtime"
```

```python
result = run_xtb(mol, method="GFN0-xTB")
print(result.metadata["engine"]["native_gfn0"])
```

Select a library explicitly with `xtb_library="/prefix/lib/libxtb.so"`.
Discovery otherwise uses `XTB_LIBRARY`, the active Python installation,
`CONDA_PREFIX`, `XTB_PREFIX`, and the system library loader. An explicit path
is respected. Static `libxtb.a` cannot be loaded with ctypes.

Parameter files are located using `XTBPATH`, then `share/xtb` associated with
the selected library. Located parameters are used explicitly and recorded by
path and SHA-256; installed text parameters can differ slightly from compiled-in
parameters. GFN0 requires its parameter file; GFN1/2 can use built-in parameters.

Inputs require a conformer marked as 3D, finite coordinates and separate
hydrogen atoms. Hydrogen counts such as `[NH4+]` must also be expanded with
`Chem.AddHs`. Charge and unpaired-electron count (`uhf`) default to the sum of
RDKit formal charges and radical electrons. Their electron parity is checked;
provide `uhf` explicitly when the RDKit radical count does not describe the
intended spin state. The adapter does not choose protonation or spin states.

`XTBResult` provides:

- Orbital energies in eV, occupations, atomic charges in elementary charges,
  dipole in atomic units, Wiberg bond orders, and total energy in Hartree.
- Metadata with the resolved conformer ID, indexed structure and geometry
  hashes, charge, spin, backend, engine/parameter identity and elapsed time.
- Native orbital coefficients with axes `[AO][MO]`; with `keep_overlap=True`,
  the verified extension also returns overlap and zero-based AO/shell/atom maps.

`homo_ev` and `lumo_ev` use an occupation threshold of `1e-7`; `gap_ev` uses
their difference when available. `fermi_energy_ev` is a HOMO-based estimate,
not the electronic chemical potential. None of these gaps is an optical
excitation energy. The native C API reports an ABI version separately from
`result.version`, which uses a conda package version when identifiable and is
otherwise `None`. The binary hash identifies the selected engine.

`annotate_mol` returns a copy by default; `copy=False` modifies the input after
validation. The original indexed topology and input conformer must match the
result. Existing conformers and bond types are preserved. WBO properties are
attached only to existing RDKit bonds; nonbonded pair values remain in the
result matrix. Atom/bond properties are not automatically portable SDF fields.

## Electronic descriptor interfaces

The descriptor API provides molecular dipole/charge statistics, neighboring frontier
energies, atomic WBO sums and explicitly indexed fragment descriptors. The
common scalar/atom/fragment extractors accept both `XTBResult` and compiled
`rdXTBTools.XTBResults` without rerunning a calculation.

```python
from rdkit_eth import calc_xtb_atom_descriptors, calc_xtb_fragment_descriptors
atoms = calc_xtb_atom_descriptors(result)
fragments = calc_xtb_fragment_descriptors(result, {"site_a": [0], "site_b": [1, 2]})
```

For nonorthogonal orbital populations, select the verified extended libxtb
installation and retain its overlap and AO mapping:

```python
from rdkit_eth import calc_orbital_populations, calc_frontier_descriptors
orbital_result = run_xtb(mol, keep_overlap=True)
populations = calc_orbital_populations(orbital_result, scheme="lowdin")
frontier = calc_frontier_descriptors(mol, orbital_result)
```

This adds Löwdin/Mulliken populations, degeneracy-aware frontier subspaces,
localization measures and atom-centroid spatial proxies. GFN0/1/2 and
radicals are supported with spin-summed populations. Official libxtb without
the extension cannot retain overlap; requesting it raises a clear error. These are ground-state descriptors, not excited-state
observables. See the [complete API, definitions and limitations](docs/descriptors.md)
and [runnable example](examples/electronic_descriptors.py).

## Native basis extension

The project-local C API exports the actual overlap matrix and basis maps in
libxtb's native AO order. No second engine supplies or reorders these data.
`keep_overlap=True` (compiled: `keepOverlap=True`) also retains MO coefficients.
Matrix storage is opt-in and costs O(NAO²). Both interfaces support the same
population and molecule-mapping helpers. See [native build details](docs/xtb-patched-build.md).

Historical tblite comparisons remain as validation history, with their tools
archived under `tools/diagnostics/archive/tblite-0.3/`. They are not a runtime
dependency. Removing tblite does not make the xTB Hamiltonian a new method;
our contribution is the native data interface, provenance, descriptor definitions
and validation.

## Batch calculations and caching

```python
from rdkit_eth import XTBOptions, XTBCache, run_xtb_batch

options = XTBOptions(method="gfn2")
cache = XTBCache()
results = run_xtb_batch([mol, mol], options=options, cache=cache,
                        return_exceptions=True)
# Successful entries are XTBResult objects; failures are Exception objects.
# options.to_dict() also works with run_xtb(mol, **...).
```

Batches execute sequentially and preserve input order; only `max_workers=1`
is supported. Native calls from Python threads are serialized within this
package, including calls through `rdXTBTools`. Other independent bindings in
the same process are outside that lock. No wall-clock timeout is enforced.

The cache is in memory. Keys include indexed structure, coordinates, resolved
charge/spin/conformer, settings, backend, engine identity, binary/parameter
hashes and relevant environment settings. Concurrent identical cache requests
share a calculation; errors are not stored. Returned metadata is copied.
Cached elapsed times refer to the original calculation. Native caching needs
a library file path that can be hashed. Clear the cache after changing external
parameters/dependencies; restart Python after replacing a loaded library.

## EHT comparison

```python
from rdkit_eth import calc_orbital_descriptors
features = calc_orbital_descriptors(mol, backend="both")
```

EHT and xTB fields have different prefixes. Their physical models and energy
references differ; the EHT total energy uses eV and xTB uses Hartree.

## Development and validation

```bash
conda activate "$PWD/.rdkit-eth/runtime"
python -m pip install --no-deps --no-build-isolation -e .
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q
python examples/rdkit_workflow.py
```

Integration tests use `XTB_LIBRARY` or normal library discovery. They skip
when required engines are absent; `python -m pytest -m 'not integration'`
runs without xTB. Shared-library tests explicitly disable subprocess execution.
See [validation notes](docs/validation.md) for the locally verified scope.

`tools/build_patched_xtb.sh PREFIX` fetches the pinned xTB commit, applies
the native basis-export patch on the upstream-fixed source, builds libxtb, installs parameters
and runs upstream and installed-library tests.
The upstream build also produces an executable used for build provenance;
it is not a calculation backend in this package. See the
[patched-build guide](docs/xtb-patched-build.md). Historical diagnostic tools
retain their upstream executable reproductions.
