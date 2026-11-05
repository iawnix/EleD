# rdkit-eth

Electronic-structure descriptors for RDKit molecules through the native xTB C API.
`rdkit-eth` keeps RDKit as the molecule and conformer model, calls `libxtb` in the
same process, and returns indexed results with explicit units and provenance.

The package provides two matching interfaces:

- a pure Python/ctypes interface, `run_xtb()`;
- an optional Boost.Python module, `rdXTBTools`, with an RDKit-style `ROMol -> RunMol()` API.

Both interfaces use the same selected `libxtb` installation. The optional native
basis extension exposes the overlap matrix and AO maps from that calculation, so
population and frontier descriptors use the actual xTB basis ordering.

## What it provides

- GFN0-xTB, GFN1-xTB and GFN2-xTB gas-phase single-point calculations;
- total energy, orbital energies and occupations, charges, dipole, Wiberg bond orders,
  orbital coefficients and optional overlap/basis maps;
- molecular, atom, bond and fragment descriptors;
- Löwdin and Mulliken orbital populations and degeneracy-aware frontier descriptors;
- fixed-geometry GFN1/2 ΔSCF IP, EA, hardness, softness and condensed Fukui descriptors;
- checked SMILES-to-3D preparation with ETKDG presets and optional MMFF94 optimization;
- ordered batches, isolated workers, in-memory/disk caching and JSONL/CSV reports;
- RDKit EHT/YAeHMOP comparison helpers. EHT and xTB remain separate physical models.

Descriptor helpers consume an existing result; calculation entry points are explicit.

## Install the complete native profile

The reproducible profile targets Linux x86-64 with Python 3.11 and a conda-forge
RDKit/Boost toolchain. From a repository checkout:

```bash
python tools/install_native.py --jobs 6
conda activate "$PWD/.rdkit-eth/runtime"
python -m rdkit_eth doctor --check compiled --method GFN2-xTB --keep-overlap
```

The installer creates separate locked runtime and build prefixes, builds the
reviewed xTB source, builds the optional `rdXTBTools` module, and verifies GFN0/1/2,
overlap/population descriptors and charge-state workflows. Keep
`.rdkit-eth/runtime` and `.rdkit-eth/toolchain` together after installation.

For an offline installation, use the Linux runtime bundle from a GitHub release:

```bash
python install.py --directory /absolute/path/rdkit-eth
/absolute/path/rdkit-eth/runtime/bin/python -m rdkit_eth doctor \
  --check compiled --method GFN2-xTB --keep-overlap
```

The bundle contains the matched Python/RDKit runtime, extended `libxtb`, native
dependencies, parameters and the compiled adapter. It is a Linux x86-64 runtime
archive rather than a platform-independent wheel. See
[the installation guide](docs/install.md) and [runtime bundle guide](docs/runtime-bundle.md).

## Use an existing libxtb installation

The pure adapter can be installed into an environment that already provides RDKit,
NumPy and a loadable xTB library:

```bash
python -m pip install .
python -m rdkit_eth doctor --check native --method GFN2-xTB
```

Library discovery checks the active Python/conda prefixes and the system loader.
Use `xtb_library="/path/to/libxtb.so"` or `XTB_LIBRARY` to select a specific file.
The selected library's path and SHA-256 are recorded in result metadata.

The compiled interface needs a matching RDKit/Boost development environment:

```bash
RDKIT_ETH_BUILD_RDXTB=1 \
RDXTB_XTB_PREFIX=/path/to/extended-xtb \
python -m pip install --no-deps --no-build-isolation .
```

The compiled wheel contains the adapter and C++ module. It uses the surrounding
RDKit, Boost and xTB runtime; the complete native profile or offline bundle is the
recommended way to obtain those matched dependencies.

## First calculation

Inputs use explicit hydrogens, a finite 3D conformer and RDKit atom ordering.
`prepare_conformers()` can create a checked conformer pool from a SMILES string or
an RDKit molecule.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_xtb, calc_xtb_descriptors, annotate_mol

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
if AllChem.EmbedMolecule(mol, randomSeed=1) != 0:
    raise RuntimeError("RDKit could not generate a conformer")
AllChem.MMFFOptimizeMolecule(mol)

result = run_xtb(mol, method="GFN2-xTB")
print(result.total_energy_eh)       # Hartree
print(result.homo_ev, result.lumo_ev, result.gap_ev)  # eV
print(result.partial_charges)      # elementary charge
print(calc_xtb_descriptors(result))

annotated = annotate_mol(mol, result)
print(annotated.GetAtomWithIdx(0).GetDoubleProp("xtb_partial_charge"))
```

`run_xtb()` leaves the molecule unchanged. `annotate_mol()` returns a copy by
default and writes validated atom, bond and molecular properties to that copy.

## Core Python API

### Single points and results

```python
from rdkit_eth import XTBOptions, run_xtb

options = XTBOptions(method="GFN2-xTB", electronic_temperature=300.0)
result = run_xtb(mol, **options.to_dict())
```

`XTBResult` contains:

| Field | Unit or shape |
| --- | --- |
| `total_energy_eh` | Hartree |
| `orbital_energies_ev`, `homo_ev`, `lumo_ev`, `gap_ev` | eV |
| `orbital_occupations` | electrons per MO |
| `partial_charges` | elementary charge per atom |
| `dipole_au` | e·Bohr, `[x, y, z]` |
| `wiberg_bond_orders` | dimensionless, `[atom, atom]` |
| `orbital_coefficients` | optional `[AO, MO]` |
| `overlap_matrix` and AO/shell maps | optional native basis data |
| `metadata` | input, settings, engine and parameter provenance |

`keep_overlap=True` retains coefficients, overlap and native AO maps. It enables
population and frontier analysis and adds O(NAO²) storage.

### Molecular, atom, bond and fragment descriptors

```python
from rdkit_eth import (
    calc_xtb_descriptors, calc_xtb_atom_descriptors,
    calc_xtb_bond_descriptors, calc_xtb_fragment_descriptors,
)

molecular = calc_xtb_descriptors(result)
atoms = calc_xtb_atom_descriptors(result)
bonds = calc_xtb_bond_descriptors(mol, result)
fragments = calc_xtb_fragment_descriptors(
    result, {"donor": [0, 1], "acceptor": [2, 3, 4]}
)
```

The scalar descriptor set includes frontier energies, dipole components and norm,
charge statistics, electron/orbital counts and Wiberg sums. Atom rows contain
charges and Wiberg sums; bond rows contain indexed Wiberg orders; fragment rows
contain charge and internal/external Wiberg summaries.

### Populations and frontier subspaces

```python
from rdkit_eth import calc_orbital_populations, calc_frontier_descriptors

result = run_xtb(mol, method="GFN2-xTB", keep_overlap=True)
lowdin = calc_orbital_populations(result, scheme="lowdin")
mulliken = calc_orbital_populations(result, scheme="mulliken")
frontier = calc_frontier_descriptors(mol, result)
```

These are spin-summed, ground-state descriptors. Near-degenerate HOMO/LUMO
orbitals are grouped using the configured energy tolerance; the output records the
indices and settings used for each group.

### Fixed-geometry reactivity

```python
from rdkit_eth import XTBOptions, run_reactivity

report = run_reactivity(
    mol,
    options=XTBOptions(method="GFN2-xTB"),
    population_scheme="lowdin",
    return_errors=True,
)
print(report["values"]["vertical_ip_eV"])
print(report["atoms"])
```

The protocol evaluates the reference, ionized and attached charge states at the
same geometry. GFN1/2 support vertical ΔSCF IP/EA and atomic/fragment Fukui
quantities using charge, Mulliken or Löwdin populations. Open-shell references
accept explicit adjacent-state spin settings.

### Batches, caching and reports

```python
from rdkit_eth import XTBCache, run_descriptor_batch, export_descriptor_reports

reports = run_descriptor_batch(
    molecules,
    record_ids=record_ids,
    cache=XTBCache(),
    max_workers=4,
    isolate=True,
)
export_descriptor_reports(reports, "features.jsonl")
```

Batch APIs preserve input order. Isolated workers provide per-molecule timeout and
crash boundaries while the parent process owns cache/checkpoint writes. Use
`XTBDiskCache` and `ReactivityCheckpoint` for datasets that span process restarts.
Reports retain units, settings, quality flags, missing fields, errors and engine
hashes. See [descriptor reports](docs/descriptor-reports.md),
[persistent batches](docs/persistent-batches.md) and
[isolated batches](docs/isolated-batches.md).

## RDKit-style C++/Python interface

`rdXTBTools` accepts an RDKit `ROMol` and returns `(success, XTBResults)`:

```python
from rdkit_eth import rdXTBTools

ok, result = rdXTBTools.RunMol(
    mol,
    method="GFN2-xTB",
    keepOverlap=True,
)
if not ok:
    raise RuntimeError(result.errorMessage)
print(result.totalEnergyHartree)       # Hartree
print(result.GetOrbitalEnergies())     # eV
print(result.GetOverlapMatrix())       # AO × AO
```

The C++ library is available as `RDKitXTB::RDKitXTB` for applications that use
RDKit's CMake package. See [the C++ API guide](docs/rdxtb-tools.md).

## Conformer preparation

```python
from rdkit_eth import prepare_conformers

prepared = prepare_conformers(
    "C[C@H](O)C(=O)O",
    num_conformers=4,
    method="ETKDGv4",
    fallback_methods=("ETKDGv3", "srETKDGv3"),
    optimize="MMFF94",
)
if not prepared.success:
    raise RuntimeError(prepared.report)
mol = prepared.molecule
```

The preparation layer records embedding attempts, stereo checks, optimization
status and coverage in `prepared.report`. It is a checked input layer for xTB
workflows, while conformer selection remains under the caller's control.

## Scope and scientific interpretation

The current release focuses on gas-phase GFN0/1/2 single points and ground-state
electronic descriptors. Energies, orbital gaps, charges, populations and Fukui
values retain their method-specific definitions and are suitable for comparative
workflows with recorded provenance. HOMO-LUMO gaps are frontier descriptors, and
vertical ΔSCF values are approximate GFN energy differences. Solvent models,
geometry optimization and excited-state solvers can be added as separate workflow
layers; the package reports the ground-state quantities returned by xTB.

## Development

```bash
conda activate "$PWD/.rdkit-eth/runtime"
python -m pip install --no-deps --no-build-isolation -e .
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q
```

The repository includes the locked build recipes, C++ tests, integration checks,
examples and validation datasets. Start with the focused guides:

- [installation](docs/install.md)
- [descriptor definitions](docs/descriptors.md)
- [descriptor reports and reactivity](docs/descriptor-reports.md)
- [compiled interface](docs/rdxtb-tools.md)
- [patched xTB build](docs/xtb-patched-build.md)
- [validation scope](docs/mainline-validation.md)

## License

`rdkit-eth` is released under the BSD 3-Clause license. xTB source and patches
retain their upstream licensing terms; see the bundled license files and
[the xTB build guide](docs/xtb-patched-build.md).
