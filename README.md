# EleD

**English** | [中文](README.zh-CN.md)

**Electronic descriptors for RDKit molecules.**

EleD extends the RDKit molecular model with electronic-structure descriptors
from xTB. An RDKit molecule remains the central object: atom and bond indices,
conformers, fragments and annotations stay connected from input to report.

The package provides a Python API, an optional RDKit-style C++ extension
(`rdXTBTools`), reusable calculation results, descriptor reports and lightweight
visualization. GFN0-xTB, GFN1-xTB and GFN2-xTB are exposed through one
consistent native calculation backend.

## What it provides

- Molecular descriptors: total energy, orbital energies, HOMO/LUMO gap,
  dipole, partial charges and Wiberg bond orders.
- Atomic, bond and fragment descriptors with stable RDKit atom indexing.
- Density analysis from retained xTB data: Mulliken and Löwdin populations,
  shell populations, Mayer bonds, frontier-bond decomposition and π descriptors.
- Fixed-geometry ΔSCF reactivity: ionization and attachment energies, global
  hardness and softness, Fukui functions and local conceptual-DFT indices.
- Conformer preparation with RDKit embedding and force-field refinement, plus
  optional xTB gradient optimization.
- Versioned reports, ordered batch calculation, persistent caching and JSONL or
  CSV export.
- RDKit molecular maps and dependency-free SVG/HTML views for fragments and
  fragment interactions.
- A separate YAeHMOP/EHT adapter for comparing electronic-structure methods.

## Install

On Linux x86-64, install conda and run the commands below. The installer sets up
RDKit, EleD, its project-built xTB calculation library (`libxtb`) and the
`rdXTBTools` C++ extension:

```bash
git clone https://github.com/iawnix/EleD.git
cd EleD
python tools/install_native.py --jobs 6
conda activate "$PWD/.eled/runtime"
python -m eled doctor --check compiled --method GFN2-xTB
```

See the [installation guide](docs/install.md) for runtime layouts, source
builds and the compiled interface.

## First calculation

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import XTBCalculator, CalcXTBDescriptors

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
if AllChem.EmbedMolecule(mol, randomSeed=1) != 0:
    raise RuntimeError("could not generate a conformer")
AllChem.MMFFOptimizeMolecule(mol)

result = XTBCalculator(method="GFN2-xTB").calculate(mol)
descriptors = CalcXTBDescriptors(result)
print(descriptors["xtb_energy_ha"])
print(descriptors["xtb_homo_ev"], descriptors["xtb_lumo_ev"])
```

`XTBCalculator` accepts an RDKit molecule with a finite 3D conformer. The
returned `XTBResult` can be reused by every descriptor layer, so a report can
combine molecular, atomic, bond, fragment and orbital features from one
single-point calculation.

## Density and orbital descriptors

Retain the overlap matrix and basis maps when density-derived features are
needed:

```python
from eled import (
    XTBCalculator, CalcXTBAtomDescriptors, CalcXTBBondDescriptors,
    CalcShellPopulations, CalcDensityMatrix, CalcXTBFragmentDescriptors,
)

result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
atoms = CalcXTBAtomDescriptors(result)
bonds = CalcXTBBondDescriptors(mol, result)
shells = CalcShellPopulations(result)
density = CalcDensityMatrix(result)
fragments = CalcXTBFragmentDescriptors(result, {
    "ethyl": [0, 1],
    "oxygen": [2],
})
```

`keep_overlap=True` makes the same calculation's coefficient, overlap, AO and
shell data available to `CalcOrbitalPopulations`, `CalcFrontierDescriptors` and
`CalcDensityDescriptors`.

## Reactivity

```python
from eled import RunReactivity, XTBOptions

report = RunReactivity(
    mol,
    options=XTBOptions(method="GFN2-xTB", electronic_temperature=300.0),
    population_scheme="lowdin",
)
print(report["values"]["ip_v_ev"])
print(report["atoms"][0]["f_plus"])
```

The reactivity workflow evaluates matched charge states at the input geometry
and returns molecular and atom-level Fukui data in a versioned report.

## Reports, batches and geometry

```python
from eled import PrepareConformers, RunDescriptorBatch, ExportDescriptorReports

prepared = PrepareConformers("CCO", num_conformers=3, optimize="MMFF94")
reports = RunDescriptorBatch(
    [prepared.molecule], record_ids=["ethanol"], max_workers=1,
)
ExportDescriptorReports(reports, "features.jsonl")
```

Use `OptimizeXTB` when the descriptor protocol calls for xTB-consistent
coordinates. `RunConformerAnalysis` summarizes descriptor variation across a
conformer pool.

## Visualization

EleD uses RDKit's `rdMolDraw2D` for molecular maps. Fragment summaries and
interaction networks are returned as SVG, and `RenderFragmentAnalysis` combines
them into an embeddable HTML fragment:

```python
from eled import (
    CalcDescriptorReport, DrawDescriptorMap, PlotFragmentSummary,
    PlotFragmentInteractions, RenderFragmentAnalysis,
)

report = CalcDescriptorReport(mol, result, fragments={
    "ethyl": [0, 1],
    "oxygen": [2],
})
molecule_svg = DrawDescriptorMap(mol, report, atom_field="q_e")
summary_svg = PlotFragmentSummary(report, field="q_e")
network_svg = PlotFragmentInteractions(report)
html_fragment = RenderFragmentAnalysis(mol, report)
```

See the [visualization guide](docs/visualization.md) for rendering and report
field details.

## Documentation

- [Documentation index](docs/README.md)
- [Installation](docs/install.md)
- [Descriptor definitions](docs/descriptors.md)
- [Descriptor catalog](docs/descriptor-catalog.md)
- [Reports and reactivity](docs/descriptor-reports.md)
- [Visualization](docs/visualization.md)
- [Geometry preparation](docs/geometry-preparation.md)
- [C++ interface](docs/rdxtb-tools.md)
- [Native xTB backend build](docs/xtb-patched-build.md)
- [Project positioning](docs/project-positioning.md)

## Development

See the [development guide](docs/development.md) for working on EleD in the
environment created by the installer.

EleD is released under the BSD 3-Clause license. The bundled xTB source and
patches retain their upstream licenses.
