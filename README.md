# EleD

**Electronic descriptors for RDKit molecules.**

EleD extends the RDKit molecular model with electronic-structure descriptors
from xTB. An RDKit molecule remains the central object: atom and bond indices,
conformers, fragments and annotations stay connected from input to report.

The package provides a Python API, an optional RDKit-style C++ extension
(`rdXTBTools`), reusable calculation results, descriptor reports and lightweight
visualization. GFN0-xTB, GFN1-xTB and GFN2-xTB are exposed through one
consistent native calculation backend.

[English](#english) · [中文](#中文)

## English

### What it provides

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

### Install

The standard installer builds a matched RDKit environment, the xTB calculation
backend and the optional C++ extension:

```bash
git clone https://github.com/iawnix/EleD.git
cd EleD
python tools/install_native.py --jobs 6
conda activate "$PWD/.eled/runtime"
python -m eled doctor --check compiled --method GFN2-xTB
```

For an existing EleD-compatible native xTB backend, install the Python layer
directly:

```bash
python -m pip install .
python -m eled doctor --check native --method GFN2-xTB
```

See the [installation guide](docs/install.md) for runtime layouts, source
builds and the compiled interface.

### First calculation

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

### Density and orbital descriptors

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

### Reactivity

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

### Reports, batches and geometry

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

### Visualization

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

### Documentation

- [Documentation index](docs/README.md)
- [Installation](docs/install.md) · [安装说明](docs/zh-CN/install.md)
- [Descriptor definitions](docs/descriptors.md) · [描述符定义](docs/zh-CN/descriptors.md)
- [Descriptor catalog](docs/descriptor-catalog.md) · [描述符目录](docs/zh-CN/descriptor-catalog.md)
- [Reports and reactivity](docs/descriptor-reports.md) · [报告与反应性](docs/zh-CN/descriptor-reports.md)
- [Visualization](docs/visualization.md) · [可视化](docs/zh-CN/visualization.md)
- [Geometry preparation](docs/geometry-preparation.md) · [构象准备](docs/zh-CN/geometry-preparation.md)
- [C++ interface](docs/rdxtb-tools.md)
- [Native xTB backend build](docs/xtb-patched-build.md)
- [Project positioning](docs/project-positioning.md) · [项目定位](docs/zh-CN/project-positioning.md)

### Development

```bash
conda activate "$PWD/.eled/runtime"
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest -q
```

EleD is released under the BSD 3-Clause license. The bundled xTB source and
patches retain their upstream licenses.

## 中文

### 项目简介

**EleD 是面向 RDKit 分子的电子结构描述符扩展。**

它以 RDKit 分子对象为中心，将原子、键、构象和片段的索引贯穿到 xTB
计算、描述符分析和结果报告中。项目提供 Python 接口、可选的
RDKit 风格 `rdXTBTools` C++ 扩展、可复用的计算结果、报告和轻量级可视化。

GFN0-xTB、GFN1-xTB 和 GFN2-xTB 通过统一的原生计算后端执行。

### 功能

- 分子级描述符：总能量、轨道能量、HOMO/LUMO 能隙、偶极矩、部分电荷和
  Wiberg 键级；
- 原子级、键级和片段级描述符，保持 RDKit 原子编号；
- 密度分析：Mulliken/Löwdin 布居、shell 布居、Mayer 键级、前线键分解和
  π 描述符；
- 固定几何下的 ΔSCF 反应性分析、电离/附着能、硬度、软度、Fukui 函数和
  局部概念 DFT 指数；
- RDKit 构象生成、MMFF94/UFF 优化和可选的 xTB 梯度优化；
- 版本化报告、批量计算、缓存以及 JSONL/CSV 导出；
- RDKit 分子图、片段汇总和片段相互作用网络可视化；
- 独立的 YAeHMOP/EHT 对照适配器。

### 安装

标准安装会构建匹配的 RDKit 环境、xTB 计算后端和可选 C++ 扩展：

```bash
git clone https://github.com/iawnix/EleD.git
cd EleD
python tools/install_native.py --jobs 6
conda activate "$PWD/.eled/runtime"
python -m eled doctor --check compiled --method GFN2-xTB
```

已有兼容的原生 xTB 后端时，可以直接安装 Python 层：

```bash
python -m pip install .
python -m eled doctor --check native --method GFN2-xTB
```

详细步骤见[安装说明](docs/zh-CN/install.md)。

### 最小示例

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import XTBCalculator, CalcXTBDescriptors

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
AllChem.MMFFOptimizeMolecule(mol)

result = XTBCalculator(method="GFN2-xTB").calculate(mol)
descriptors = CalcXTBDescriptors(result)
print(descriptors["xtb_gap_ev"])
```

需要轨道、密度或 shell 描述符时，在计算中设置 `keep_overlap=True`，再
调用相应的 `Calc...` 函数。需要反应性分析时使用 `RunReactivity`。

### 文档与开发

从[双语文档目录](docs/README.md)开始查看描述符定义、报告格式、可视化、
构象准备和 C++ 接口。

```bash
conda activate "$PWD/.eled/runtime"
python -m pip install --no-deps --no-build-isolation -e .
python -m pytest -q
```

EleD 采用 BSD 3-Clause 许可证；随项目提供的 xTB 源码和补丁遵循其上游许可证。
