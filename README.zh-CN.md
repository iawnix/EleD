# EleD

[English](README.md) | **中文**

**EleD 是面向 RDKit 分子的电子结构描述符扩展。**

它以 RDKit 分子对象为中心，将原子、键、构象和片段的索引贯穿到 xTB
计算、描述符分析和结果报告中。项目提供 Python 接口、可选的
RDKit 风格 `rdXTBTools` C++ 扩展、可复用的计算结果、报告和轻量级可视化。

GFN0-xTB、GFN1-xTB 和 GFN2-xTB 通过统一的原生计算后端执行。

## 功能

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

## 安装

在 Linux x86-64 上准备好 conda，然后执行以下命令。安装器会统一安装
RDKit、EleD、项目构建的 xTB 计算库（`libxtb`）和 `rdXTBTools` C++ 扩展：

```bash
git clone https://github.com/iawnix/EleD.git
cd EleD
python tools/install_native.py --jobs 6
conda activate "$PWD/.eled/runtime"
python -m eled doctor --check compiled --method GFN2-xTB
```

详细步骤见[安装说明](docs/zh-CN/install.md)。

## 最小示例

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

需要轨道、密度或 shell 描述符时，在计算中设置 `keep_overlap=True`，再
调用相应的 `Calc...` 函数。需要反应性分析时使用 `RunReactivity`。

## 密度和轨道描述符

计算密度相关特征时，保留重叠矩阵和基组映射：

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

`keep_overlap=True` 会保留同一次计算中的轨道系数、重叠矩阵、AO 和 shell
映射，供 `CalcOrbitalPopulations`、`CalcFrontierDescriptors` 和
`CalcDensityDescriptors` 使用。

## 反应性分析

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

该流程在输入几何上计算匹配的电荷态，返回分子级反应性指标和原子级
Fukui 数据。

## 报告、批量计算与构象

```python
from eled import PrepareConformers, RunDescriptorBatch, ExportDescriptorReports

prepared = PrepareConformers("CCO", num_conformers=3, optimize="MMFF94")
reports = RunDescriptorBatch(
    [prepared.molecule], record_ids=["ethanol"], max_workers=1,
)
ExportDescriptorReports(reports, "features.jsonl")
```

需要 xTB 优化坐标时使用 `OptimizeXTB`；需要汇总多个构象的描述符变化时，
使用 `RunConformerAnalysis`。

## 可视化

EleD 使用 RDKit 的 `rdMolDraw2D` 绘制分子图，以 SVG 返回片段汇总和
相互作用网络；`RenderFragmentAnalysis` 将这些视图组合成可嵌入的 HTML：

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

绘图选项和报告字段见[可视化指南](docs/zh-CN/visualization.md)。

## 文档

- [文档目录](docs/README.md)
- [安装说明](docs/zh-CN/install.md)
- [描述符定义](docs/zh-CN/descriptors.md)
- [描述符目录](docs/zh-CN/descriptor-catalog.md)
- [报告与反应性](docs/zh-CN/descriptor-reports.md)
- [可视化](docs/zh-CN/visualization.md)
- [构象准备](docs/zh-CN/geometry-preparation.md)
- [C++ 接口（英文）](docs/rdxtb-tools.md)
- [原生 xTB 后端构建（英文）](docs/xtb-patched-build.md)
- [项目定位](docs/zh-CN/project-positioning.md)

## 开发

在完整安装创建的环境中修改 EleD，具体步骤见[开发指南](docs/zh-CN/development.md)。

EleD 采用 BSD 3-Clause 许可证；随项目提供的 xTB 源码和补丁遵循其上游许可证。
