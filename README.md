# rdkit-eth

**RDKit-native electronic-structure descriptors powered by an extended xTB C API.**

`rdkit-eth` keeps RDKit as the molecule, atom-order and conformer model, then
uses one project-built `libxtb` installation for gas-phase GFN0, GFN1 and GFN2
calculations. It exposes the result through Python/ctypes and an optional
RDKit-style C++/Boost.Python module. Units, atom indexing, provenance and
missing-value reasons are part of the public API.

[English](#english) · [中文](#中文)

## English

### What it provides

- GFN0-xTB, GFN1-xTB and GFN2-xTB single-point calculations;
- molecular, atom, bond and fragment descriptors;
- HOMO/LUMO, dipole, partial-charge and Wiberg-bond features;
- same-calculation Mulliken/Löwdin orbital and shell populations;
- degeneracy-aware frontier and signed frontier-bond descriptors;
- Mayer and local π descriptors when the required basis data are defined;
- fixed-geometry GFN1/2 ΔSCF IP, EA, hardness, softness and Fukui descriptors;
- checked RDKit conformer preparation with ETKDG, MMFF94/UFF and optional xTB
  gradient optimization;
- ordered batches, isolated workers, caching and JSONL/CSV reports;
- RDKit YAeHMOP/EHT helpers kept separate from the xTB model.

Density and frontier descriptors require the reviewed extended `libxtb` build
and `keep_overlap=True`. They are extracted from the retained result and do not
invoke a second engine calculation. The package does not provide an excited-state
solver or a solvent model.

### Installation

The complete native profile targets Linux x86-64 with Python 3.11 and a
conda-forge RDKit/Boost toolchain:

```bash
git clone https://github.com/iawnix/rdkit-eth.git
cd rdkit-eth
python tools/install_native.py --jobs 6
conda activate "$PWD/.rdkit-eth/runtime"
python -m rdkit_eth doctor --check compiled --method GFN2-xTB --keep-overlap
```

This builds the reviewed xTB source, its parameters and the optional
`rdXTBTools` module in matched prefixes. Keep `.rdkit-eth/runtime` and
`.rdkit-eth/toolchain` together. See the [installation guide](docs/install.md)
and [offline runtime bundle guide](docs/runtime-bundle.md).

If an extended `libxtb` is already installed beside RDKit, install the Python
adapter and select the library explicitly when needed:

```bash
python -m pip install .
python -m rdkit_eth doctor --check native --method GFN2-xTB
```

The ordinary upstream xTB library does not provide the basis-export extension
used by shell and density descriptors. The selected library path, build hash,
parameter hash and extension capabilities are recorded in every result.

### First calculation

Inputs need explicit hydrogens and a finite 3D conformer. `run_xtb()` leaves the
input molecule unchanged.

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_xtb, calc_xtb_descriptors

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
if AllChem.EmbedMolecule(mol, randomSeed=1) != 0:
    raise RuntimeError("RDKit could not generate a conformer")
AllChem.MMFFOptimizeMolecule(mol)

result = run_xtb(mol, method="GFN2-xTB")
features = calc_xtb_descriptors(result)
print(features["xtb_total_energy_hartree"])
print(features["xtb_homo_ev"], features["xtb_lumo_ev"])
```

### Density and atom-level descriptors

```python
from rdkit_eth import (
    run_xtb, calc_xtb_atom_descriptors, calc_xtb_bond_descriptors,
    calc_shell_populations, calc_density_matrix,
)

result = run_xtb(mol, method="GFN2-xTB", keep_overlap=True)
atoms = calc_xtb_atom_descriptors(result)
bonds = calc_xtb_bond_descriptors(mol, result)
shells = calc_shell_populations(result)
density = calc_density_matrix(result)
```

`keep_overlap=True` retains the native `C`, `S`, AO-atom map and shell map. The
analysis checks `trace(P S)` against the modeled electron count and preserves
`None` plus a reason when a descriptor is undefined. See the
[descriptor catalog](docs/descriptor-catalog.md) and
[density descriptor guide](docs/density-descriptors.md).

### Fixed-geometry reactivity

```python
from rdkit_eth import run_reactivity, XTBOptions

report = run_reactivity(
    mol,
    options=XTBOptions(method="GFN2-xTB", electronic_temperature=300.0),
    population_scheme="lowdin",
    return_errors=True,
)
print(report["values"]["vertical_ip_ev"])
print(report["atoms"][0]["f_plus"])
```

The workflow calculates Q, Q+1 and Q−1 at one fixed geometry. GFN1/2 support
vertical ΔSCF descriptors; GFN0 is not used for this analysis. Electronic
temperature controls orbital occupation smoothing. It is not a nuclear
temperature or an excited-state calculation.

### Reports and batches

```python
from rdkit_eth import run_descriptor_batch, export_descriptor_reports

reports = run_descriptor_batch(molecules, record_ids=ids, max_workers=4,
                               isolate=True)
export_descriptor_reports(reports, "features.jsonl")
```

Reports retain units, analysis settings, engine and parameter identity, quality
flags, missing-value reasons and failures. Scalar columns are stable; atom,
bond, population and density sections remain indexed nested data.

### Geometry preparation

`prepare_conformers()` is an input layer around RDKit embedding and force-field
optimization. `optimize_xtb()` is an opt-in xTB gradient optimizer for workflows
that need xTB-consistent coordinates. Geometry optimization supports descriptor
preparation; it is not a transition-state, solvent or excited-state module.

### Documentation and research

- [Documentation index](docs/README.md)
- [Installation](docs/install.md) · [安装说明](docs/zh-CN/install.md)
- [Descriptor definitions](docs/descriptors.md) · [描述符定义](docs/zh-CN/descriptors.md)
- [Reports and reactivity](docs/descriptor-reports.md) · [报告与反应性](docs/zh-CN/descriptor-reports.md)
- [Descriptor catalog](docs/descriptor-catalog.md) · [描述符目录](docs/zh-CN/descriptor-catalog.md)
- [Density descriptors](docs/density-descriptors.md) · [密度描述符](docs/zh-CN/density-descriptors.md)
- [EAS research plan](docs/research-plan-eas.md) · [EAS 研究方案](docs/zh-CN/research-plan-eas.md)
- [Compiled RDKit interface](docs/rdxtb-tools.md)
- [Extended xTB build](docs/xtb-patched-build.md)
- [Frozen EAS protocol](protocol.json)

### Development

```bash
conda activate "$PWD/.rdkit-eth/runtime"
python -m pip install --no-deps --no-build-isolation -e .
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q
```

Contributions should include a focused test, documented units and invariance,
explicit missing-value semantics, and provenance changes when the public report
schema changes. See [CONTRIBUTING.md](CONTRIBUTING.md) if present in the checkout.

### License

The project code is released under the BSD 3-Clause license. xTB source and
patches retain their upstream licenses; see the bundled license files.

## 中文

### 项目简介

`rdkit-eth` 是一个面向 RDKit 的电子结构描述符库。RDKit 负责分子图、原子
编号和构象；项目构建的扩展版 `libxtb` 负责 GFN0/GFN1/GFN2 气相电子结构
计算。Python/ctypes 接口和可选的 RDKit 风格 C++/Boost.Python 接口使用同一
后端。单位、原子索引、计算来源、引擎哈希和缺失值原因都属于公开接口。

### 已支持功能

- GFN0、GFN1、GFN2 单点计算；
- 分子级、原子级、键级和片段级描述符；
- HOMO/LUMO、能隙、偶极矩、原子电荷和 Wiberg 键级；
- 同一次计算得到的 Mulliken/Löwdin 轨道和 shell 布居；
- 简并不变的前线轨道描述符和前线键分解；
- 满足条件时的 Mayer 键级和局部 π 描述符；
- 固定几何下的 GFN1/2 ΔSCF、IP、EA、硬度、软度和 Fukui 描述符；
- 基于 RDKit 的 ETKDG 构象准备、MMFF94/UFF 优化和可选 xTB 梯度优化；
- 批量计算、隔离 worker、缓存以及 JSONL/CSV 报告；
- 与 xTB 分开的 RDKit YAeHMOP/EHT 辅助接口。

密度和前线描述符需要扩展版 `libxtb` 以及 `keep_overlap=True`。这些描述符
直接消费同一次计算保留的结果，不会再次调用量子引擎。项目当前不提供激发态
求解器或溶剂模型。

### 安装

完整原生环境面向 Linux x86-64、Python 3.11 和 conda-forge RDKit/Boost：

```bash
git clone https://github.com/iawnix/rdkit-eth.git
cd rdkit-eth
python tools/install_native.py --jobs 6
conda activate "$PWD/.rdkit-eth/runtime"
python -m rdkit_eth doctor --check compiled --method GFN2-xTB --keep-overlap
```

已有扩展版 `libxtb` 时，也可以只安装 Python 适配层：

```bash
python -m pip install .
python -m rdkit_eth doctor --check native --method GFN2-xTB
```

普通上游 xTB 库没有本项目密度分析所需的 basis 导出扩展。每个结果都会记录
实际使用的库路径、构建哈希、参数哈希和扩展能力。

### 最小示例

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import run_xtb, calc_xtb_descriptors

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
AllChem.MMFFOptimizeMolecule(mol)

result = run_xtb(mol, method="GFN2-xTB")
features = calc_xtb_descriptors(result)
print(features["xtb_gap_ev"])
```

需要 shell、Mulliken/Löwdin 或密度描述符时：

```python
result = run_xtb(mol, method="GFN2-xTB", keep_overlap=True)
```

### 描述符和研究文档

从[双语文档目录](docs/README.md)开始。核心中文页面包括[安装说明](docs/zh-CN/install.md)、
[描述符定义](docs/zh-CN/descriptors.md)、[报告与反应性](docs/zh-CN/descriptor-reports.md)
和[EAS 研究方案](docs/zh-CN/research-plan-eas.md)。

### 开发与许可

```bash
conda activate "$PWD/.rdkit-eth/runtime"
python -m pip install --no-deps --no-build-isolation -e .
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python -m pytest -q
```

项目代码采用 BSD 3-Clause 许可证；xTB 源码和补丁遵循其上游许可证。
