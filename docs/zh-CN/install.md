# 安装 rdkit-eth

## 完整原生环境

推荐使用 Linux x86-64、Python 3.11 和 conda-forge RDKit/Boost：

```bash
git clone https://github.com/iawnix/rdkit-eth.git
cd rdkit-eth
python tools/install_native.py --jobs 6
conda activate "$PWD/.rdkit-eth/runtime"
python -m rdkit_eth doctor --check compiled --method GFN2-xTB --keep-overlap
```

安装器会创建两个相互匹配的目录：

- `.rdkit-eth/runtime`：Python、RDKit、扩展版 `libxtb`、参数文件和
  `rdXTBTools`；
- `.rdkit-eth/toolchain`：编译器和原生构建依赖。

两个目录需要一起保留在原始路径。安装完成后，`doctor` 会检查 GFN0/1/2、
overlap/basis 导出、描述符和固定几何反应性流程。

## 已有扩展版 libxtb

如果环境中已经有项目要求的扩展版 `libxtb`，可以只安装 Python 适配层：

```bash
python -m pip install .
python -m rdkit_eth doctor --check native --method GFN2-xTB
```

也可以显式指定库：

```python
from rdkit_eth import XTBCalculator

result = XTBCalculator(xtb_library="/absolute/path/to/libxtb.so").calculate(mol)
```

上游未修改的 xTB 库不包含本项目的 basis 导出扩展，因此不能支持
Mulliken/Löwdin shell 描述符。程序会检查扩展能力并记录库路径和 SHA-256。

## 最小计算

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import XTBCalculator, calc_xtb_descriptors

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
AllChem.MMFFOptimizeMolecule(mol)

result = XTBCalculator(method="GFN2-xTB").calculate(mol)
print(calc_xtb_descriptors(result))
```

输入需要显式氢、有限的 3D 构象和 RDKit 原子顺序。`XTBCalculator` 不修改输入
分子。

## 构建 C++ 接口

```bash
RDKIT_ETH_BUILD_RDXTB=1 \
RDXTB_XTB_PREFIX=/absolute/path/to/extended-xtb \
python -m pip install --no-deps --no-build-isolation .
```

该模块需要与当前 RDKit/Boost 环境匹配。完整安装、离线 runtime bundle 和
故障排查见英文 [install.md](../install.md) 与 [runtime-bundle.md](../runtime-bundle.md)。
