# 安装 EleD

EleD 的安装器统一配置 RDKit、EleD、项目构建的 xTB 计算库（`libxtb`）和
`rdXTBTools` C++ 扩展，使用包含基组导出和 SCC 诊断扩展的 xTB 版本。

## 完整安装

当前安装流程支持 Linux x86-64，安装前请准备好 conda。安装器会配置 Python 3.11
和匹配的 RDKit/Boost：

```bash
git clone https://github.com/iawnix/EleD.git
cd EleD
python tools/install_native.py --jobs 6
conda activate "$PWD/.eled/runtime"
python -m eled doctor --check compiled --method GFN2-xTB --keep-overlap
```

安装器会创建两个相互匹配的目录：

- `.eled/runtime`：Python、RDKit、项目构建的 `libxtb`、参数文件和
  `rdXTBTools`；
- `.eled/toolchain`：编译器和原生构建依赖。

两个目录需要一起保留在原始路径。安装完成后，`doctor` 会检查 GFN0/1/2、
overlap/basis 导出、描述符和固定几何反应性流程。

## 最小计算

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import XTBCalculator, CalcXTBDescriptors

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
AllChem.MMFFOptimizeMolecule(mol)

result = XTBCalculator(method="GFN2-xTB").calculate(mol)
print(CalcXTBDescriptors(result))
```

输入需要显式氢、有限的 3D 构象和 RDKit 原子顺序。`XTBCalculator` 不修改输入
分子。

## 构建 C++ 接口

```bash
ELED_BUILD_RDXTB=1 \
RDXTB_XTB_PREFIX=/absolute/path/to/extended-xtb \
python -m pip install --no-deps --no-build-isolation .
```

该模块需要与当前 RDKit/Boost 环境匹配。完整安装、离线 runtime bundle 和
故障排查见英文 [install.md](../install.md) 与 [runtime-bundle.md](../runtime-bundle.md)。

在已安装环境中修改 Python 代码，见[开发指南](development.md)。
