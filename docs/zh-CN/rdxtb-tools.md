# 编译接口 rdXTBTools

`eled.rdXTBTools` 是一个可选的 Boost.Python 扩展，底层是 C++ 库。它接受 RDKit 原生的
`ROMol`，直接调用 xTB C API，并仿照 `rdEHTTools.RunMol` 的风格返回
`(success, XTBResults)`。C++ 命名空间为 `RDKit::XTBTools`。它与 YAeHMOP 无关，也不属于
上游 RDKit 发行版。

推荐的完整安装方式是先运行 `python tools/install_native.py`，再执行
`conda activate "$PWD/.eled/runtime"`。它会把扩展版 libxtb 安装在同一前缀中；下文的手动
构建方式仍可用于开发。

## 构建与安装

支持的开发平台是 Linux x86-64、Python 3.11，以及配套的 conda-forge RDKit 2026.09.1 /
Boost 1.90 SDK。完整的锁定安装见 [install.md](install.md)。SDK 和 xTB 源码工具链使用
不同的前缀。先按该文档构建或选定 xTB，然后在项目根目录执行：

```bash
conda create --yes --prefix "$PWD/.rdxtb-build-env" --file locks/linux-64/rdxtb-sdk.txt
conda activate "$PWD/.rdxtb-build-env"
ELED_BUILD_RDXTB=1 RDXTB_XTB_PREFIX="$PWD/.xtb-native" python -m pip install --no-deps --no-build-isolation -e .
```

该环境会安装 RDKit 及其 C++ 开发头文件/库、Boost、NumPy、OpenSSL、C++ 编译器、CMake 和
Ninja。xTB 从它自己的前缀中选取。不设置 `ELED_BUILD_RDXTB=1` 的普通 pip 安装仍然只装纯
Python 适配器。仅有 pip 安装的 RDKit wheel 并不提供这里需要的 C++ SDK。构建和导入扩展请
使用同一个环境；即使版本号相同，混用不同构建的 RDKit 或 Boost 也可能破坏 C++ ABI。导入时
会检查 RDKit 版本，但无法保证所有 ABI 组合都兼容。

请始终通过 `RDXTB_XTB_PREFIX` 指定要使用的 xTB 安装，例如
`RDXTB_XTB_PREFIX=/prefix/with/xtb`；该前缀必须提供 `include/xtb.h` 和 `lib/libxtb.so`。
`RDXTB_RDKIT_PREFIX` 可覆盖 SDK 前缀。运行时选择库的顺序是：`xtbLibrary` 参数、
`XTB_LIBRARY`、构建时记录的 xTB 路径；如果构建时未配置路径，再交给系统加载器查找
`libxtb.so`。显式指定的路径加载失败时直接报错。该模块支持 xTB C API 主版本 1 和 2 共有
的调用。这个 ABI 版本不是 xTB 程序版本。

启动 Python 前先选定扩展安装：

```bash
source .xtb-native/activate.sh
```

C++ 核心会在**加载任何参数之前**校验完整的构建清单、审查过的上游修复身份、基组导出和
SCC 诊断等扩展，以及所选二进制的摘要。对于 GFN0/1/2，未打补丁的官方库会被拒绝。
`rdXTBTools` 始终使用扩展版原生库。另一个 `XTBCalculator` 接口默认通过 ctypes 调用
libxtb，并执行相同的校验；两个接口都不会回退到其他引擎或进程。

在同一个已激活的环境中构建 wheel：

```bash
ELED_BUILD_RDXTB=1 python -m pip wheel --no-deps --no-build-isolation -w dist .
```

这些是与 SDK 和 xTB 安装绑定的本地开发二进制，不打包 RDKit、Boost、libxtb 及其
Fortran/数值依赖。请保留这些前缀。该 wheel 不是可移植的 manylinux 发行包。macOS 和
Windows 需要平台加载器支持，目前尚未实现。

## Python API 与单位

```python
from rdkit import Chem
from rdkit.Chem import AllChem, rdEHTTools
from eled import rdXTBTools

mol = Chem.AddHs(Chem.MolFromSmiles("CCO"))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
ok, xtb = rdXTBTools.RunMol(mol, method="GFN2-xTB")
if not ok:
    raise RuntimeError(xtb.errorMessage)
print(xtb.totalEnergyHartree)
print(xtb.GetOrbitalEnergies())
print(xtb.GetWibergBondOrders())
print(xtb.chargeModel)  # GetAtomicCharges() 使用的电荷模型

eht_ok, eht = rdEHTTools.RunMol(mol)
if eht_ok:
    print(eht.totalEnergy)  # eV，EHT 模型自身的能量
```

| 物理量 | `rdXTBTools.XTBResults` | 与 EHT 的区别 |
| --- | --- | --- |
| 总能量 | `totalEnergyHartree`，Hartree | EHT 的 `totalEnergy` 单位是 eV，定义和能量零点也不同 |
| 轨道能量 | `GetOrbitalEnergies()`，eV，`[MO]` | 都用 eV，但哈密顿量和基组不同 |
| 占据数 | `GetOrbitalOccupancies()`，电子，`[MO]` | 有限电子温度下 xTB 的占据数可以是分数 |
| 原子电荷 | `GetAtomicCharges()`，元电荷，`[atom]` | **GFN1/2：** 收敛的 SCC Mulliken 电荷（`chargeModel == "xtb_scc_mulliken"`）；**GFN0：** PEEQ/EEQ（`"gfn0_eeq"`） |
| CM5 电荷 | `GetCM5Charges()`，元电荷，`[atom]` | 对 GFN1/2 SCC Mulliken 电荷的几何相关修正；需要 GFN1/2 SCC 电荷和计算几何 |
| 采样 ESP | `GetESPPoints()`（Å，`[point,xyz]`）、`GetESPValues()`（Hartree/e，`[point]`） | 向 `RunMol` 传入 `espPointsAngstrom`；`espSampled` 表示是否可用。共享的 [ESP/RESP 拟合](esp-api.md) 可接受此结果 |
| 采样实空间性质 | `realSpaceExportVerified`、`realSpaceSampled`、`GetRealSpacePoints()`（Å，`[point,xyz]`）、`GetRealSpaceDensity()`（`e/Bohr^3`）、`GetRealSpaceALIE()`（`eV`）、`GetRealSpaceVirtualDensity()`（`1/Bohr^3`）、`GetRealSpaceLEA()`（`eV`） | 向 `RunMol` 传入 `realSpacePointsAngstrom`、`realSpaceOccupationThreshold`、`realSpaceLEAWindowEV`；共享的 `CalcRealSpaceProperties` 和 `CalcRealSpaceSurfaceDescriptors` 可接受此结果 |
| 偶极矩 | `GetDipole()`，e·Bohr，`[x,y,z]` | 带电分子的偶极依赖坐标原点；使用传入的坐标 |
| 核梯度 | `GetGradient()`，Hartree/Bohr，`[atom,xyz]` | 供 `OptimizeXTB` 使用；不是激发态梯度 |
| 键级 | `GetWibergBondOrders()`，无量纲，`[atom,atom]` | 不是 EHT 的重叠布居矩阵 |
| MO 系数 | `GetOrbitalCoefficients()`，`[AO,MO]` | 可选；按 xTB 基组顺序，不是 EHT 基组顺序 |
| 重叠矩阵/映射 | `GetOverlapMatrix()`、`GetAOAtomMap()`、`GetAOShellMap()`、`GetShellAtomMap()`、`GetShellAngularMomenta()` | 可选的原生基组数据；不是 EHT 矩阵 |
| 能量分解与费米能级 | `energyDecompositionAvailable`、`energyDecompositionFlags`（引擎给出的位掩码）、`electronicEnergyHartree`、`repulsionEnergyHartree`、`electrostaticEnergyHartree`、`anisotropicElectrostaticEnergyHartree`、`anisotropicXCEnergyHartree`、`dispersionEnergyHartree`、`halogenBondEnergyHartree`、`solvationEnergyHartree`（及 Born/SASA/HB/shift 分项）、`fermiLevelAlphaHartree`/`EV`、`fermiLevelBetaHartree`/`EV` | 方法不适用、气相或 `electronicTemperature <= 0.1 K` 时为 `None` |
| D4 极化率与色散系数 | `polarizabilityAvailable`、`molecularPolarizabilityAU`（Bohr³）、`molecularC6AU`（Ha·Bohr⁶）、`molecularC8AU`（Ha·Bohr⁸）、`GetAtomicPolarizabilities()`、`GetAtomicC6()` | 仅 GFN2 和 GFN0；GFN1（DFT-D3）不可用 |
| CAMM 多极矩 | `multipolesAvailable`、`GetAtomicDipoles()`（e·Bohr，`[atom,3]`）、`GetAtomicQuadrupoles()`（e·Bohr²，`[atom,6]`）、`GetMolecularQuadrupole()`（e·Bohr²，`[6]`） | 无迹，分量顺序 `(xx,xy,yy,xz,yz,zz)`；仅 GFN2/GFN1（GFN0 不可用） |
| 原子核处静电势（EPN） | `nuclearESPAvailable`、`GetNuclearESP()`（Hartree/e，`[atom]`） | 价电子与其他原子核（芯电子屏蔽后的电荷 $Z-n_\mathrm{core}$）产生的精确电势，不含自身核项 $Z_A/r_{iA}$ |
| 隐式溶剂化溯源 | `solvent`、`solvationModel`、`solvationGrid` | 气相为 `None`；设置 `solvent` 时为 `"alpb"` 或 `"gbsa"` |
| 模型电子数 | `numElectrons`，占据数之和取整 | 价电子/模型电子；水为 8 |
| 总电子数 | `numTotalElectrons`，核电荷之和减去分子电荷 | 包含芯电子；水为 10 |

这里没有与 EHT 兼容的 `totalEnergy`、`fermiEnergy`、哈密顿量或重叠布居别名。把 Hartree
换算成 eV 并不能让两种总能量变得可比。轨道能隙不是激发能。在经过验证的原生扩展下，
`keepOverlap=True` 会保留同一次计算的 S 和 C，以及从零开始编号的 AO/壳层/原子映射，从而
可以使用共享的布居分析函数。轨道符号和简并轨道之间的旋转仍是任意的。

调用签名：

```python
RunMol(mol, confId=-1, method="GFN2-xTB", charge=None,
       numUnpairedElectrons=None, keepOrbitalCoefficients=False,
       xtbLibrary="", parameterDirectory="", maxIterations=250,
       accuracy=1.0, electronicTemperature=300.0, keepOverlap=False,
       espPointsAngstrom=None, solvent=None, solvationModel="alpb",
       solvationGrid=None, realSpacePointsAngstrom=None,
       realSpaceOccupationThreshold=1e-7, realSpaceLEAWindowEV=0.0)
```

`confId` 是 RDKit 构象的 **ID**，不是它在构象列表中的位置；`-1` 表示第一个构象。输入坐标
单位为 Å，必须有限且标记为 3D。请用 `Chem.AddHs` 把所有氢计数展开成原子。电荷和未成对
电子数默认取 RDKit 的形式电荷和自由基计数，并校验电子数与奇偶性。电子温度单位为开尔文；
accuracy 是 xTB 无量纲的数值精度设置。`solvent` 默认为 `None`（气相）；传入非空的溶剂名
会通过 `solvationModel`（`"alpb"` 或 `"gbsa"`）和可选的 Lebedev 网格 `solvationGrid`
启用 xTB 隐式溶剂化（`None` 使用 xTB 默认的 230 点网格；独立 C++ 的
`RDKit::XTBTools::XTBOptions` 用 `""` 和 `0` 表示）。C++ 模块只做单点计算，不内置优化
方法；Python 包通过 `OptimizeXTB` 提供 xTB 梯度优化，通过
`CalcVibrationalThermochemistry` 提供谐振热化学，`rdXTBTools` 只提供底层的单点数据。

参数文件的查找顺序：先 `parameterDirectory`，再以冒号分隔的 `XTBPATH`，最后是所选库的
`share/xtb`。显式指定的目录必须包含所需的参数文件。找不到文件时，GFN1/2 可以使用编译时
内置的参数；GFN0 必须有参数文件。元数据记录实际使用的库路径/哈希、参数文件路径/哈希
（内置参数时为空）、C API 版本、所选方法/构象、电荷/自旋，以及可用时的 GFN0 校验信息和
源码提交。哈希只标识文件本身，不覆盖所有传递依赖。替换已加载的库之后请重启进程。

无效的输入或选项抛出 `ValueError`；库和 ABI 错误抛出 `RuntimeError`。xTB 初始化、SCF 或
结果提取出错时返回 `False`，错误信息在 `result.errorMessage` 中。对失败结果调用数值
getter 会抛出 `RuntimeError`。可选的 MO 系数未保留时会抛出 `ValueError`，提示用
`keepOrbitalCoefficients=True` 重新计算。

结果对象拥有自己的存储。每个 NumPy getter 都返回副本；标量元数据在 Python 中只读。原子
顺序与输入分子一致，输入分子不会被修改。包装层在持有 GIL 时对 `ROMol` 做快照，计算期间
释放 GIL。C++ 调用是串行化的；Python 调用还与 ctypes 适配器共用同一把锁。独立的其他 xTB
绑定以及单独链接的静态 C++ 库副本不共享这把锁。这并不说明 libxtb 在一般意义上是线程
安全的。

底层的 `rdXTBTools.XTBResults` 对象与 Python 的 `XTBResults` 数据类具有对应的
`camelCase` 标量属性和返回 NumPy 副本的 `Get*()` 方法，都可以用于共享的描述符函数和
`AnnotateMol`。两者经由同一个归一化函数转换为相同的字段和溯源信息，因此标量描述符表在
两种接口下完全一致。批处理/缓存 API 通过 ctypes 生成 `XTBResults` 数据类。独立 C++ 的
`RDKit::XTBTools::XTBResults` 也会通过 `RDKit::XTBTools::calcCM5Charges` 直接为
`GFN1/2` 填充 `cm5Charges`。

## 共享的描述符提取

```python
from eled import CalcXTBDescriptors, CalcOrbitalPopulations, AnnotateMol
ok, xtb = rdXTBTools.RunMol(mol, keepOverlap=True)
if not ok:
    raise RuntimeError(xtb.errorMessage)
features = CalcXTBDescriptors(xtb)
populations = CalcOrbitalPopulations(xtb)
annotated = AnnotateMol(mol, xtb)
```

这些函数从不重新计算。包装层记录 C++ 实际使用的分子快照的溯源信息；之后若原子顺序、
拓扑或几何发生变化，映射会被拒绝。独立 C++ 结果的 `inputProvenanceJSON` 为空。
`basisExportVerified` 表示经过验证的 ABI/清单能力；`overlapRetained` 表示这个结果是否
真的持有基组数据。即使只做基本单点计算，包装层也会拒绝缺少完整审查扩展集的上游库。
定义见[描述符说明](descriptors.md)。

## 独立使用 C++

C++ 核心不需要 Python 运行时。单独配置并安装：

```bash
cmake -S . -B build/core -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="$CONDA_PREFIX" \
  -DRDXTB_XTB_PREFIX="$PWD/.xtb-native" \
  -DRDXTB_BUILD_PYTHON=OFF \
  -DCMAKE_INSTALL_PREFIX="$PWD/build/rdxtb-sdk"
cmake --build build/core --parallel 4
ctest --test-dir build/core --output-on-failure
cmake --install build/core
```

只有当 `$CONDA_PREFIX` 中包含本项目构建的扩展库时，才使用
`-DRDXTB_XTB_PREFIX="$CONDA_PREFIX"`；上游发布版 xTB 不满足运行时约定。外部 CMake 项目
可以这样使用：

```cmake
find_package(RDKitXTB CONFIG REQUIRED)
target_link_libraries(my_program PRIVATE RDKitXTB::RDKitXTB)
```

配置该项目时，把 SDK 安装前缀和 conda 环境都加入 `CMAKE_PREFIX_PATH`。其 C++ 代码可以
这样调用：

```cpp
#include <XTBTools.h>

RDKit::XTBTools::XTBOptions options;
options.method = "GFN2-xTB";
RDKit::XTBTools::XTBResults result;
if (!RDKit::XTBTools::runMol(mol, result, options)) {
  throw std::runtime_error(result.errorMessage);
}
// result.totalEnergyHartree; result.orbitalEnergiesEV
```

无效输入抛出 `std::invalid_argument`；库的问题抛出 `std::runtime_error`。读取 C++ 数值
字段之前请先检查返回的布尔值；与 Python getter 不同，公开的 C++ 字段不会对失败结果做
保护。矩阵是自有的行主序 vector，各轴含义见 `XTBTools.h`。运行时通过 `dlopen` 加载，
链接时不依赖 libxtb。

## 验证

```bash
source .xtb-native/activate.sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 ELED_REQUIRE_RDXTB=1 \
  python -m pytest -q
python examples/rdxtb_tools.py
```

设置这个"必须存在模块"的标志后，缺少扩展不会让编译接口测试被悄悄跳过；扩展存在但导入
失败时，测试收集总会失败。测试在乙醇、水、苯、铵离子和甲基自由基上比较三种方法；验证
单位、系数布局、数组所有权、输入校验、SCF 失败与恢复、与 EHT 共存以及混合的原生线程；
并在初始化之前拒绝未经验证的 GFN0。C++ 测试还在不依赖 Python 的情况下执行多次调用。
[编译接口 CI 工作流](../../.github/workflows/rdxtb-tools.yml) 覆盖发布版和打补丁的
依赖。见[本地验证记录](../validation.md)（英文）。

## 共享的数值设置（0.6）

`maxIterations`、`accuracy` 和 `electronicTemperature` 在 ctypes 中有对应的、经过校验的
选项，并与 `effectiveElectronicTemperature` 一起保留在编译结果中。升级后请重新编译 C++
使用方。取值范围、单位、批量导出以及可接受编译结果的 GFN1/2 电荷态流程见
[描述符报告](descriptor-reports.md)。

## 电荷态重启协议（0.7）

独立 C++ 的 `runChargeStates(mol, options, ionizedUnpaired, attachedUnpaired,
reuseReference=true)` 依次返回参考态、电离态和加电子态结果，在第一次失败处停止（失败的
结果仍保留在 vector 中）。只支持 GFN1/2；相邻电荷态的自旋数需显式给出。两个相邻电荷态
都从未改变的参考态深拷贝初始化。引擎句柄的作用域限于这次调用，调用期间持有计算锁。
`XTBResults.initialization` 记录 `independent` 或 `reference`。升级后请重新编译使用方。

Python 用户调用 `RunReactivity(..., interface="compiled")`；其内部的编译接口桥接在释放
GIL 的同时也持有共享的 ctypes 锁。
