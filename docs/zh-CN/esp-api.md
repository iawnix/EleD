# 内存中的 xTB 静电势（ESP）接口

两种电荷拟合函数与 `CalcESP` 使用同一套采样点和电势校验。拟合还检查源几何与原子序数是否逐一对应，并在记录了原子数时核对坐标数量。

扩展版 `libxtb` 为非周期 GFN0、GFN1 和 GFN2-xTB 单点计算提供 ESP C ABI。
上游已经把 ESP 积分算法编入共享库，但没有提供对应的公共 C 接口。
原有 Fortran 绘图入口会读写文件并使用固定临时文件名，不适合直接作为外部接口。

EleD 新增可选的结果快照，复用原生 `xtb --esp` 的积分算法和筛选阈值，
在内存中计算用户指定空间点上的静电势，不再启动自洽计算或命令行进程。
EleD 和 rdXTBTools 接口会在单点计算中采样 ESP，返回可序列化的独立数值数组。
ESP/RESP 原子电荷拟合使用这些数组，作为独立的分析步骤。

## Python 使用流程

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import (XTBCalculator, GenerateESPGrid, CalcESP,
                  FitESPCharges, FitRESPCharges)

mol = Chem.AddHs(Chem.MolFromSmiles("O"))
AllChem.EmbedMolecule(mol, randomSeed=42)
points = GenerateESPGrid(mol, scales=(1.4, 1.6, 1.8, 2.0), points_per_atom=100)
result = XTBCalculator(method="GFN2-xTB").calculate(mol, esp_points_angstrom=points)
potential = CalcESP(result)  # 元组，Hartree/e
esp = FitESPCharges(result)
resp = FitRESPCharges(result, equivalence_groups=[(1, 2)])
print(esp.charges_e, resp.charges_e)  # 原始 RDKit 原子顺序，单位 e
record = resp.to_dict()  # 可写入 JSON，包含拟合结果、诊断、设置和来源
```

`esp_points_angstrom` 接受有限实数 `[npoints, 3]` 数组或序列；`None` 或空序列表示
不采样。`XTBOptions` 会保存不可变的独立副本。结果中的 `esp_points_angstrom` 和
`esp_values_au` 分别保存采样坐标与静电势，与其他描述符来自同一次单点计算。
`CalcESP` 读取已有采样值；需要新位置的静电势时，要在计算时提供这些位置。
Python 返回前会释放所有原生句柄，磁盘缓存和工作进程不依赖存活的引擎对象。
采样坐标参与缓存及断点续算配置标识。批处理选项会把同一组坐标用于所有输入；
需要分子表面网格时，应分别为每个分子生成网格。

未采样结果调用 `CalcESP` 会抛出 `DescriptorUnavailableError`；已存数据格式错误
会抛出 `ValueError`。计算采样值需要通过校验的 ESP 扩展库，不会用原子点电荷势代替。

编译接口的等价用法：

```python
from eled import rdXTBTools
ok, compiled = rdXTBTools.RunMol(mol, method="GFN2-xTB", espPointsAngstrom=points)
if not ok:
    raise RuntimeError(compiled.errorMessage)
xyz = compiled.GetESPPoints()  # 独立 NumPy 副本，Angstrom [point, xyz]
v = compiled.GetESPValues()   # 独立 NumPy 副本，Hartree/e [point]
resp = FitRESPCharges(compiled, equivalence_groups=[(1, 2)])
```

`compiled.espSampled` 表示是否有采样结果。公共的 `CalcESP`、`FitESPCharges` 和
`FitRESPCharges` 接受两种后端的结果。独立 C++ 接口使用
`XTBOptions::espPointsAngstrom` 和 `XTBResults::espPointsAngstrom/espValuesAU`
（按点展开的向量）。已有 Mulliken/EEQ/CM5 字段保持其定义；拟合电荷通过独立的
`ESPChargeFit` 返回，并标明 `charge_model`。

## 采样和拟合定义

`GenerateESPGrid(mol, conf_id=-1, scales=(1.4,1.6,1.8,2.0), points_per_atom=100)`
使用 RDKit 范德华半径，以及固定在实验室坐标系中的 Fibonacci 球面方向。
每个原子的每层球面生成候选点，移除落入同层其他原子球体内部的点。
输入必须有显式氢原子和有限三维坐标。这是一套明确的采样规则，不称为 CHELPG
或 Merz–Kollman。应保存采样坐标及网格设置。分子旋转后重新生成的网格会改变；
验证刚体不变性时，应让同一组采样点随分子一起变换。

拟合矩阵为 `A[k,i] = 1 / distance_bohr(point_k, atom_i)`，目标值是计算得到的
分子静电势，总电荷来自计算来源信息，拟合电荷单位为 e。接口签名为：

```python
FitESPCharges(result, *, weights=None, equivalence_groups=(), fixed_charges=None,
              rcond=1e-12)
FitRESPCharges(result, *, restraint_strength=0.0005, restraint_width=0.1,
               restrain_hydrogens=False, weights=None, equivalence_groups=(),
               fixed_charges=None, max_iterations=1000, tolerance=1e-8, rcond=1e-12)
```

ESP 最小化 `0.5 * sum(w * (Aq - V)**2)`。单阶段 RESP 另外加入
`sum(a_i * (sqrt(q_i**2 + b**2) - b))`，其中 `b` 为 `restraint_width`（单位 e），
`a_i` 为 `restraint_strength`，默认不约束氢原子，所以氢的 `a_i=0`。
强度参数是上述原子单位数值目标函数中的系数。权重乘在残差平方上，默认均为 1，
允许零，必须非负且至少有一个正值。权重**不会归一化**；改变权重尺度或采样密度
会改变惩罚项的相对强度。

总电荷、用户指定的等价关系 `equivalence_groups` 和固定电荷
`fixed_charges={原子编号: 电荷_e}` 通过约束零空间求解精确满足。编号遵循输入顺序，
不会自动根据分子图推断等价原子。约束冲突时失败。ESP 要求剩余自由参数方向上的
拟合矩阵满秩；RESP 只在惩罚项消除剩余歧义时才能接受欠秩数据，否则也会失败。
`rcond` 控制秩判定。RESP 在最大电荷变化不超过 `tolerance` e 时收敛，达到迭代
上限仍未收敛时抛出 `RuntimeError`，不返回成功拟合结果。

`ESPChargeFit` 包含原子电荷、拟合静电势、残差（拟合减参考）、加权 RMS、最大残差、
约束误差、数据秩、自由参数数量、条件数、迭代次数和目标函数值。
数据矩阵欠秩时 `condition_number=None`，即使 RESP 惩罚项能够确定唯一解也如此。
结果还包含完整拟合设置、引擎/参数/计算来源，以及数值输入与拟合设置的 SHA256。
协议 ID 为 `eled-esp-v1` 和 `eled-resp-single-stage-v1`。

此 RESP 接口只实现一个双曲惩罚阶段，不宣称实现完整 Amber 两阶段协议、规定的
QM 层级/采样网格或经过验证的力场电荷。可以分别分析多个构象；多构象联合拟合
和自动化学等价识别不在此接口范围内。基组密度、几何、采样和惩罚都会影响电荷；
解释结果前应检查残差及条件数。

## C 接口

```c
int xtb_rdkit_getESPAPIVersion(void);  /* 返回 1 */
void xtb_rdkit_setKeepESP(xtb_TEnvironment env, xtb_TResults res, bool keep);
int xtb_rdkit_getESP(xtb_TEnvironment env, xtb_TResults res,
                    int npoints, const double *points,
                    int capacity, double *output);
```

在创建环境、分子、已加载参数的计算器和结果对象后：

```c
xtb_rdkit_setKeepESP(env, res, true);
xtb_singlepoint(env, mol, calc, res);
/* 读取前先检查 xtb_checkEnvironment(env)。 */
double points[2][3] = {{4.1, 2.3, 3.7}, {-3.2, 4.5, 2.1}};
double potential[2];
int ok = xtb_rdkit_getESP(env, res, 2, &points[0][0], 2, potential);
/* 只有 ok == 1 时才能读取 potential；否则查看 xtb_getError。 */
```

采样点为连续存储的 `double[npoints][3]`，单位 **Bohr**，使用计算时的坐标系。
输出为 `double[npoints]`，单位为**静电势原子单位 Hartree/e**；正负号对应正单位
试探电荷的势能符号。调用方负责分配两个缓冲区，且必须满足 `capacity == npoints > 0`。
输出顺序与输入采样点严格一致。

静电势由 xTB 价电子基组中的电子密度和分子有效核芯贡献组成。即使计算时溶剂或
外部点电荷极化了波函数，输出也不包含溶剂反应场或外部点电荷的直接贡献。
它不是核附近的全电子静电势，也不是周期体系静电势。

成功返回 1。缺少快照、不支持的计算、空指针、错误的点数或容量、非有限坐标、
距原子核小于 1e-10 Bohr 的采样点，以及计算得到非有限值时返回 0，错误原因通过
`xtb_getError` 读取。失败时不改写输出缓冲区，不用零表示缺失值。
环境指针为空时返回 0，无法提供环境诊断。

## 生命周期与数值范围

默认不保留快照，必须在 `xtb_singlepoint` **之前**开启；计算后开启不会补建旧结果。
快照独立保存几何、基组、占据数和笛卡尔轨道系数，不依赖分子与计算器对象的寿命。
`xtb_copyResults` 会深复制快照。每次调用保留开关、新的单点尝试（包括初始化失败
和 SCF 失败）、Hessian 或 CPCM-X 接口，都会使旧快照失效。保留开关仍开启时，
后续受支持的成功单点会重新生成快照。额外存储量为 O(NAO²)。
GFN-FF 和周期体系不生成 ESP 快照。

筛选阈值与 CLI 默认值一致：占据数 1e-6、打包密度矩阵 1e-4、原始高斯函数对
指数截断 10、密度与收缩系数乘积 1e-5、加权原始函数对 5e-7。求和顺序与 CLI 不同。
笛卡尔积分筛选会限制严格的旋转不变性，这是数值近似，不是改变电荷划分方式。
此接口不保证整个 libxtb 的线程安全；调用方仍应像 EleD 一样串行执行引擎操作。

数值验证覆盖分子表面尺度的采样点。在极远位置（例如 10,000 Bohr），上游积分
例程可能打印收敛警告，同时仍返回有限数值。本扩展保留该行为；此时返回有限值
本身不能证明结果精度。

## 构建和验证

复用标准安装中的固定运行时和工具链：

```bash
conda run --no-capture-output -p "$PWD/.eled/toolchain" \
  python tools/build_patched_xtb.py "$PWD/.eled/runtime" \
  --dependency-prefix "$PWD/.eled/toolchain" --jobs 6
.eled/runtime/bin/python tools/diagnostics/verify_esp.py .eled/runtime
```

构建器检查 ESP 补丁与修改后源码的哈希，运行上游 CTest 和已有基组/SCC 验证，
再比较安装后的 ESP ABI 与 CLI 数值。通过后才在构建记录中写入
`capabilities.esp_api = 1` 和 `esp_checks`。已加载旧动态库的进程需要重启。

验证覆盖三种方法下的水、硫化氢（包含 d 轨道笛卡尔变换）、铵离子和羟基自由基，
以及结果复制、关闭保留、错误输入、初始化/SCF 失败、无 ESP 文件输出和水分子的
刚体变换、原子重排、采样点重排。`tests/test_esp.py` 也运行同一套集成验证。
标准构建未启用 CPCM-X。
