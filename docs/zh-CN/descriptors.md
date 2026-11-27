# 描述符（1.4）

公开函数/方法采用 RDKit `PascalCase`，参数与结果属性采用 `camelCase`，导出的描述符字段名统一采用带显式单位后缀的小写 `snake_case`。基于结果对象的描述符函数只分析已有的 `XTBResults`；需要一步从 `Mol` 提取分子标量或接入 `rdMolDescriptors.Properties` 时，可使用 `CalcMolXTBDescriptors()` 与 `RegisterRDKitProperties()`；只有 `XTBCalculator.Calculate()`、`CalcMolXTBDescriptors()`、`RunReactivity()` 等入口会调用引擎。

## 单点描述符

```python
from eled import (
    XTBCalculator, CalcMolXTBDescriptors, CalcXTBDescriptors,
    CalcXTBAtomDescriptors, CalcXTBBondDescriptors, AnnotateMol,
)

result = XTBCalculator(method="GFN2-xTB").Calculate(mol)
molecular = CalcXTBDescriptors(result)
atoms = CalcXTBAtomDescriptors(result)
bonds = CalcXTBBondDescriptors(mol, result)
annotated = AnnotateMol(mol, result, privateProps=True)  # 写入 _ 前缀的 RDKit 私有属性
# 或单步提取标量字典：CalcMolXTBDescriptors(mol)
```

如需将数值型分子描述符注册到 RDKit 的 `rdMolDescriptors.Properties`（同一构象的多次属性查询共享一次单点计算，缺失标量返回 `NaN`）：

```python
from rdkit.Chem import rdMolDescriptors
from eled import RegisterRDKitProperties

names = RegisterRDKitProperties(method="GFN2-xTB")
props = rdMolDescriptors.Properties(["xtb_energy_ha", "xtb_gap_ev", "xtb_dip_debye"])
energy_ha, gap_ev, dipole_debye = props.ComputeProperties(mol)
```

主要字段包括：

- `xtb_energy_ha`：总电子能，Hartree；
- `xtb_homo_ev`、`xtb_lumo_ev`、`xtb_gap_ev`：eV；
- `xtb_dip_*_ebohr`、`xtb_dip_debye`：偶极分量和模长；
- `xtb_q_*_e`：xTB 原子电荷统计；GFN1/2 是 SCC Mulliken，GFN0 是 PEEQ/EEQ，准确模型写入 provenance；
- `xtb_cm5_q_*_e`、原子 `cm5_q_e`：由 GFN1/2 SCC Mulliken 电荷和计算几何显式推导的 CM5 电荷；
- `q_e`、`wbo_sum`：原子级特征；
- `wbo`：对应 RDKit 键的 Wiberg 键级。

CM5 在 SCC Mulliken 电荷上加入与距离有关的反对称原子对校正，单位为元电荷。
H/C/N/O 元素对使用 xTB 的特定参数，其他元素对使用按原子顺序计算的原子参数差。
校正保持总电荷守恒；原子重新编号后，各原子的电荷随映射保持一致。坐标重合时
抛出 `ValueError`。本次修复纠正了非特定元素对的额外符号反转；涉及这些元素对的
历史 CM5 归档应使用保存的几何和 SCC 电荷重新生成。SCC 电荷和 Wiberg 键级
不受此校正影响。

CM5 缺失时，报告为 GFN0 记录 `cm5_unsupported_method`，为 GFN1/2 记录
`cm5_not_reported`。每组已提供的 CM5 电荷必须与原子逐一对应。

电荷数组沿用计算时的原子顺序。直接提取、Fukui 分析和报告都会检查电荷模型及原子数。缺少重叠矩阵、轨道系数或 AO 到原子的映射时，报告保留标量结果并记录轨道部分缺失。已提供的矩阵维度或原子索引有误时，接口报告输入错误。

HOMO/LUMO 由轨道占据阈值定义。能隙是基态前线轨道差值，不是光学激发能。

## 密度和 shell 描述符

```python
from eled import XTBCalculator, CalcShellPopulations, CalcDensityMatrix

result = XTBCalculator(method="GFN2-xTB").Calculate(mol, keepOverlap=True)
shells = CalcShellPopulations(result)
density = CalcDensityMatrix(result)
```

扩展版 `libxtb` 会返回同一次计算中的 `C`、`S`、AO atom map 和 shell map。
程序构造

```text
P = C · diag(occupation) · Cᵀ
```

并检查 `trace(P S)` 是否等于电子数。可用字段包括：

- `mul_s_e`、`mul_p_e`、`mul_d_e`；
- `low_s_e`、`low_p_e`、`low_d_e`；
- `mul_ps`、`low_ps`；
- `tr_ps_e`、`n_err_e`。

这些是基组相关的密度分区量，不是唯一的实验可观测量。GFN0 的 EEQ 电荷
与 Mulliken/Löwdin 布居不是同一个定义。

## 前线与 π 描述符

```python
from eled import (
    CalcFrontierDescriptors,
    CalcFrontierBondDescriptors, CalcPiDescriptors,
)

frontier = CalcFrontierDescriptors(mol, result)
frontier_bonds = CalcFrontierBondDescriptors(mol, result)
pi = CalcPiDescriptors(mol, result)
```

这些接口提供简并轨道组布居、participation ratio、熵、带符号的 HOMO/LUMO
键重叠、π 电子数和 π 键级。局部 π 特征需要定义明确的局部平面法向。条件不满足时返回 `None` 和原因。

## 反应性描述符

```python
from eled import RunReactivity

report = RunReactivity(mol, populationScheme="lowdin", returnErrors=True)
```

固定几何计算 Q、Q+1、Q−1，得到 `ip_v_ev`、`ea_v_ev`、
`mu_ev`、`chi_ev`、`eta_ev`、`s_1ev`、`omega_ev`、Gázquez 给/受电子能力
`omega_minus_ev`/`omega_plus_ev`，以及原子级
`f_plus`、`f_minus`、`f_radical`、`dual`、`s_loc_1ev`、
`omega_loc_ev`、`n_loc_1ev`、`omega_minus_loc_ev`、`omega_plus_loc_ev`、
`delta_omega_loc_ev` 和局域超软度 `s2_loc_1ev2`。该流程支持 GFN1/2，
不自动搜索自旋态，也不进行电荷态几何优化。

## Grimme FOD 静态关联、凝聚线性响应矩阵与键中点局域电场

```python
from eled import (
    RunFOD, CalcFODFromResult, CalcLinearResponseMatrix,
    CalcBondElectricField, CalcDensityIsosurface,
    CalcRealSpaceProperties, CalcRealSpaceSurfaceDescriptors,
)

fod = RunFOD(mol, electronicTemperature=5000.0)         # 与 xtb --fod --etemp 5000 对齐
chi = CalcLinearResponseMatrix(result, scheme="lowdin") # 单点闭壳层线性响应矩阵 chi_AB
ef = CalcBondElectricField(mol, stepAngstrom=1e-3)      # 单次 SCF 批量计算所有键中点局域电场
iso = CalcDensityIsosurface(mol, isovalueAU=0.001)      # 单次 SCF 精确 0.001 a.u. 等值面 + ESP/ALIE/LEA
```

- `RunFOD` / `CalcFODFromResult` 实现 Grimme–Hansen 有限温 Fermi–Dirac 分数占据电子密度（FOD）诊断，返回 `n_fod_e`、自旋分辨费米能级 `fermi_alpha_ev`/`fermi_beta_ev`、每轨道 `orbital_fod_e` 以及原子与 `s/p/d` 壳层 `fod_pop_e`/`fod_s_e`/`fod_p_e`/`fod_d_e`。
- `CalcLinearResponseMatrix` 在 `keepOverlap=True` 的闭壳层单点结果上计算独立粒子原子凝聚线性响应矩阵 $\chi_{AB}^{\mathrm{IP}}$（`chi_e_per_ev`、`chi_e_per_ha`、对角 `chi_diag_e_per_ev`、`atomic_softness_e_per_ev`、`nonlocal_coupling_e_per_ev`），精确满足 $\sum_B \chi_{AB} = 0$，并对开壳层、未知自旋、分数占据或零能隙体系显式返回 `available=False` 与 `missing` 原因。
- `CalcBondElectricField` / `CalcBondElectricFieldFromESP` 在单次 `libxtb` 单点计算中对全部 RDKit 键中点构造 6 点中心有限差分模板求取 $\mathbf{E} = -\nabla V$，输出键向投影 `e_parallel_au`/`e_parallel_mv_cm`、垂直分量 `e_perp_au`/`mv_cm`、模长 `e_norm_au`/`mv_cm` 及三维向量 `e_vec_au`（`1 a.u. = 5142.2067476326 MV/cm`）。
- `CalcRealSpaceProperties`、`CalcRealSpaceSurfaceDescriptors` 与 `CalcDensityIsosurface` 通过扩展版 `libxtb` 的 `realspace_extension` ABI v1 在内存中精确求值笛卡尔 GTO 价电子密度 $\rho(\mathbf{r})$、实空间平均局域电离能 $\mathrm{ALIE}(\mathbf{r})$、虚轨道密度和局域电子亲和能 $\mathrm{LEA}(\mathbf{r})$，并在单次 SCF 内给出精确 `0.001 a.u.` 等值面的分子体积、表面积、球度以及表面 ESP/ALIE/LEA 统计；低于 `densityCutoffAU` 的点显式返回 `None`，绝不填零。

## 轨道加权反应性、ALIE/PDOS、COOP/ECOHP 与 Pipek–Mezey 局域轨道

在同一次 `keepOverlap=True` 单点结果上，可直接提取以下单点轨道与键分析描述符，不再调用引擎：

```python
from eled import (
    XTBCalculator, CalcOrbitalWeightedFukui, CalcAtomicOrbitalEnergies,
    CalcBondOverlapDescriptors, CalcLocalizedOrbitals,
)

result = XTBCalculator(method="GFN2-xTB").Calculate(mol, keepOverlap=True)
owf = CalcOrbitalWeightedFukui(result, deltaEV=0.5, scheme="lowdin")
aoe = CalcAtomicOrbitalEnergies(result, scheme="lowdin", windowEV=2.0)
bond_ovlp = CalcBondOverlapDescriptors(mol, result)
lmo = CalcLocalizedOrbitals(mol, result, scheme="lowdin")
```

- `CalcOrbitalWeightedFukui` 计算以 $\mu = (\mathrm{HOMO}+\mathrm{LUMO})/2$ 为中心的高斯轨道加权
  `f_w_plus`、`f_w_minus`、`f_w_radical` 和 `dual_w`（精确满足求和规则）。
- `CalcAtomicOrbitalEnergies` 返回原子占据/虚轨道布居（`occ_n_e`、`vir_n`）、平均占据轨道能
  `occ_e_mean_ev`（凝聚态 ALIE 代理量的相反数：$\mathrm{ALIE}_A = -\bar{\varepsilon}_{\mathrm{occ},A}$）、
  `vir_e_mean_ev`、前线高斯加权轨道能以及 HOMO/LUMO PDOS 窗口（`pdos_homo_win_e`、`pdos_lumo_win`）。
- `CalcBondOverlapDescriptors` 为每条 RDKit 键计算 `coop_e`（成键为正）、`ecohp_ev`（成键稳定化为负）
  以及引擎/CM5 键极性（`delta_q_e`、`signed_delta_q_e`、`q_prod_e2` 及对应 `*_cm5_*` 字段；GFN0 的 CM5 字段为 `None`）。
- `CalcLocalizedOrbitals` 使用 2×2 Jacobi 旋转执行 Pipek–Mezey 局域化（`"lowdin"` 或对称化 `"mulliken"`），
  输出 `lone_pair`/`bond`/`delocalized` 分类、旋转 Fock 对角能 `lmo_energy_ev`、质心 `centroid_a`
  以及原子孤对数 `n_lone_pairs`、局域键数 `n_localized_bonds` 和最高孤对能 `max_lone_pair_energy_ev`；
  开壳层、未知自旋或分数占据体系显式返回 `None` 及 `missing` 原因。

## 扩展 `libxtb` 性质（能量分解、多极矩、极化率与原子核静电势）

扩展性质 ABI v1 在每次成功的 GFN0/1/2 单点计算中自动导出能量分解、费米能级、原子核静电势（EPN）、
D4 极化率（GFN2/GFN0）与 CAMM 多极矩（GFN2/GFN1），原生与编译 `rdXTBTools` 结果完全等价：

```python
from eled import (
    XTBCalculator, CalcEnergyDecomposition, CalcMultipoleDescriptors,
    CalcPolarizabilityDescriptors, CalcNuclearESP,
)

result = XTBCalculator(method="GFN2-xTB").Calculate(mol)
decomp = CalcEnergyDecomposition(result)
multipoles = CalcMultipoleDescriptors(result, mol=mol)
pol = CalcPolarizabilityDescriptors(result, mol=mol)
epn = CalcNuclearESP(result, mol=mol)
```

不适用的方法项绝不填零：GFN0/1 的 `e_aes_ha`/`e_axc_ha` 为 `None`，GFN2 的 `e_xb_ha` 为 `None`，
气相的 `g_solv_ha` 为 `None`，`electronicTemperature <= 0.1 K` 时费米能级为 `None`，
GFN0 的多极矩返回 `available=False`，GFN1（使用 DFT-D3）的 D4 极化率/C6/C8 返回 `available=False`。
有限电场响应极化率不在 `libxtb` 6.7.1（未启用 `tblite`/PTB）的 GFN0/1/2 支持范围内
（外点电荷仅耦合到单极壳层电荷，无法描述原子内与离平面极化），请直接使用解析 D4 的
`CalcPolarizabilityDescriptors(result)`。

## 隐式溶剂化（ALPB/GBSA）与振动热化学

`XTBOptions`、`XTBCalculator`、`RunXTBBatch` 和 `rdXTBTools.RunMol` 均支持可选的
`solvent`（默认 `None` 保持气相行为与缓存键不变）、`solvationModel`（`"alpb"` 或 `"gbsa"`）和
`solvationGrid`。`CalcSolvationDescriptors` 在固定几何下运行配对的气相与溶剂单点；
`CalcVibrationalThermochemistry` 通过解析梯度 `GetGradient()` 的中心有限差分计算笛卡尔 Hessian 与偶极导数：

```python
from eled import CalcSolvationDescriptors, CalcVibrationalThermochemistry

solv = CalcSolvationDescriptors(mol, "water", method="GFN2-xTB", solvationModel="alpb")
thermo = CalcVibrationalThermochemistry(
    mol, method="GFN2-xTB", temperatureK=298.15, pressureAtm=1.0, symmetryNumber=1
)
```

`CalcSolvationDescriptors` 返回 `delta_g_solv_ha`/`kcalmol`、直接溶剂能 `g_solv_direct_ha`、
溶质极化能 `solute_polarization_energy_ha`、Born/SASA/HB/shift 分项（GFN1/2）、偶极与前线能级位移
以及原子电荷位移（`delta_q_e`、`delta_cm5_q_e`）。`CalcVibrationalThermochemistry` 返回谐振频率
`frequencies_cm1`、红外强度 `ir_intensities_km_mol`、转动常数 `rotational_constants_cm1` 及
RRHO（`rotorCutoffCm1=0`）或 Grimme qRRHO（默认 `rotorCutoffCm1=100.0`）的 ZPE/U/H/S/G，
并通过 `is_stationary_point=False` 和 `warnings` 明确标记非驻点几何或虚频。

## 缺失值和不变性

`None` 表示不可用或定义域不满足，不能当成零。报告中的 `missing` 字段会
记录原因，例如 `basis_not_retained` 或
`spin_density_unavailable_for_open_shell`。标量密度描述符在刚体平移/旋转下保持不变；带电分子的偶极模长依赖坐标原点。
密度描述符在同时变换 `C`、`S`、AO 映射和 p 分量方向时保持 AO 排列和轨道
相位不变。

完整字段、单位、成熟度和缺失语义见 [descriptor-catalog.md](../descriptor-catalog.md)。

π 布居、π 参考电荷和投影键指标的公式及适用范围见
[密度描述符](density-descriptors.md)。构象汇总和敏感性比较见
[构象统计](conformer-sensitivity.md)。
