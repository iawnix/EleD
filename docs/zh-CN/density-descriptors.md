# 密度描述符

EleD 从保留的 xTB 结果中提取 shell 布居、键级和局部 π 特征，不增加引擎调用。
`keep_overlap=True` 保存同一次计算的重叠矩阵 `S`、系数 `C`、占据数 `n`、
AO/shell 映射及 p 轨道分量方向。

```python
from eled import XTBCalculator, CalcDensityDescriptors, CalcPiDescriptors

result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
density = CalcDensityDescriptors(mol, result)
pi = CalcPiDescriptors(mol, result, scheme="lowdin", plane_tol_a=0.1)
```

原子和键沿用 RDKit 输入编号。接受 `mol` 的分析函数会核对拓扑、原子顺序和
坐标。保留 overlap 后，`CalcDescriptorReport` 默认包含这些分析区段。

## 密度与 shell 布居

对自旋求和的占据数定义：

```text
P = C diag(n) Cᵀ
M = P S
D = S½ P S½
N = tr(P S)
```

分析前检查 `Cᵀ S C = I`。`P` 的矩阵元素随 AO 基组变换，`tr(P S)` 保持
不变，单位为电子数。`n_err_e = tr(P S) - n_e`，其中 `n_e` 优先采用引擎
报告的模型电子数；缺少该计数时采用轨道占据之和。

一个 shell 的 Mulliken 布居是其在 `M` 中的对角元素之和；Löwdin 布居是
其在 `D` 中的对角元素之和。原子行按角动量汇总所属 shell。

| 字段 | 含义 | 单位 |
| --- | --- | --- |
| `mul_s_e`、`mul_p_e`、`mul_d_e` | Mulliken s/p/d 布居 | e |
| `low_s_e`、`low_p_e`、`low_d_e` | Löwdin s/p/d 布居 | e |
| `mul_ps`、`low_ps` | p 布居 / s 布居 | 无量纲 |
| `mul_qrec_e`、`low_qrec_e` | 参考价电子数减去相应总布居 | e |
| shell 行的 `mul_n_e`、`low_n_e` | 单个 shell 的布居 | e |
| `tr_ps_e`、`n_e`、`n_err_e` | 密度电子数、参考电子数、两者差值 | e |
| `shell_err_e` | 两种定义各自的 shell 到原子求和最大偏差 | e |

不存在的 s/p/d 通道为零；s 布居为零时，p/s 比为 `None`。此接口支持至 d
壳层。电荷重建的参考价电子数覆盖 H、B、C、N、O、F、Si、P、S、Cl、Ge、
As、Se、Br、I；其他元素返回 `None` 和 `reference_valence_unavailable`。
GFN1/2 的 Mulliken 重建电荷可以与引擎电荷核对；Löwdin 电荷采用另一种分区。
GFN0 的 EEQ 电荷来自不同模型，因此两个重建电荷字段均为 `None`。

## Mayer 与前线轨道键描述符

闭壳层 RDKit 键 A–B 的 Mayer 键级为：

```text
mbo(A,B) = Σ[μ∈A,ν∈B] Mμν Mνμ
```

单位为无量纲，与引擎的 `wbo` 并列输出。项目固定版本的 xTB 在闭壳层下
使用同样的矩阵收缩计算 Wiberg，回归测试检查两者相等；该关系依赖这一
具体实现约定。开壳层或自旋未知时，由于缺少自旋密度，`mbo` 为 `None`。

对 HOMO/LUMO 轨道组 G，定义：

```text
PG = C_G C_Gᵀ / |G|
mul_ovlp_e(A,B) = 2 Σ[μ∈A,ν∈B] PGμν Sμν
```

这是**按一个电子归一化的组平均值**，单位为 e，不是整个已占据轨道组的
总电子布居。正值表示成键重叠，负值表示反键重叠。`frontier_group` 和
`orb_idx` 记录轨道组，能量及占据容差与 `CalcFrontierDescriptors` 共用。
轨道相位变化和选定子空间内的旋转不改变组投影结果。缺少轨道组时，仍保留
对应键行，重叠值为 `None`，原因是 `frontier_group_unavailable`。

## 局部 π 投影

`CalcPiDescriptors` 适用于具有一个完整 p 壳层、且局部平面明确的 sp²/芳香
原子。环原子的拟合平面使用所有包含该原子的环的原子集合；其他原子使用
自身及相邻原子。至少需要三个不共线点，所有点到拟合平面的最大距离不得
超过 `plane_tol_a`，默认 0.1 Å。

结果的 `ao_p_directions` 将 p 分量方向附在 AO 行上。局部单位法向据此
映射为仅在原子 A 的 p 壳层上非零的向量 `vA`。AO 排列或符号变化时，同步
变换方向信息即可保持投影定义。`scheme` 支持默认的 `"lowdin"` 和 `"mulliken"`：

| 字段 | Löwdin | Mulliken | 单位 |
| --- | --- | --- | --- |
| `pi_n_e` | `vAᵀ D vA` | `vAᵀ sym(M) vA` | e |
| `pi_q_e` | `nπ,ref − pi_n_e` | `nπ,ref − pi_n_e` | e |
| `pi_bo` | `(vAᵀ D vB)²` | `(vAᵀ M vB)(vBᵀ M vA)` | 无量纲 |

π 参考电子数默认仅对**中性、无自由基的 sp²/芳香碳取 1 e**。其他原子的
`pi_n_e` 仍可计算；`pi_q_e` 则需要显式参考值，例如 `pi_ref={atom_index: 2.0}`。
参考值范围为 0 至 2 e。原子的总价电子数不能作为 π 参考电子数。
返回的 `pi_ref` 使用 JSON 字符串键记录显式覆盖值。`CalcDensityDescriptors`、
报告、批处理和构象分析均支持 `pi_scheme`、`pi_plane_tol_a`、`pi_ref`，
详见[报告参数](descriptor-reports.md#π-分析参数)。

`pi_bo` 是投影键指标，不是重叠布居的绝对值，也不是总键级的完整分解。
它要求闭壳层且两端投影均有效。任一法向反号均不改变结果。Löwdin 表达式
非负，Mulliken 的乘积可以带符号。

所有原子和 RDKit 键均保留行。不可用字段为 `None`，并在 `missing` 中逐项
记录原因，例如 `local_plane_degenerate`、`local_plane_nonplanar`、
`p_component_directions_unavailable`、`pi_reference_unspecified` 或
`spin_density_unavailable_for_open_shell`。返回的 `normal` 是坐标元数据，
不作为标量描述符。

## 不变性与验证

标量密度描述符在刚体运动、MO 相位变化，以及同步变换 `C`、`S`、映射和分量
方向的 AO 排列/符号变化下不变。`P` 矩阵自身随基组变换；相同占据轨道间的
旋转保持密度不变；完整前线轨道组内部旋转保持组平均重叠不变。局部法向的
符号不影响 π 数值。

测试覆盖非正交解析密度、电子数和 shell 守恒、无效平面、原子/坐标来源检查、
简并旋转，以及 GFN0/1/2 的 ctypes/compiled 真实计算。苯、吡啶、乙烯、呋喃
还用于验证 π 描述符的刚体运动不变性。

密度分析基于 Mulliken（1955）、Löwdin（1950）和 Mayer（1983）的定义。
`DescriptorCatalog()` 提供机器可读的定义和单位；
[文献定位](../literature-density-ml-audit.md) 单独记录。

## 敏感性诊断

`DescriptorSensitivity(reference, gfn1=..., gfn2=..., temperature_300=...,
temperature_1000=...)` 读取已有的单点或构象报告，分别输出构象标准差、
GFN1/GFN2 绝对差值和 300 K/1000 K 电子温度绝对差值，不重新计算。

比较要求原子顺序、几何、电荷/自旋、引擎、分析设置及失败构象位置一致。
方法比较保持数值设置相同；温度比较保持方法和参数相同。逐字段检查有效
构象子集，覆盖不一致时返回 `None` 和 `unmatched_conformer_coverage`。
未提供的比较保留明确缺失信息。

输出包含分子级 `values`、按编号排列的 `atoms`/`bonds`，以及 `units`、
`atom_units`、`bond_units`。不统计笛卡尔坐标/偶极分量、轨道索引和带电
分子的偶极模长。这些量用于描述协议敏感性，不合成为单一可信度评分。
