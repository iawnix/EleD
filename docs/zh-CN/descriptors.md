# 描述符

所有描述符都使用小写 `snake_case` 字段名，并明确单位。描述符函数只分析
已有的 xTB 结果；只有 `XTBCalculator.calculate()`、`RunReactivity()` 等入口会调用引擎。

## 单点描述符

```python
from eled import (
    XTBCalculator, CalcXTBDescriptors, CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors,
)

result = XTBCalculator(method="GFN2-xTB").calculate(mol)
molecular = CalcXTBDescriptors(result)
atoms = CalcXTBAtomDescriptors(result)
bonds = CalcXTBBondDescriptors(mol, result)
```

主要字段包括：

- `xtb_energy_ha`：总电子能，Hartree；
- `xtb_homo_ev`、`xtb_lumo_ev`、`xtb_gap_ev`：eV；
- `xtb_dip_*_ebohr`、`xtb_dip_debye`：偶极分量和模长；
- `xtb_q_*_e`：原子部分电荷统计；
- `q_e`、`wbo_sum`：原子级特征；
- `wbo`：对应 RDKit 键的 Wiberg 键级。

YAeHMOP 适配器也使用同一套短字段规则，字段为
`eht_homo_ev`、`eht_lumo_ev`、`eht_gap_ev`、`eht_fermi_ev`、
`eht_energy_ev`、`eht_n_orb` 和 `eht_n_e`。`eht_` 前缀表示这是 RDKit
EHT 计算的结果。

HOMO/LUMO 由轨道占据阈值定义。能隙是基态前线轨道差值，不是光学激发能。

## 密度和 shell 描述符

```python
from eled import XTBCalculator, CalcShellPopulations, CalcDensityMatrix

result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
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

## 前线、Mayer 和 π 描述符

```python
from eled import (
    CalcFrontierDescriptors, CalcMayerBondDescriptors,
    CalcFrontierBondDescriptors, CalcPiDescriptors,
)

frontier = CalcFrontierDescriptors(mol, result)
mayer = CalcMayerBondDescriptors(mol, result)
frontier_bonds = CalcFrontierBondDescriptors(mol, result)
pi = CalcPiDescriptors(mol, result)
```

这些接口提供简并轨道组布居、participation ratio、熵、Mayer 键级、带符号
的 HOMO/LUMO 键重叠、π 电子数和 π 键级。Mayer 使用闭壳层的自旋求和密度；开壳层需要当前未提供的自旋密度。局部
π 特征需要定义明确的局部平面法向。条件不满足时返回 `None` 和原因。

## 反应性描述符

```python
from eled import RunReactivity

report = RunReactivity(mol, population_scheme="lowdin", return_errors=True)
```

固定几何计算 Q、Q+1、Q−1，得到 `ip_v_ev`、`ea_v_ev`、
`mu_ev`、`chi_ev`、`eta_ev`、`s_1ev`、`omega_ev`，以及原子级
`f_plus`、`f_minus`、`f_radical`、`dual`，以及 `s_loc_1ev`、
`omega_loc_ev` 和 `n_loc_1ev`。该流程支持 GFN1/2，
不自动搜索自旋态，也不进行电荷态几何优化。

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
