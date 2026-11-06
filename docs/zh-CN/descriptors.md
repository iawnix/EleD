# 描述符

所有描述符都使用小写 `snake_case` 字段名，并明确单位。描述符函数只分析
已有的 xTB 结果；只有 `run_xtb()`、`run_reactivity()` 等入口会调用引擎。

## 单点描述符

```python
from rdkit_eth import (
    run_xtb, calc_xtb_descriptors, calc_xtb_atom_descriptors,
    calc_xtb_bond_descriptors,
)

result = run_xtb(mol, method="GFN2-xTB")
molecular = calc_xtb_descriptors(result)
atoms = calc_xtb_atom_descriptors(result)
bonds = calc_xtb_bond_descriptors(mol, result)
```

主要字段包括：

- `xtb_total_energy_hartree`：总电子能，Hartree；
- `xtb_homo_ev`、`xtb_lumo_ev`、`xtb_gap_ev`：eV；
- `xtb_dipole_*_e_bohr`、`xtb_dipole_norm_debye`：偶极分量和模长；
- `xtb_charge_*_e`：原子部分电荷统计；
- `partial_charge_e`、`wiberg_sum`：原子级特征；
- `wiberg_order`：对应 RDKit 键的 Wiberg 键级。

HOMO/LUMO 由轨道占据阈值定义。能隙是基态前线轨道差值，不是光学激发能。

## 密度和 shell 描述符

```python
from rdkit_eth import run_xtb, calc_shell_populations, calc_density_matrix

result = run_xtb(mol, method="GFN2-xTB", keep_overlap=True)
shells = calc_shell_populations(result)
density = calc_density_matrix(result)
```

扩展版 `libxtb` 会返回同一次计算中的 `C`、`S`、AO atom map 和 shell map。
程序构造

```text
P = C · diag(occupation) · Cᵀ
```

并检查 `trace(P S)` 是否等于电子数。可用字段包括：

- `mulliken_s/p/d_electron_count_e`；
- `lowdin_s/p/d_electron_count_e`；
- `mulliken_p_to_s_ratio`、`lowdin_p_to_s_ratio`；
- `trace_PS_e`、`electron_conservation_error_e`。

这些是基组相关的密度分区量，不是唯一的实验可观测量。GFN0 的 EEQ 电荷
与 Mulliken/Löwdin 布居不是同一个定义。

## 前线、Mayer 和 π 描述符

```python
from rdkit_eth import (
    calc_frontier_descriptors, calc_mayer_bond_descriptors,
    calc_frontier_bond_descriptors, calc_pi_descriptors,
)

frontier = calc_frontier_descriptors(mol, result)
mayer = calc_mayer_bond_descriptors(mol, result)
frontier_bonds = calc_frontier_bond_descriptors(mol, result)
pi = calc_pi_descriptors(mol, result)
```

这些接口提供简并轨道组布居、participation ratio、熵、Mayer 键级、带符号
的 HOMO/LUMO 键重叠、π 电子数和 π 键级。Mayer 需要闭壳层自旋密度；局部
π 特征需要定义明确的局部平面法向。条件不满足时返回 `None` 和原因。

## 反应性描述符

```python
from rdkit_eth import run_reactivity

report = run_reactivity(mol, population_scheme="lowdin", return_errors=True)
```

固定几何计算 Q、Q+1、Q−1，得到 `vertical_ip_ev`、`vertical_ea_ev`、
`hardness_ev`、`softness_per_ev`、`electrophilicity_ev`，以及原子级
`f_plus`、`f_minus`、`f_radical`、`dual` 和局部软度。该流程支持 GFN1/2，
不自动搜索自旋态，也不进行电荷态几何优化。

## 缺失值和不变性

`None` 表示不可用或定义域不满足，不能当成零。报告中的 `missing` 字段会
记录原因，例如 `basis_not_retained` 或
`spin_density_unavailable_for_open_shell`。标量描述符在刚体平移/旋转下
保持不变；密度描述符在同时变换 `C`、`S` 和 AO 映射时保持 AO 排列和轨道
相位不变。

完整字段、单位、成熟度和缺失语义见 [descriptor-catalog.md](../descriptor-catalog.md)。
