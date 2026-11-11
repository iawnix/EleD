# 构象统计与敏感性（1.3）

`RunConformerAnalysis` 对输入分子的指定构象计算描述符，默认等权统计。
原子和键编号在构象间固定。该函数读取已有三维构象；构象生成或优化由
[构象准备](geometry-preparation.md) 完成。

```python
from eled import RunConformerAnalysis

report = RunConformerAnalysis(mol, boltzmann_weighted=True,
                              boltzmann_temperature=300.0)
gap = report["statistics"]["xtb_gap_ev"]
charge = report["atom_statistics"]["0"]["atoms.q_e"]
homo = report["bond_statistics"]["0"]["frontier_bonds.homo.mul_ovlp_e"]
lumo = report["bond_statistics"]["0"]["frontier_bonds.lumo.mul_ovlp_e"]
pi_bond = report["bond_statistics"]["0"]["pi.pi_bo"]
```

π 分析通过 `pi_scheme`、`pi_plane_tol_a` 和 `pi_ref` 配置，详见
[报告参数](descriptor-reports.md#π-分析参数)。同一组参数应用于所有输入构象，
并记录在汇总报告和逐构象报告中；等权和电子能加权统计均使用该 π 定义。

## 等权统计

分子级 `statistics`、原子级 `atom_statistics`、键级 `bond_statistics`
提供 `mean`、`std`、`min`、`max` 及覆盖计数。分子级还提供 `range`。
原子和键的 JSON 索引键为字符串；`units`、`atom_units`、`bond_units`
分别记录三个层级的单位。

每个构象对一个索引字段最多贡献一个值，HOMO/LUMO 键重叠分开统计，π 原子
布居和 π 键指标分开统计。`n_total = n_valid + n_missing + n_failed`。
全部缺失时统计值为 `None`；只有一个有效观测时，总体标准差为零（`ddof=0`），
这不代表已经覆盖构象空间。重复输入的构象仍分别计权。

`frontier_group_composition_changed` 表示前线轨道组的索引集合发生变化。
相同索引不等于完成了跨构象轨道追踪；每个构象的报告保留实际分组信息。
笛卡尔分量、轨道索引和带电分子的偶极模长不进入默认统计。

## 电子能加权

可选的 `boltzmann_weighted=True` 增加独立的 `boltzmann` 区段，保留默认的
等权统计。权重为：

```text
w_i ∝ exp[−(E_i − E_min)/(k_B T)]
```

`E_i` 来自同一计算协议下的 xTB 电子能。这是**电子能权重，不是自由能权重**；
未计入振动熵、零点能或溶剂效应。`boltzmann_temperature` 是构象汇总温度，
与控制轨道占据平滑的 xTB `electronic_temperature` 分开设置。

`boltzmann.statistics`、`boltzmann.atom_statistics` 和
`boltzmann.bond_statistics` 给出加权均值、标准差及覆盖计数。
`weight_sum` 保存该字段可用观测的权重总和，再以该总和归一化。
权重下溢为零的有限观测仍计入有效数；如果该字段可用观测的总权重为零，
加权结果为 `None`，原因是 `no_positive_weight_for_available_values`。
单位与等权统计一致。

## 方法和电子温度比较

`DescriptorSensitivity` 从已有报告提取构象标准差、GFN1/GFN2 绝对差值，
以及 300 K/1000 K 电子温度绝对差值。它核对有序几何和计算协议，并逐字段
要求有效构象子集匹配，避免将不同子集的均值解释为方法差异。
详细约定见[密度描述符](density-descriptors.md#敏感性诊断)。
