# 描述符目录

项目的机器可读目录由 `eled.descriptor_catalog()` 提供，成熟度由
`eled.descriptor_status()` 提供。字段统一使用小写 `snake_case`。

## 推荐特征

适合作为第一版模型输入：

- `xtb_total_energy_hartree`、`xtb_homo_ev`、`xtb_lumo_ev`、`xtb_gap_ev`；
- 偶极分量和模长；
- 原子 `partial_charge_e`、`wiberg_sum`；
- RDKit 键的 `wiberg_order`；
- Mulliken/Löwdin s/p/d 布居和 p/s 比；
- GFN1/2 固定几何 ΔSCF 的 IP、EA、硬度、软度和 Fukui 指数。

## 诊断字段

诊断字段用于记录可靠性，不默认作为模型输入：

- 前线轨道索引、占据阈值和电子/轨道计数；
- 构象离散度；
- GFN1 与 GFN2 的差异；
- 300 K 与 1000 K 电子温度差异。

## 探索性字段

Mayer 键级、带符号的前线键重叠、π 布居、π 键级和局部亲核性代理量属于
探索性字段。它们定义完整，但目前尚未证明在跨骨架任务中具有稳定预测增益。

## 单位和缺失值

- 能量：Hartree 或 eV；
- 电荷和布居：`e`；
- 偶极：`e*Bohr` 或 Debye；
- 键级、比例和 Fukui 指数：无量纲；
- 原子和键索引：RDKit 原始编号。

`None` 表示不可用或定义域不满足，不能当成零。报告的 `missing` 区段保存
原因，例如 `basis_not_retained`、`local_plane_normal_undefined` 或
`spin_density_unavailable_for_open_shell`。

密度描述符在同时变换 AO 映射、`C` 和 `S` 时对 AO 排列和轨道相位不变；
标量能量和布居在刚体平移/旋转下不变。完整英文目录见
[descriptor-catalog.md](../descriptor-catalog.md)。
