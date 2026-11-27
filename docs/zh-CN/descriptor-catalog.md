# 描述符目录

项目的机器可读目录由 `eled.DescriptorCatalog()` 提供，成熟度由
`eled.DescriptorStatus()` 提供。公开描述符函数与类方法采用 RDKit 风格的
`PascalCase`（`Calc...`、`Run...`、`Calculate`、`Get*`），Python 参数与结果属性采用
`camelCase`（`confId`、`keepOverlap`、`totalEnergyHartree`）；导出的报告/CSV 特征列名统一采用带单位后缀的小写 `snake_case`。

## 推荐特征

适合作为第一版模型输入：

- `xtb_energy_ha`、`xtb_homo_ev`、`xtb_lumo_ev`、`xtb_gap_ev`；
- 偶极分量和模长；
- 原子 `q_e`、`cm5_q_e`、`wbo_sum`；`q_e` 的确切模型写入 `charge_model`，`cm5_q_e` 表示几何校正后的 CM5 电荷；
- RDKit 键的 `wbo`；
- Mulliken/Löwdin s/p/d 布居和 p/s 比；
- GFN1/2 固定几何 ΔSCF 的 IP、EA、硬度、软度和 Fukui 指数；
- Politzer–Murray 表面静电势统计（`esp_surface`）与原子核静电势 `epn_au`；
- D4 极化率/C6/C8、能量分解（`properties`）、ALPB/GBSA 隐式溶剂化位移（`solvation`）及谐振 RRHO/qRRHO 热化学（`thermochemistry`）。

## 诊断字段

诊断字段用于记录可靠性，不默认作为模型输入：

- 前线轨道索引、占据阈值和电子/轨道计数；
- 费米能级与直接溶剂能/溶质极化能拆解；
- 构象离散度；
- GFN1 与 GFN2 的差异；
- 300 K 与 1000 K 电子温度差异。

## 探索性字段

Mayer 键级、带符号的前线键重叠、π 布居、π 键级、局部亲核性代理量、轨道加权 Fukui、凝聚态 ALIE/PDOS、COOP/ECOHP/键极性、Pipek–Mezey 局域轨道（`localized_orbitals`）及 CAMM 原子/分子多极矩属于探索性字段。它们定义完整，但目前尚未证明在跨骨架任务中具有稳定预测增益。

`DescriptorCatalog()` 的顶层区段包括 `molecular`、`atom`、`fragment`、`bond`、`density`、`reactivity`、`local_geometry`、`esp_surface`、`localized_orbitals`、`properties`、`solvation`、`thermochemistry` 与 `diagnostic`。

## 单位和缺失值

- 能量：Hartree 或 eV；
- 电荷和布居：`e`；
- 偶极：`e*Bohr` 或 Debye；
- 键级、比例和 Fukui 指数：无量纲；
- 原子和键索引：RDKit 原始编号。

`None` 表示不可用或定义域不满足，不能当成零。报告的 `missing` 区段保存
原因，例如 `basis_not_retained`、`local_plane_nonplanar` 或
`spin_density_unavailable_for_open_shell`。

密度描述符在同时变换 AO 映射、p 分量方向、`C` 和 `S` 时对 AO 排列和轨道相位不变；
标量能量和布居在刚体平移/旋转下不变。完整英文目录见
[descriptor-catalog.md](../descriptor-catalog.md)。
