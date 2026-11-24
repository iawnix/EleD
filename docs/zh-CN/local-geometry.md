# 局部位阻与环几何描述符

EleD 分析具有显式氢和三维构象的 RDKit 分子。所有索引沿用原始原子编号。
五个接口都接受 `conf_id=-1`，表示使用 RDKit 默认构象，并记录实际构象 ID。
计算只需要 RDKit 和 NumPy，使用已有坐标。

```python
from eled import PrepareConformers, CalcHOMA

prepared = PrepareConformers('c1ccc2ccccc2c1')
if not prepared.success:
    raise RuntimeError(prepared.report)
mol = prepared.molecule
rings = mol.GetRingInfo().AtomRings()
reports = [CalcHOMA(mol, ring) for ring in rings]
print([r['values']['homa'] for r in reports])
```

稠环可在不同调用中共享原子和键。环输入是按相邻键依次排列的、原子各出现一次的
闭合路径，末原子必须与首原子成键。由调用者选择分析哪些环。

## 接口和返回记录

| 接口 | 选择对象 | 主要结果 |
| --- | --- | --- |
| `CalcSterimol(mol, atom_indices, axis=(base, attached))` | 取代基和有向连接键 | `l_a`、`b1_a`、`b5_a`，Å |
| `CalcBuriedVolume(mol, atom_indices, center=i)` | 占据空间的原子和被排除的中心原子 | `vbur_percent`，百分数 |
| `CalcPlaneDescriptors(mol, atom_indices, other_atom_indices=None)` | 一个或两个至少含三个原子的集合 | `plane_rms_a`、`plane_max_a`，可选 `centroid_distance_a`，Å；`plane_angle_deg`，度 |
| `CalcHOMA(mol, ring_atom_indices)` | 有序环路径 | 无量纲 `homa` 和逐键 `length_a`，Å |
| `CalcRingPuckering(mol, ring_atom_indices)` | 至少四元的有序环 | `q_total_a`、偶数环 `q_alternating_a`，Å；六元环 `theta_deg`；各模式的 `q_a` 和 `phi_deg` |

报告使用 `schema_version=1`、对应的 `report_type` 和 `protocol_version='1.0'`，
包含 `values`、`units`、`missing`、`selection`、`analysis_settings`、`provenance`。
来源信息记录实际构象 ID、原子顺序与坐标哈希、原子序数和坐标。
未定义值为 `None` 并附原因，输入错误抛出 `ValueError`。
用 `json.dumps(report, allow_nan=False)` 保存 JSON 或 JSONL。
`ExportDescriptorReports` 使用电子描述符报告 schema，局部几何报告直接保存为 JSON。

位阻默认使用 RDKit van der Waals 半径，`radius_scale=1.0`。
其他半径体系通过完整原子长度的 `radii_angstrom` 向量传入。
报告保存有效半径和 RDKit 版本。氢是否参与由选择集合显式决定。
比较结果时应保持半径、氢处理、几何准备和描述符参数一致。

## Sterimol

轴由 base 指向 attached，两原子必须成键。集合必须包含 attached 并排除 base。
每个被选原子是一个给定半径的球，相对位置以 base 为原点。

- `L=max(轴向位置+半径)+l_correction_angstrom`，默认修正 0.4 Å。
- `B1=min_方向 max_原子(横向位置·方向+半径)`。
- `B5=max_原子(到轴的距离+半径)`。

B1 默认沿横截面整圆均匀采样 `angular_samples=7200` 个方向，离散计算相对精确
极小值的正偏差不超过 `最大横向距离*pi/采样数`，记录在
`b1_error_bound_angstrom`。增加采样数可缩小误差界。
一个轴上原子球给出 B1=B5=球半径。与文献结果比较时需核对半径体系和 L 的定义。

## 埋藏体积

以 `center` 原子为球心，`sphere_radius_angstrom=3.5`，计算所选原子球并集占据的
球内比例。中心原子必须从集合中排除，所选氢原子参与体积计算。

采用以球心为中心的 Cartesian 中点网格，默认请求间距 `spacing_angstrom=0.1`。
球直径分成整数个网格单元，记录有效间距、球内点数和占据点数，百分数为
`100*占据点数/球内点数`，重叠部分只计一次。有限网格下旋转不变性有数值误差，
可减小间距检查结果敏感性。支持直径方向 4–256 个单元，按平面处理控制内存。
与配体 %Vbur 文献比较时还需要匹配半径缩放及氢处理。

## 平面与环起伏

平面采用不加权最小二乘。`plane_rms_a` 和 `plane_max_a` 分别是到平面的均方根
和最大绝对距离。两个无向平面的夹角为 0–90 度。
原子共线或最小特征值简并时，法向不唯一，返回 `plane_normal_not_unique`。

环起伏采用 Cremer–Pople 定义（JACS 1975, 97, 1354–1358）。对中心化坐标 r_j，
j=0..N−1，法向平行于
`sum(r_j*sin(2*pi*j/N)) × sum(r_j*cos(2*pi*j/N))`，沿法向投影得到 z_j。
m=2..floor((N−1)/2) 模式的余弦和正弦分量分别为
`sqrt(2/N)*sum(z_j*cos(2*pi*m*j/N))` 和
`-sqrt(2/N)*sum(z_j*sin(2*pi*m*j/N))`。
偶数环另有带符号模式 `sum((-1)^j*z_j)/sqrt(N)`。
总起伏 Q 为 `sqrt(sum(z_j**2))`。六元环还返回
`theta=atan2(q2,q3)`，范围 0–180 度。相位为 0–360 度，依赖环的起点和方向。
幅度小于 1e-10 Å 的模式相位为 `None`，平面六元环的 theta 也为 `None`。
法向无法定义时，各相关字段都记录缺失原因。

## HOMA

`HOMA=1-mean(alpha_AB*(R_AB-Ropt_AB)**2)`。
默认 C–C 参数为 Ropt=1.388 Å、alpha=257.7 Å⁻²，来源是
Krygowski, J. Chem. Inf. Comput. Sci. 1993, 33, 70–78。
保留负值，计算包括闭合边在内的全部环键。计算几何指标无需 RDKit 芳香性标签。

杂原子环可指定完整的 `parameters={(6,6):(Ropt_CC,alpha_CC),
(6,7):(Ropt_CN,alpha_CN)}` 和非空 `parameter_source`。
自定义表替代内置参数表。任何元素对缺少参数时，整个 HOMA 返回 `None`，
原因是 `unsupported_element_pair`，同时保留全部测量键长。
HOMA 依赖几何准备和参考参数，反映芳香性的几何特征，需结合具体化学体系解释。

## 多构象

用相同原子集合和参数，逐一传入实际构象 ID，保留各构象报告。例如取代基的
B5 可按下式计算等权均值和总体标准差。

```python
import numpy as np
from eled import CalcSterimol

rows = [CalcSterimol(mol, substituent, axis=axis, conf_id=c.GetId())
        for c in mol.GetConformers()]
b5 = np.array([row['values']['b5_a'] for row in rows])
mean, std = float(b5.mean()), float(b5.std())
# 给定构象权重 w 时：
# mean = float(np.average(b5, weights=w))
# std = float(np.sqrt(np.average((b5-mean)**2, weights=w)))
```

汇总其他指标时保留缺失原因和覆盖率。相位角需要圆统计，Boltzmann 权重需要
适用的构象能量与明确的温度。

[可运行示例](../../examples/local_geometry.py) 展示苯甲醚的位阻、萘的两个环和
环己烷起伏。测试使用单球 Sterimol、两球相交体积、理想键长及独立傅里叶环模式，
并检查旋转、平移和原子编号映射。
