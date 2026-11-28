# 局部位阻与环几何描述符

EleD 分析具有显式氢和三维构象的 RDKit 分子。所有索引沿用原始原子编号。
局部几何接口都接受 `confId=-1`，表示使用 RDKit 默认构象，并记录实际构象 ID。
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

多个局部描述符可以使用同一套选择对象。

```python
from eled import AnalyzeLocalGeometry, RunLocalGeometryBatch

selections = {
    'atoms': {'substituent': [0, 1, 2], 'ring_atoms': [3, 4, 5, 6, 7, 8]},
    'reactionSites': {'electrophilic_carbon': [3]},
    'bonds': {'attachment': [(2, 3)]},
    'rings': {'aromatic_ring': [3, 4, 5, 6, 7, 8]},
}
one = AnalyzeLocalGeometry(mol, selections=selections,
    sterimol={'substituent': {'axis': (3, 2)}},
    homa={'aromatic_ring': {}})
many = RunLocalGeometryBatch([mol], selections,
    homa={'aromatic_ring': {}})
```

同一选择格式也可通过 `NormalizeSelections` 使用，索引沿用原始原子编号。
反应位点是命名的原子集合，可以和 Fukui 描述符一起保存。批量结果保留成功的
描述符，并在 `errors` 中记录单个描述符的错误。

稠环可在不同调用中共享原子和键。环输入是按相邻键依次排列的、原子各出现一次的
闭合路径，末原子必须与首原子成键。由调用者选择分析哪些环。

## 接口和返回记录

| 接口 | 选择对象 | 主要结果 |
| --- | --- | --- |
| `CalcSterimol(mol, atomIndices, axis=(base, attached))` | 取代基和有向连接键 | `l_a`、`b1_a`、`b5_a`，Å |
| `CalcBuriedVolume(mol, atomIndices, center=i, zAxisAtom=None, xzPlaneAtom=None)` | 占据空间的原子、被排除的中心原子及可选分子对齐坐标系 | `vbur_percent`，以及可选的四象限 `vbur_q1_percent`..`vbur_q4_percent` 与半球 `vbur_z_pos_percent`/`vbur_z_neg_percent`，百分数 |
| `CalcPyramidalization(mol, atomIdx, neighborIndices=None)` | 三配位中心原子（或显式指定 3 个相邻原子） | 无量纲 `pyramidalization_p`、`pyramidalization_angle_deg`、`poav_angle_deg`、`angle_sum_deg`（度）及 `pyramid_height_a`（Å） |
| `CalcConeAngle(mol, atomIndices, center=i, axisAtom=None)` | 配体原子、被排除的顶点 `center` 及可选固定轴 `axisAtom` | `cone_angle_deg`、`half_angle_deg`（度），以及 `cone_axis` 与 `tangent_atom_indices` |
| `CalcPlaneDescriptors(mol, atomIndices, otherAtomIndices=None)` | 一个或两个至少含三个原子的集合 | `plane_rms_a`、`plane_max_a`，可选 `centroid_distance_a`，Å；`plane_angle_deg`，度 |
| `CalcHOMA(mol, ringAtomIndices)` | 有序环路径 | 无量纲 `homa` 和逐键 `length_a`，Å |
| `CalcRingPuckering(mol, ringAtomIndices)` | 至少四元的有序环 | `q_total_a`、偶数环 `q_alternating_a`，Å；六元环 `theta_deg`；各模式的 `q_a` 和 `phi_deg` |
| `CalcMCI(mol, result, ringAtomIndices)` | 有序环（3–9 元）和保留的闭壳层密度数据 | 无量纲 `mci` 与 `i_ring` |

报告使用统一的 `schema_version=2`、`descriptor_version`、对应的 `report_type` 和
`protocol_version='1.0'`，包含 `values`、`units`、`missing`、`selection`、
`analysis_settings`、`provenance`。
来源信息记录实际构象 ID、原子顺序与坐标哈希、原子序数和坐标。
未定义值为 `None` 并附原因，输入错误抛出 `ValueError`。
用 `json.dumps(report, allow_nan=False)` 保存 JSON 或 JSONL。
`ExportDescriptorReports` 使用电子描述符报告 schema，局部几何报告直接保存为 JSON。

位阻默认使用 RDKit van der Waals 半径，`radiusScale=1.0`。SambVca 文献中的 %Vbur 使用 Bondi 半径 ×1.17，与之比较时请传入 `radiusScale=1.17`。
其他半径体系通过完整原子长度的 `radiiAngstrom` 向量传入。
报告保存有效半径和 RDKit 版本。氢是否参与由选择集合显式决定。
比较结果时应保持半径、氢处理、几何准备和描述符参数一致。

## Sterimol

轴由 base 指向 attached，两原子必须成键。集合必须包含 attached 并排除 base。
每个被选原子是一个给定半径的球，相对位置以 base 为原点。

- `L=max(轴向位置+半径)+lCorrectionAngstrom`，默认修正 0.4 Å。
- `B1=min_方向 max_原子(横向位置·方向+半径)`。
- `B5=max_原子(到轴的距离+半径)`。

B1 默认沿横截面整圆均匀采样 `angularSamples=7200` 个方向，离散计算相对精确
极小值的正偏差不超过 `最大横向距离*pi/采样数`，记录在
`b1_error_bound_angstrom`。增加采样数可缩小误差界。
一个轴上原子球给出 B1=B5=球半径。与文献结果比较时需核对半径体系和 L 的定义。

## 埋藏体积

以 `center` 原子为球心，`sphereRadiusAngstrom=3.5`，计算所选原子球并集占据的
球内比例。中心原子必须从集合中排除，所选氢原子参与体积计算。

采用以球心为中心的 Cartesian 中点网格，默认请求间距 `spacingAngstrom=0.1`。
球直径分成整数个网格单元，记录有效间距、球内点数和占据点数，百分数为
`100*占据点数/球内点数`，重叠部分只计一次。有限网格下旋转不变性有数值误差，
可减小间距检查结果敏感性。支持直径方向 4–256 个单元，按平面处理控制内存。
与配体 %Vbur 文献比较时还需要匹配半径缩放及氢处理。传入 `zAxisAtom`（及可选 `xzPlaneAtom`）时，`CalcBuriedVolume` 同时输出以 `center -> zAxisAtom` 为 +z 轴的 SambVca 四象限埋藏体积（`vbur_q1_percent`..`vbur_q4_percent`）和正负半球占比（`vbur_z_pos_percent`、`vbur_z_neg_percent`）。

## 锥体化与配体锥角

`CalcPyramidalization(mol, atomIdx, *, neighborIndices=None, confId=-1)` 计算三配位中心原子的 Haddon POAV1 锥体化指标：`pyramidalization_p = sin(theta_p)`（平面为 `0`，理想四面体为 `1/3`）、`pyramidalization_angle_deg`（平面 `0°`，理想四面体 `19.4712°`）、`poav_angle_deg`、`angle_sum_deg` 及中心到邻居平面的距离 `pyramid_height_a`。退化几何返回 `None` 并写入 `missing` 原因。

`CalcConeAngle(mol, atomIndices, *, center, axisAtom=None, ...)` 从顶点 `center` 计算包围所选配体原子范德华球面的精确 Bilbrey–Allen（2013）最小外接圆锥角（传入 `axisAtom` 时按 `center -> axisAtom` 固定轴计算 Tolman 轴向锥角）。若任一所选原子球面包围顶点（`r_i >= d_i`），显式返回 `cone_angle_deg=None` 及 `apex_inside_atom_sphere`。二者均可通过 `AnalyzeLocalGeometry(..., pyramidalization=..., coneAngle=...)` 统一调度。

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
内置表包含 C–C `(1.388,257.7)`、C–N `(1.334,93.52)`、C–O `(1.265,157.38)`
和 C–S `(1.677,94.09)`，单位分别为 Å 和 Å⁻²，来源是 Krygowski, J. Chem.
Inf. Comput. Sci. 1993, 33, 70–78。
保留负值，计算包括闭合边在内的全部环键。计算几何指标无需 RDKit 芳香性标签。

内置参数包括 C-C `(1.388,257.7)`、C-N `(1.334,93.52)`、C-O `(1.265,157.38)`
和 C-S `(1.677,94.09)`，单位分别为 Å 和 Å⁻²，参数来源为 Krygowski 1993。
其他参数体系可指定完整的 `parameters={(6,6):(Ropt_CC,alpha_CC),
(6,7):(Ropt_CN,alpha_CN)}` 和非空 `parameterSource`。
自定义表替代内置参数表。任何元素对缺少参数时，整个 HOMA 返回 `None`，
原因是 `unsupported_element_pair`，同时保留全部测量键长。
HOMA 依赖几何准备和参考参数，反映芳香性的几何特征，需结合具体化学体系解释。

`CalcMCI` 使用同一次计算保留的闭壳层密度和重叠矩阵，按 Giambiagi–Bultinck
定义对 `M = P*S`（Mulliken）或 `M = S^(1/2)*P*S^(1/2)`（Löwdin）执行 AO 子块矩阵
连乘并求迹，计算 `i_ring` 与首原子固定的全部环排列之和 `mci = (1/4) * sum Tr(M_{A1,A2}...M_{An,A1})`。
该矩阵迹收缩严格满足刚体旋转/平移、AO 正交/符号/置换、MO 相位以及环起点与方向不变性；
开壳层或未知自旋体系显式返回 `None` 及 `missing` 原因。局部 `pi_bo` 是单键描述符，
不能通过取平均得到 MCI。

## 多构象

用相同原子集合和参数，逐一传入实际构象 ID，保留各构象报告。例如取代基的
B5 可按下式计算等权均值和总体标准差。

```python
import numpy as np
from eled import CalcSterimol

rows = [CalcSterimol(mol, substituent, axis=axis, confId=c.GetId())
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
