# 可视化

EleD 的可视化层消费已有的 RDKit 分子以及描述符报告（`CalcDescriptorReport`、
`RunReactivity`、`CalcDensityDescriptors`）或直接传入的 `XTBResults`，不会重新启动
xTB 计算。RDKit 提供二维布局、原子/键编号与分子绘图；EleD 提供电子描述符数值、
连续色阶与片段分析视图。

## 分子与键级描述符视图

`DrawDescriptorMap` 既支持带片段划分的报告，也支持无片段的标准单分子报告或
`XTBResults`。当输入分子带有三维坐标时，EleD 会在内部副本上生成清晰的二维绘图坐标，
不会修改调用方的分子对象。默认隐藏显式氢原子，并通过原子属性精确追踪原始编号（即使
RDKit 保留同位素氢、带电氢或立体化学氢也不会发生编号错位）；数值默认以 RDKit
`atomNote` / `bondNote` 形式标注，保留 `O`、`N`、`Cl` 等元素符号。

```python
from eled import DrawDescriptorMap

# 直接从 XTBResults 或无片段报告绘制电荷连续色阶图
charge_svg = DrawDescriptorMap(mol, result, atomField="q_e")

# 绘制键级或密度派生键描述符（如 wbo、mbo、pi_bo、frontier_bonds.homo.mul_ovlp_e）
bond_svg = DrawDescriptorMap(mol, report, bondField="wbo", colorBy="value")

# 传入 RDKit MolDrawOptions 回调自定义绘图样式
styled_svg = DrawDescriptorMap(
    mol,
    report,
    atomField="cm5_q_e",
    colorBy="value",
    labelStyle="note",      # "note"（默认）、"replace" 或 "none"
    addAtomIndices=True,
    includeHydrogens=False,
    showColorbar=True,
    drawOptions=lambda opts: setattr(opts, "bondLineWidth", 2.0),
)
```

支持的字段包括顶层原子字段（`q_e`、`cm5_q_e`、`wbo_sum`、`f_plus`、`f_minus`、
`f_radical`、`dual`、`s_loc_1ev`、`omega_loc_ev`、`n_loc_1ev`）、前线/壳层/π 嵌套
原子字段（`homo_population`、`lumo_population`、`mul_p_e`、`low_p_e`、`low_ps`、
`pi_n_e`、`pi_q_e`）以及键级字段（`wbo`、`mbo`、`pi_bo`、
`frontier_bonds.homo.mul_ovlp_e`、`frontier_bonds.lumo.mul_ovlp_e`），同时也支持带区段
前缀的路径写法（如 `pi.pi_q_e`、`mayer_bonds.mbo`）。当报告包含 `provenance` 时，
`BuildFragmentVisualData` 与 `DrawDescriptorMap` 会核对原子数和有序结构哈希。

`colorBy` 控制着色模式：

- `"auto"`（默认）：包含片段且未指定 `bondField` 时按片段着色；否则按数值连续着色（`"value"`）。
- `"value"`：带正负号字段（如 `q_e`、`cm5_q_e`、`dual`、`pi_q_e`、`mul_ovlp_e`）使用零对称
  红/蓝发散色阶，非负字段使用单向蓝色色阶，并默认附带底部色条（可用 `showColorbar` 控制）。
- `"fragment"`：按片段分配离散颜色。
- `"none"`：不填充高亮颜色。

## 片段汇总与组合 HTML 视图

```python
from eled import (
    CalcDescriptorReport,
    DrawDescriptorMap,
    PlotFragmentSummary,
    PlotFragmentInteractions,
    RenderFragmentAnalysis,
)

fragments = {
    "acetyl": [0, 1, 2],
    "phenoxy": [3, 4, 5, 6, 7, 8, 9],
}
report = CalcDescriptorReport(mol, result, fragments=fragments)

molecule_svg = DrawDescriptorMap(mol, report, atomField="q_e")
summary_svg = PlotFragmentSummary(report, field="q_e")
network_svg = PlotFragmentInteractions(report)
html_fragment = RenderFragmentAnalysis(mol, report, atomField="q_e")
```

- `PlotFragmentSummary` 绘制任意片段数值字段（如 `q_e`、`cm5_q_e`、`f_plus`、`dual`、
  `s_loc_1ev`）的双向条形图，动态适配标签宽度，并将缺失值显式标为 `missing`。
- `PlotFragmentInteractions` 使用 `report["fragments"]["pairs"]` 中的片段对 Wiberg 总和，
  以自适应宽度的胶囊节点绘制相互作用网络，不会重复叠加双向 external sum。
- `RenderFragmentAnalysis` 组合分子图与片段视图（无片段时自动仅展示分子图），默认将
  `summaryField` 与 `atomField` 联动，并提供与内嵌 SVG 一致的亮/暗色主题样式。

核心包不增加 matplotlib、Plotly 或 NetworkX 依赖，所有函数返回自包含的 SVG/HTML
字符串，可直接嵌入 notebook、报告或网页。

片段必须由用户显式提供。未分配原子仍会显示在分子图中，并由 `BuildFragmentVisualData`
返回。目前不会绘制实空间电子密度或轨道等值面，因为 xTB 结果包含 AO 矩阵，而不是体网格。
