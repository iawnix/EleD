# 可视化

EleD 的可视化层消费已有的 RDKit 分子和描述符报告，不会重新启动 xTB
计算。RDKit 提供分子图、原子编号和二维绘图；EleD 提供电子描述符数值和
片段映射。

```python
from eled import (
    DrawDescriptorMap, PlotFragmentSummary, PlotFragmentInteractions,
    RenderFragmentAnalysis,
)

fragments = {
    "acetyl": [0, 1, 2],
    "phenoxy": [3, 4, 5, 6, 7, 8, 9],
}
report = calc_descriptor_report(mol, result, fragments=fragments)

molecule_svg = DrawDescriptorMap(mol, report, atom_field="partial_charge_e")
summary_svg = PlotFragmentSummary(report, field="charge_e")
network_svg = PlotFragmentInteractions(report)
html_fragment = RenderFragmentAnalysis(mol, report)
```

`DrawDescriptorMap` 使用 RDKit 的 `rdMolDraw2D`，保留原始 RDKit 原子编号。
`PlotFragmentSummary` 绘制带正负号的片段数值，例如电荷、Fukui 指数或局部
软度。`PlotFragmentInteractions` 使用
`report["fragments"]["pairs"]` 中明确给出的片段对 Wiberg 值，不会再次把两个
方向的 external sum 相加。

CamelCase 函数是文档中的 RDKit 风格入口；同时提供
`draw_descriptor_map`、`plot_fragment_summary` 等小写别名。核心包不增加
matplotlib、Plotly 或 NetworkX 依赖，函数返回 SVG/HTML 字符串，可嵌入
notebook、报告或网页。

片段必须由用户显式提供。未分配原子仍会显示在分子图中，并由
`BuildFragmentVisualData` 返回。目前不会绘制真正的电子密度或轨道等值面，
因为 xTB 结果包含 AO 矩阵，而不是实空间体网格。
