# EAS 区域选择性研究方案

## 科学问题

在严格按分子骨架和相似性聚类划分的原子级 EAS 任务中，来自同一次 xTB
计算的 shell-resolved density 特征，是否能在 RDKit 结构特征和 GFN1
partial-charge baseline 之外提供稳定的位点排序信息？

这个问题比 QM8 的分子级激发能回归更适合检验原子描述符，因为标签直接对应
反应原子。项目不声称发明 Mulliken、Löwdin 或 Fukui 定义，贡献集中在：

1. RDKit 原子编号与扩展版 `libxtb` 的 `P`、`S`、AO/shell map 的一致接口；
2. 原子级 EAS 任务上的严格特征消融和跨骨架测试；
3. 构象、方法和电子温度三种可靠性诊断。

## 特征组

- A：RDKit 原子、键和局部拓扑；
- B：A + GFN1 电荷、shell 和 Wiberg 特征；
- C：GFN2 的对应特征；
- D：C + Mulliken/Löwdin s/p/d 和 p/s 特征；
- E：D + Mayer、前线键和 π 特征；
- F：D + Fukui、局部软度和局部亲电性；
- G：D/F + 多构象 mean/std/min/max。

A--D 是主要结论，E--G 单独作为探索性扩展。缺失值保留 mask 和原因，
不能把缺失填成零。

## 计算和评估

这里的模型和评估属于独立的 benchmark 层，不是描述符 API 的一部分。
机器可读的模型协议见
[`benchmarks/eas/model_protocol.json`](../../benchmarks/eas/model_protocol.json)，
`protocol.json` 只负责特征生成和数据清单。

使用公开 RegioML/RegioSQM EAS 原子标签，沿用 `protocol.json` 的筛选、
严格 ETKDGv3（不启用其它 preset 的静默回退）、GFN2 几何和扩展版 `libxtb`。
GFN1 与 GFN2 在同一坐标上做单点。主分析为 300 K；同时记录 GFN2/1000 K
电子占据敏感性诊断，不代表热构象或激发态。输入清单保存协议 SHA256 和
泄漏安全的五折 `cv_fold`。

模型使用正则化 logistic regression 和梯度提升树。测试集按 scaffold 加
Morgan cluster 隔离，所有标准化、特征选择和缺失处理只在训练折拟合。
主要指标为分子等权的 macro ROC-AUC、MRR、top-1 和 top-2，使用分子级
2,000 次 bootstrap 计算候选组相对 GFN1 baseline 的差值和 95% CI。

如果没有稳定增益，负结果仍然保留，并明确描述适用边界。完整预注册方案见
[English research plan](../research-plan-eas-en.md)。
