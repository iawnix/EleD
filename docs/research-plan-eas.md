# 可发表研究方案：RDKit 原生 xTB 密度描述符与 EAS 区域选择性

## 研究定位

拟定题目：**RDKit-native, basis-consistent xTB density descriptors for
atom-level electrophilic aromatic substitution regioselectivity**。

这项工作的科学问题不是重新实现 xTB，也不是声称发明了 Mulliken、Löwdin
或 Fukui 定义。核心问题是：

> 在严格按分子骨架划分、避免同系物泄漏的原子级 EAS 任务上，来自同一次
> xTB 计算的 shell-resolved density 信息，是否能在 partial charge 和
> RDKit 结构特征之外提供稳定的区域选择性信息？

这个问题有明确的化学意义：EAS 标签属于原子级，而常用的 QM8 激发能等
分子级标量无法检验原子位点排序。它也有明确的软件意义：现有接口通常只
暴露电荷、能量或轨道标量；本项目把同一次计算的 `P`、`S`、AO atom/shell
map 和 RDKit 原子编号绑定起来，使原子级 density feature 可以复现、检查和
直接进入 RDKit 工作流。

## 可主张的创新点

论文中只主张以下三点，避免把成熟的量子化学定义写成新理论：

1. **接口创新**：为扩展版 `libxtb` 提供与 RDKit 原子顺序一致的 overlap、
   AO/shell 映射和密度分析接口；Mulliken/Löwdin shell populations 在同一
   计算结果上提取，不重新调用量子引擎。
2. **验证创新**：首次在本项目固定的 atom-labelled EAS 协议中，对
   charge、shell、density、frontier 和 Fukui 特征做严格逐组消融，并用
   scaffold 加 fingerprint-cluster 划分测试跨骨架泛化。
3. **可靠性创新**：把构象离散度、GFN1/GFN2 差异和电子温度扰动作为三个
   分开的诊断量，检验“特征有效”是否只是几何或数值协议敏感性的结果。

因此，论文可以在新增特征没有稳定提升时仍然成立：负结果将回答哪些
看似有物理意义的 density feature 在严格外推设置下没有额外价值。

## 预注册假设

- **H1（主要假设）**：GFN2 shell/density 特征相对于 GFN1 partial-charge
  baseline，在 EAS 原子排序上提高至少一个主要指标，且分子级 bootstrap
  置信区间不跨零。
- **H2（来源假设）**：shell-resolved s/p 特征的增益大于只增加分子级
  HOMO/LUMO 标量的增益；这一区分直接检验“原子信息”而不是“更多轨道数值”。
- **H3（稳健性假设）**：使用多个构象的均值和离散度比单一最低能构象更稳健，
  但电子能加权的结果不解释为自由能结果。
- **H4（边界假设）**：Mayer、前线键分解和局部 π 投影的效果依赖化学子集，
  不能在没有独立验证时作为普适特征。

H1--H4 在查看测试集结果前冻结。若 H1 不成立，结论写成“接口和验证协议
完成，但没有观察到稳定增益”，不能改写为方法成功。

## 数据与标签

主要数据使用公开的 RegioML/RegioSQM EAS 原子标签；`protocol.json` 规定
输入筛选、原子映射、排除规则和泄漏约束。每条记录保留原始 ID、来源行、
SMILES、canonical SMILES、骨架和被标记的反应原子。

纳入规则：单一反应分子、可重现的原子映射、明确的 EAS 位点标签、RDKit
能够完成 sanitization 和 3D 坐标准备。混合物、金属体系、重复反应、歧义
标签和无法映射的标签全部排除，并在 `excluded` 中记录原因。CYP 代谢位点
只作为预先声明的备选任务，不与 EAS 结果合并。

有效样本数、元素分布、每个分子候选原子数和正负标签比例在特征计算前
冻结并报告。若有效分子数不足以支持独立 scaffold 测试集，研究标记为
pilot，不能包装成泛化结论。

## 计算协议

使用现有扩展版 `libxtb`，不使用 tblite 或 ASE。RDKit 负责分子图、氢原子、
构象生成和原子编号；xTB 只负责电子结构。

- 构象：沿用 `protocol.json` 的 ETKDGv3 和固定随机种子；失败构象保留原因，
  不用零值填补。
- 几何：每个分子只用预先声明的 GFN2 优化坐标作为比较几何；GFN1 和 GFN2
  在完全相同的坐标上做单点计算，避免“方法差异”和“几何差异”混淆。
- 电子温度：主分析使用 300 K；300 K 与 1000 K 的差异只作为诊断。这里的
  温度是 xTB 的 Fermi 占据平滑参数，不是核坐标的热力学温度，不生成热构象，
  也不代表激发态。300 K 是统一的常规基线；1000 K 是有意放大的稳定性扰动，
  用来暴露小能隙、开壳层或近简并体系对分数占据的敏感性。
- 密度：所有需要 shell 或 population 的计算使用 `keep_overlap=True`，并
  检查 `trace(P S)`、shell-to-atom 加和和 ctypes/compiled 一致性。
- 电荷态：Fukui/局部概念 DFT 使用同一几何的 Q、Q+1、Q-1 GFN1 或 GFN2
  结果；不做几何弛豫和自动自旋态搜索。GFN0 不参与 ΔSCF 结论。

当前 `protocol.json` 已冻结主要 EAS 规则；在真正生成特征前还需要把“配对
的 GFN1/GFN2 计算、模型种子、缺失值处理和统计重复次数”写入协议哈希，
作为论文的 protocol artifact。

## 特征组与公平比较

每个原子保留 RDKit 原子编号。构象级数值按 atom index 聚合，模型中只使用
训练折拟合的标准化和缺失值处理。

| 组别 | 内容 | 目的 |
| --- | --- | --- |
| A | RDKit 原子/键、芳香性、环、价态、局部拓扑 | 结构基线 |
| B | A + GFN1 partial charge、Wiberg 和 shell baseline | 现有 xTB 基线 |
| C | B 中同样的 GFN2 特征 | 方法替换效应 |
| D | C + Mulliken/Löwdin s/p/d、p/s、密度守恒诊断 | 主要候选增量 |
| E | D + Wiberg/Mayer、前线键和 π 特征 | 机制探索组 |
| F | D + 三态 Fukui、局部软度和局部亲电性 | 反应性扩展组 |
| G | D 或 F + 构象 mean/std/min/max | 集成与可靠性效应 |

A--D 是主要结论所需的最小闭环；E--G 必须单独报告，不能把所有字段
拼接后只报一个最高分。探索性特征的缺失保留为显式 mask 和原因，训练时
只在训练折中估计中位数；测试折不能参与填补或特征筛选。

分子级能量、HOMO/LUMO 和偶极只能作为分子上下文，不能复制后伪装成独立
的原子观测。所有原子级结果报告覆盖率和缺失率。

## 模型、划分和指标

主划分将共享 Bemis--Murcko scaffold 或 Morgan Tanimoto 0.7 cluster 的
记录放在同一 partition；训练集内部再做分组交叉验证和超参数选择。测试集
在任何特征选择、阈值选择和模型调参前锁定。

模型使用两个互补的固定实现：

- 正则化 logistic regression：检验线性增量和系数方向；
- 梯度提升树：检验非线性、局部相互作用和缺失指示的增量。

主要指标按分子等权计算，而不是让原子数多的分子支配结果：macro ROC-AUC、
MRR、top-1 和 top-2 命中率。补充报告 PR-AUC、Brier/calibration、每个
分子候选位点数分层结果，以及每个元素和芳香环子集的结果。

每个主要指标使用以分子为单位的至少 2,000 次 bootstrap，报告候选组减去
GFN1 baseline 的差值和 95% CI。预先定义的成功条件是：至少一个主要指标的
差值 CI 完全高于零，且其余主要指标没有显著恶化；否则报告为未观察到稳定
增益。随机标签、只用几何/拓扑的模型和原子编号打乱作为负对照。

## 必须完成的技术验证

在进入模型训练前完成并归档：

- AO 排列、轨道相位、简并子空间旋转、刚体平移/旋转不变性；
- `trace(P S)` 与电子数守恒；shell 布居与原子布居守恒；
- ctypes 与 compiled 两条路径在 GFN0/1/2 上的数值一致性；
- 闭壳层 Mayer 与 Wiberg 并列结果；开壳层明确返回缺失原因；
- 构象之间 atom/bond index 固定，前线轨道组组成变化单独标记；
- 每条记录的引擎 SHA256、参数文件 SHA256、温度、坐标哈希、失败原因和
  特征覆盖率。

这些是方法可信度检查，不作为预测性能的额外样本。

## 结果、风险和发表出口

论文结果至少包含：数据流图、接口/密度守恒验证、特征消融表、跨骨架测试、
构象和方法敏感性、运行成本、缺失率以及所有不利结果。不要把基态 HOMO/LUMO
描述成激发能，也不要把 GFN0 EEQ 电荷描述成 Mulliken 或 CM5 电荷。

如果 D 相对 B/C 有稳定增益，主要论文可投 **Journal of Chemical Information
and Modeling** 或相关计算化学期刊；软件接口、测试和可复现运行包可同时作为
**JOSS** 软件论文。若没有预测增益，仍可发表为“可复现的 RDKit--xTB
atom-level descriptor benchmark”，但必须把负结果和适用边界放在摘要和结论中。

交付物固定为：协议 JSON 及哈希、筛选后数据清单、特征 JSONL/Parquet、失败
清单、训练脚本、锁定环境、完整测试日志、模型和 bootstrap 结果。任何测试集
结果都不能只保存在 notebook 中。

## 执行顺序

1. 在特征计算前补全并哈希 protocol artifact，确认 GFN1/GFN2 配对计算和
   模型/随机种子设置。
2. 用 20--50 个代表性分子完成端到端 dry run，检查映射、覆盖率、守恒和
   运行成本；dry run 不进入最终测试集。
3. 锁定完整数据清单和 scaffold/cluster 划分，生成 A--G 特征组。
4. 先运行技术验证，再运行模型；所有选择只发生在训练折。
5. 生成主表、负结果表、敏感性图和可复现归档，最后才撰写性能结论。
