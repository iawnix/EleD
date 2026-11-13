# 应用案例调研

本页记录 EleD 的候选应用案例和选择依据。它是研究决策记录，不是模型
基准测试。数据集、模型代码、训练权重和数值结果仍放在外部研究工程中。

## 案例需要证明什么

合适的案例应同时满足以下条件：

1. 标签附着在原子或键上，RDKit 的原子身份是问题的一部分；
2. 能记录数据源、版本和数据使用权；
3. 能按分子划分，避免近邻类似物跨越测试边界；
4. 同一份 EleD 归档可以交给小型表格模型和图模型，而不必重新计算电子结构。

案例要证明的是电子描述符接口可用。它不把 EleD 变成反应预测器，也不要求
某个 EleD 特征组超过已发表模型。

## 2026-11-13 核查的候选案例

| 案例 | 标签和规模 | 数据源与权利 | 与 EleD 的匹配度 | 决定 |
| --- | --- | --- | --- | --- |
| **RegioSQM20** | EAS 反应位点标签；随附 `compounds.smiles` 有 535 个分子、604 个原始标注原子条目 | 代码和随附仓库带 MIT 许可证；论文为 CC BY 4.0。导出反应数据的再分发权仍需单独记录 | 很高：标签粒度与电荷、Fukui、shell、π 和前线描述符一致；输入保留 RDKit 原子顺序 | **第一案例** |
| **RegioML** | EAS 原子分数；公开模型和描述符代码，完整训练导出在 ERDA | GitHub 代码为 MIT；ERDA 数据包及再分发条件需单独核对 | 很高，适合作为相关工作和互操作性检查；外部数据不适合做第一个独立数据源 | **参考 / 可选复核** |
| **GNN-SOM** | KEGG RPAIR 的原子 SOM 和键 SOM 标签；论文报告 21,023 个分子 | RPAIR 按 KEGG 学术许可下载；GitHub 仓库没有整理后的训练数据 | 可检验 EleD 的节点和边归档，但许可和反应中心重建增加风险 | **第二阶段** |
| **ESNUEL / ESNUEL-ML** | QM 得到的甲基阳离子亲和能和甲基阴离子亲和能；分别有 47,921 和 47,440 个分子，原子级回归目标 | 代码为 MIT；数据和原始 QM 归档通过 ERDA 分发，使用前需核对 | 适合大规模数值回归和归档吞吐量研究；与实验选择性联系较弱 | **可选规模案例** |
| **SET-LSF** | Minisci/P450 后期官能化标签；公开文献子集和受限的完整数据集 | GitHub 代码和非专有子集为 MIT；完整 Pfizer 数据集受限 | 有化学意义且是原子标签，但机制混合、依赖反应条件，公开子集也较小 | **未来可选** |

## 为什么先用 RegioSQM20

RegioSQM20 论文描述了 535 个文献 EAS 溴化反应，并给出实验观察位点的
原子级标签。此次核查的支持数据仓库快照为
`jensengroup/SI_RegioSQM20@6f9539594b9a07babd76db347ca6d2f550c3ad97`。
其中 `RegioSQM20/compounds.smiles` 每行是
`record_id SMILES site_indices`；逗号分隔的编号表示多个观察位点。使用
RDKit 2026.09.1 解析得到 535 条有效记录、604 个原始正标签原子，没有越界
编号。随附 SDF 是二维结构文件，不包含这些反应标签，因此适配器必须把
SMILES 文件和行身份一起保留。

这个案例适合 EleD，因为目标和报告特征处于同一原子粒度。特征行可以同时
保留原始 RDKit 原子编号、标签、缺失原因和计算来源；同一份冻结归档随后可
交给原子排序表格模型或节点/边 GNN。论文使用了 xTB/GFN1、质子亲和能、互变
异构体枚举和溶剂模型。这些属于相关工作或可选消融，不要求 EleD 复刻
RegioSQM20 预测器。

计算前需要明确处理原始导出：

- 把逗号分隔的标签保留为原子编号集合；
- SMILES 清洗后验证标签仍是芳香碳原子；
- 多个位点记录保持多标签语义，不压成单一类别；
- 记录来源提交、文件哈希和标签解释规则；
- 将仓库代码许可证与底层反应观测数据的权利分开记录。

RegioSQM20 也被后续反应选择性和量子/机器学习工作作为参考。例如 2024 年
的 ESNUEL 论文以及多项近期位点选择性研究都引用了它。这个案例因此容易
被研究者理解，同时可以清楚呈现 EleD 的贡献：可复用、保留 RDKit 原子/键
身份的电子特征归档，而不是又一个 RegioSQM 实现。

## 其他候选案例的作用

**RegioML** 是最接近的历史对照。它用 GFN1-xTB CM5 电荷和排序后的
GraphChargeShell 向量，再交给 LightGBM。EleD 可以通过公开接口产生同类
输入，但不复制其模型或解析器。完整数据包链接在 ERDA，因此在确认数据条款
和原子映射前，不把它作为第一个冻结数据源。

**GNN-SOM** 把目标扩展到原子和键。论文从 KEGG RPAIR 的反应中心标签构建
数据，报告 21,023 个分子和十组 80/10/10 划分。此次核查的 GitHub 快照为
`HassounLab/GNN-SOM@dbe88c5aea3cb45f416d4d11d8bc882b72ce3c57`，其中有推理
notebook、模型状态和配置，没有整理后的训练集。KEGG 的学术许可要求在受控
研究环境中生成特征归档，能否发布归档要按许可条款决定。

**ESNUEL** 适合做规模和目标质量研究。ESNUEL-ML 报告了 47,921 个分子上
650,857 个原子 MCA 值，以及 47,440 个分子上 534,119 个原子 MAA 值；输入
描述符使用 GFN1-xTB CM5，目标来自更高层级 QM。代码公开且为 MIT，但数据和
原始 QM 输出在外部 ERDA；原始计算结果约 29 GB 压缩、355 GB 解压。因此应
在 EAS 归档合约稳定后再做，不作为第一个端到端案例。

**SET-LSF** 可作为自由基反应性扩展。公开仓库有 MIT 许可的非专有文献子集
和 Minisci/P450 位点标签代码，完整 Pfizer 数据集受限。混合反应机制和条件
向量使其不适合作为第一个通用电子描述符 schema 测试，但以后可以检验在外部
研究提供反应条件向量时 EleD 特征是否仍然有用。

## 第一个案例的冻结计划

外部研究工程在任何 xTB 计算前冻结：

1. RegioSQM20 来源提交和 `compounds.smiles` 的 SHA256；
2. 逗号分隔、零基原子编号的解析规则；
3. 分子筛选和重复记录规则；
4. ETKDGv3/构象和扩展版 libxtb GFN2/300 K 计算设置；
5. `protocol.json` 中的特征组和缺失值语义；
6. scaffold 加 Morgan 聚类的分子划分，并在模型选择前封存测试集。

第一轮只用 20–50 个分子做 dry run，报告标签映射、描述符覆盖率、失败构象、
运行时间和归档大小。冻结 manifest 后再计算完整归档，并交给外部模型工程。
外部工程可以训练紧凑的排序/分类模型和小型 GNN，但 EleD 不包含这些实现。

案例报告应分别给出结构特征、电荷特征和完整 EleD 特征组，同时报告覆盖率和
成本。失败和不利结果必须保留；与 GFN1、RegioML 或 RegioSQM20 的比较只能
作为背景或消融，不能成为 EleD 的验收条件。

## 来源记录

- RegioSQM20 论文：[Ree, Göller and Jensen (2021)](https://doi.org/10.1186/s13321-021-00490-7)，文章为 CC BY 4.0。
- RegioSQM20 代码：[jensengroup/RegioSQM20](https://github.com/jensengroup/RegioSQM20)。
- RegioSQM20 数据：[jensengroup/SI_RegioSQM20](https://github.com/jensengroup/SI_RegioSQM20)。
- RegioML 论文：[Ree, Göller and Jensen (2022)](https://doi.org/10.1039/D1DD00032B)。
- RegioML 代码：[jensengroup/RegioML](https://github.com/jensengroup/RegioML)。
- GNN-SOM 论文：[Porokhin et al. (2023)](https://doi.org/10.1093/bioinformatics/btad089)。
- GNN-SOM 代码：[HassounLab/GNN-SOM](https://github.com/HassounLab/GNN-SOM)。
- ESNUEL 工作流：[jensengroup/ESNUEL](https://github.com/jensengroup/ESNUEL)。
- ESNUEL-ML：[jensengroup/ESNUEL_ML](https://github.com/jensengroup/ESNUEL_ML)。
- SET-LSF：[emmaking-smith/SET_LSF_CODE](https://github.com/emmaking-smith/SET_LSF_CODE)。
