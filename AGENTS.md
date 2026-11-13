# EleD 工作规则（草稿）

> 这是供审查的简短草稿，确认后再作为正式规则使用。

## 项目边界

- 统一在 `/home/iaw/project/EleD` 工作。
- EleD 包含 RDKit 电子结构描述符、扩展版 `libxtb`、`rdXTBTools`、报告、协议、示例和核心测试。
- 数据集、模型训练代码、GNN 代码、训练模型和 benchmark 结果放在 EleD 仓库之外。
- `research/eas_model/` 是本地外部研究目录，被 Git 忽略，不得由 `eled` 导入，也不进入核心 CI。
- 不提交运行时、构建目录、wheel、缓存和大规模生成报告。

## Conda 环境规则

- 只使用固定路径，不按任务、日期或随机名称创建环境。
- 创建环境前先检查固定路径是否已有 `conda-meta/history`；存在时直接复用。
- 锁定文件不匹配时不要创建第二个环境。先报告差异，只有明确要求时才清理并重建。
- 纯 Python 环境固定为 `.runtime-env`。
- 完整 EleD 环境固定为 `.eled/runtime`，对应构建工具链固定为 `.eled/toolchain`。
- C++/xTB 源码构建只使用 `.rdxtb-build-env`、`.xtb-build-env` 和 `.xtb-rdkit`；不创建带版本号或日期的新路径。
- 外部模型研究与 EleD 共用同一个基础运行时版本，但使用独立的固定 overlay 环境；不要把模型依赖安装进 EleD 环境。
- 默认研究环境为 `research/.eas-env`，通过 `.eled/runtime/bin/python -m venv --system-site-packages` 创建，只额外安装研究依赖。
- 如果研究依赖与 EleD 版本冲突，才创建独立环境，并记录 EleD 版本、归档哈希和研究依赖锁定文件。
- 优先使用 `conda run -p <固定路径>`，避免反复激活和创建环境。
- 测试完成后不要自动删除环境；清理必须由用户明确要求。

## 测试环境

- 纯 Python、schema 和报告测试使用固定的 `.runtime-env`。仅在该路径不存在时创建：

  ```bash
  test -f .runtime-env/conda-meta/history || \
    conda create --yes --prefix "$PWD/.runtime-env" \
      --file locks/linux-64/runtime.txt
  conda run -p "$PWD/.runtime-env" python -m pytest -q \
    -m 'not integration and not benchmark'
  ```

- 涉及 xTB、密度描述符、几何优化、native/compiled 或 `rdXTBTools` 时，复用 `.eled/runtime`：

  ```bash
  test -f .eled/runtime/conda-meta/history || \
    python tools/install_native.py --jobs 6
  .eled/runtime/bin/python -m pytest -q -m 'not benchmark'
  ```

- 先运行相关的定向测试，再运行完整测试。C++ 改动还要运行对应的 CTest 和 Python 一致性测试。
- 如果系统 Python 没有 RDKit 或编译运行时，不要用它判断项目测试失败。

- 外部研究使用固定的 `research/.eas-env`，不重复安装 EleD 基础依赖：

  ```bash
  test -x research/.eas-env/bin/python || \
    .eled/runtime/bin/python -m venv --system-site-packages research/.eas-env
  research/.eas-env/bin/python -m pip install -r research/eas-requirements.txt
  research/.eas-env/bin/python -m pytest -q research/eas_model/tests
  ```

## 开发规则

- 保持 RDKit 原有的原子和键编号顺序。
- 每个描述符都要写明单位、定义和缺失值含义。
- 未定义值不能编码成零。
- 保持 native/compiled 结果一致，并检查数值不变性。
- 中英文文档同步更新。
- 提交前运行定向测试和 `git diff --check`。

## 外部研究规则

- 模型研究只消费冻结的 EleD 特征归档，并在独立研究工程中维护。
- 计算或训练前冻结数据源、协议、分子级划分和特征 schema。
- 保留失败记录、覆盖率、来源信息以及原始原子/键映射。
- 超参数只能在训练数据内部选择，测试集在最终评估前保持未使用。
- 报告负结果和限制。案例研究不能改变 EleD API，也不要求 EleD 超过某个基线。

## 外部研究版本与记录

- 外部研究工程使用独立 Git 仓库；每次正式运行前提交代码并保持工作树干净。
- 每个运行都记录 EleD commit/tag、特征协议 ID 和 SHA256、数据源版本和许可证、模型代码 commit、环境锁定文件、随机种子和完整命令。
- 描述符或特征 schema 改变时，增加特征协议版本；模型或超参数改变时，增加模型协议版本。
- 运行结果使用唯一的 run ID 保存，不覆盖已有结果；大数据和模型文件放在外部归档，并保存 SHA256。
- 实验记录至少包含输入 manifest、配置、失败与覆盖率、运行时间、预测结果、指标和备注。
- 发布结果时同时提供代码 commit、特征归档哈希、环境信息和运行记录，保证结果可以追溯到具体版本。
