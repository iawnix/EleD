# EleD documentation / EleD 文档

Project overview / 项目首页：[English](../README.md) · [中文](../README.zh-CN.md)

EleD is an RDKit extension for electronic-structure descriptors backed by xTB.
The documentation follows the same path as a typical analysis: install the
runtime, prepare a molecule, calculate descriptors, build a report and inspect
or export the result.

EleD 是基于 xTB 的 RDKit 电子结构描述符扩展。文档按实际使用流程组织：安装
运行环境、准备分子、计算描述符、生成报告，最后进行可视化或导出。

## User guides / 用户指南

| Stage / 阶段 | English / 英文 | 中文 / Chinese |
| --- | --- | --- |
| Install / 安装 | [Install](install.md) | [安装说明](zh-CN/install.md) |
| Prepare geometry / 准备构象 | [Geometry preparation](geometry-preparation.md) | [构象准备](zh-CN/geometry-preparation.md) |
| Calculate descriptors / 计算描述符 | [Descriptors](descriptors.md) | [描述符定义](zh-CN/descriptors.md) |
| ESP and fitted charges / 静电势与拟合电荷 | [ESP and RESP](esp-api.md) | [静电势与拟合电荷](zh-CN/esp-api.md) |
| Local sterics and rings / 局部位阻与环 | [Local geometry](local-geometry.md) | [局部几何](zh-CN/local-geometry.md) |
| Build reports / 生成报告 | [Reports and reactivity](descriptor-reports.md) | [报告与反应性](zh-CN/descriptor-reports.md) |
| Explore results / 查看结果 | [Visualization](visualization.md) | [可视化](zh-CN/visualization.md) |
| Export and scale / 导出与批量计算 | [Persistent batches](persistent-batches.md) · [Isolated batches](isolated-batches.md) | [持久化批处理](zh-CN/persistent-batches.md) · [隔离批处理](zh-CN/isolated-batches.md) |

## Reference / 参考

- [Development](development.md) · [开发指南](zh-CN/development.md)
- [Descriptor catalog](descriptor-catalog.md) · [描述符目录](zh-CN/descriptor-catalog.md)
- [Conformer statistics](conformer-sensitivity.md) · [构象统计](zh-CN/conformer-sensitivity.md)
- [Density descriptors](density-descriptors.md) · [密度描述符](zh-CN/density-descriptors.md)
- [C++ interface](rdxtb-tools.md) · [C++ 接口](zh-CN/rdxtb-tools.md)
- [Native xTB backend](xtb-patched-build.md) · [原生 xTB 后端构建](zh-CN/xtb-patched-build.md)
- [Runtime bundle](runtime-bundle.md) · [离线运行时包](zh-CN/runtime-bundle.md)
- [Project positioning](project-positioning.md) · [项目定位](zh-CN/project-positioning.md)

## Validation / 验证

Application studies, datasets and model training are maintained outside the
EleD repository.

应用研究、数据集和模型训练在 EleD 仓库之外维护。

- [SCF diagnostics](scf-diagnostics.md) · [SCC 收敛诊断](zh-CN/scf-diagnostics.md)
- [Validation history](validation.md)

### Engineering records (English only) / 工程记录（仅英文）

These pages record past investigations and migrations. They are kept for
traceability, are not user guides, and are not translated.

以下页面记录历史排查与迁移过程，仅为可追溯性保留，不是用户指南，不提供中文版。

- [Migration to 0.2](migration-0.2.md)
- [tblite comparison](tblite-comparison.md)
- [GFN0 native state](gfn0-native-state.md) · [GFN0 re-verification](gfn0-reverification.md)
- [Orbital alignment](orbital-alignment.md)
- [xTB 6.7.1 basis differences](xtb-671-basis-differences.md)
