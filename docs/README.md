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
| Build reports / 生成报告 | [Reports and reactivity](descriptor-reports.md) | [报告与反应性](zh-CN/descriptor-reports.md) |
| Explore results / 查看结果 | [Visualization](visualization.md) | [可视化](zh-CN/visualization.md) |
| Export and scale / 导出与批量计算 | [Persistent batches](persistent-batches.md) · [Isolated batches](isolated-batches.md) | — |

## Reference / 参考

- [Descriptor catalog](descriptor-catalog.md) · [描述符目录](zh-CN/descriptor-catalog.md)
- [Density descriptors](density-descriptors.md) · [密度描述符](zh-CN/density-descriptors.md)
- [C++ interface](rdxtb-tools.md)
- [Native xTB backend](xtb-patched-build.md)
- [Runtime bundle](runtime-bundle.md)
- [YAeHMOP/EHT comparison](eht-comparison.md) · [YAeHMOP/EHT 对照](zh-CN/eht-comparison.md)
- [Project positioning](project-positioning.md) · [项目定位](zh-CN/project-positioning.md)

## Research and validation / 研究与验证

The research pages document the EAS feature protocol, benchmark preparation and
numerical validation of the descriptor layers. They are useful when reproducing
experiments or extending the feature set.

研究页面记录 EAS 特征协议、数据准备和描述符层的数值验证，适合复现实验或
扩展描述符功能。

- [EAS research plan](research-plan-eas-en.md) · [EAS 研究方案](zh-CN/research-plan-eas.md)
- [EAS protocol audit](protocol-audit.md)
- [EAS benchmark](../benchmarks/eas/README.md)
- [Validation scope](mainline-validation.md)
- [SCF diagnostics](scf-diagnostics.md)
- [Literature and ML positioning](literature-density-ml-audit.md)
