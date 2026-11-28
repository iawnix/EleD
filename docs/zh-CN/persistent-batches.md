# 可持久化的反应性数据集（0.8）

`RunReactivityBatch` 在固定几何下计算三个电荷态，支持显式记录 ID、逐分子的自旋/选项/
片段设置、成功结果缓存和可断点续算的检查点。它沿用现有的 libxtb 协议（默认 300 K，
并从参考态复制波函数）。不做结构优化，也不搜索自旋态。

```python
from eled import (
    ReactivityJob, XTBOptions, XTBDiskCache, ReactivityCheckpoint,
    RunReactivityBatch, ExportDescriptorReports,
)

# 提供带显式氢原子和有限 3D 构象的分子。
jobs = [
    ReactivityJob("water", water),
    ReactivityJob("radical", methyl_radical, ionizedUnpaired=0, attachedUnpaired=0),
    ReactivityJob("custom", another_mol, options=XTBOptions(accuracy=0.1),
                  fragments={"first_atom": [0]}),
]
cache = XTBDiskCache("results/cache.sqlite")
checkpoint = ReactivityCheckpoint("results/dataset.sqlite")
rows = RunReactivityBatch(jobs, cache=cache, checkpoint=checkpoint)
ExportDescriptorReports(rows, "results/descriptors.csv")
ExportDescriptorReports(rows, "results/descriptors.jsonl")
```

中断后，用相同的任务和检查点再次调用同一函数即可继续。ID 必须是非空且唯一的字符串；
输出顺序与输入顺序一致。任务自己的 `options` 会整体替换批处理的 `options` 对象，而不是
逐字段合并。`populationScheme`、`interface` 和 `initialization` 是整个批次的设置，含义与
`RunReactivity` 相同。开壳层任务必须同时指定相邻两个电荷态的未成对电子数。每个输入分子
都会做快照。

## 缓存与检查点

**缓存只保存成功的结果**，可以在不同数据集之间共享。记录 ID 不参与缓存身份，因此同一个
计算可以服务于不同的 ID。只缓存完整的反应性三态计算；失败的报告和不完整的波函数不会
缓存。`XTBDiskCache` 也可用于 `RunXTBBatch` 和 `RunDescriptorBatch`，并在规范化的单点
结果中保留元组/矩阵结构。不使用 pickle、原生句柄，也不引入新的 Python 依赖。

**检查点记录完整有序的数据集以及每次已提交的尝试**，包括失败。完整的三态计算成功或报告
计算失败后，该行以原子方式提交。三态计算中途中断时，该分子从参考态重新开始。如果中断前
已成功写入缓存，即使对应的检查点行尚未提交，这次计算仍可复用。

已完成的检查点行默认会被复用，包括失败行。使用 `retryFailed=True` 可以在相同协议下显式
重试失败项。`checkpoint.GetHistory(recordId)` 返回所有已保存的尝试；用成功的重试替换
当前失败行时，之前的诊断信息会保留。`checkpoint.GetReports()` 按数据集顺序返回已完成的
行。没有产生完整报告的中断尝试不计入已提交的尝试。

更改 ID、顺序、几何或任何计算协议都会导致现有检查点被拒绝。新增数据、修正输入或使用新
设置时，请使用**新的检查点路径**并保留共享缓存；未改变的成功计算仍可复用。缓存与检查点
必须使用不同的文件。

## 协议身份

缓存键和检查点清单会对以下内容做指纹：

- 带索引的拓扑、几何、构象、分子电荷和参考态自旋；
- GFN 方法、迭代上限、精度、电子温度和基组保留设置；
- 相邻电荷态的自旋、初始化策略、布居方案和片段索引；
- libxtb 二进制和参数的身份，以及相关的环境变量/线程设置；
- Python 实现源码、使用时的编译适配器二进制、Python 主次版本、RDKit 和 NumPy 版本；
- 存储/缓存、描述符/schema 和反应性协议的版本。

这些检查刻意保守：即使数值结果相同，更改路径/环境或重新构建编译模块也可能使复用失效。
续算需要一个可用且兼容的引擎；缓存文件不能离线替代运行时校验。修改代码或库之后请重启
Python，并在一次运行期间保持环境不变。离线分析请直接打开已导出的报告。

`provenance.batch_protocol` 记录参与指纹的配置。`execution` 包含 `source`
（`calculation`、`cache`、`checkpoint` 或 `input_validation`）、`key` 和已提交的尝试
编号。CSV 除已有的 provenance/settings/error 列外还包含 `execution_json`；JSONL 保留
所有嵌套的电荷态和原子数据。

## 错误策略与存储保证

所有任务配置和运行时能力都在计算前检查。无效的分子几何（包括缺少构象、未展开的氢计数或
非有限坐标）以及 `mol=None` 会生成描述符为空的输入错误行。`None` 可以表示上游解析失败；
调用方应另外保留原始记录。真正的 xTB 单点失败会保留三个电荷态的诊断信息，派生描述符
为空。

缺少开壳层自旋、自旋奇偶性无效、片段定义错误、ID 重复、方法不受支持、库或基组能力缺失
/不兼容，都会直接抛出异常。意外的编程或分析错误同样会抛出。从 0.8 开始，
`RunReactivity(returnErrors=True)` 不再吞掉任意异常：分配/初始化/提取失败抛出
`XTBEnvironmentError`，真正的单点失败仍以 `ChargeStateCalculationError` 报告。C++ 结果
通过 `failureKind` 区分这两类；升级后请重新构建编译接口的使用方。

SQLite 提交保证已完成的行不受中断影响。JSON 内容带 SHA256 校验和；数据损坏、存储
schema 不兼容和 I/O 失败会直接中止，而不会悄悄变成缓存未命中或分子计算失败。校验和用于
发现意外的记录损坏。

检查点运行时持有一个非阻塞的操作系统文件锁：另一个写入同一检查点的进程会明确失败。共享
缓存的写入使用 SQLite 事务；同时发生的相同未命中可能各算一次。默认按顺序执行；
[隔离工作进程](isolated-batches.md) 支持并行任务和墙钟时间限制。缓存和检查点的写入始终
在父进程中、在本地文件系统上进行；这不是分布式调度器。磁盘容量、保留和清理策略由调用方
管理。批处理会把输入和返回的报告全部放在内存中；大数据集可以拆成多块，各自使用独立的
检查点并共享同一个缓存。

## 复现示例与验证

```bash
python examples/reactivity_batch.py --output build/my-dataset
python examples/reactivity_batch.py --output build/my-dataset  # 续算
python tools/validate_reactivity_dataset.py --output build/corpus-check
python tools/compare_fukui_cli.py \
  --executable .eled/runtime/bin/xtb \
  --library .eled/runtime/lib/libxtb.so \
  --cases tests/data/reactivity_cases.json --output build/cli-check
```

冻结的测试集包含 37 个分子，覆盖常见官能团、阳离子/阴离子、卤素、含磷/含硫化合物以及
四个开壳层体系。数据集验证脚本为每种方法额外加入两个低温算例和一个迭代上限算例，然后
检查 GFN1/2、两种初始化方式和两种接口（共 320 组三态计算），包括不重新计算的续算。原有
的命令行对比覆盖 74 个标准分子/方法组合，受限于其 Fukui 输出只有三位小数。所有失败都计入
数据集统计。未优化的几何和接口之间的一致性，不能说明实验精度或激发态的有效性。
