# 进程隔离的并行 libxtb 批处理（0.11）

`RunXTBBatch`、`RunDescriptorBatch` 和 `RunReactivityBatch` 接受
`maxWorkers=1`、`timeout=None` 和 `isolate=False`。默认值保持原有的进程内行为。
当工作进程数大于 1、设置了超时或 `isolate=True` 时，会启动全新的 spawn 方式
Python 进程，每个进程直接调用同一个 libxtb 库。不会搜索或调用 xTB 可执行程序。

```python
from eled import RunDescriptorBatch, XTBOptions

def main():
    # 分子必须带显式氢原子并已准备好 3D 坐标。
    reports = RunDescriptorBatch(
        molecules, options=XTBOptions(keepOverlap=True),
        maxWorkers=4, timeout=60.0,
    )

if __name__ == '__main__':
    main()
```

按照 Python spawn 模式的要求，请使用可导入的脚本并加上 `__main__` 保护。无法导入
主模块的交互环境应改为运行脚本。启动 Python 之前请设置 `OMP_NUM_THREADS=1` 和
`OPENBLAS_NUM_THREADS=1`，避免各工作进程叠加数值库线程。适配器不会在调用方不知情时
修改全局线程设置。每个工作进程从新的解释器启动，处理完一个分子（反应性计算为一组
完整的三个电荷态）后退出。这种设计优先保证隔离，而不是启动吞吐量；小分子使用默认的
顺序接口可能更快。

## 执行与失败

超时是有限的正秒数，从任务启动工作进程开始计时，涵盖进程启动、计算和结果序列化，
不包括排队时间；它也不是对每个电荷态单独计时。超时的工作进程会先被终止，必要时强制
结束并回收。原生代码崩溃或进程未返回结果即退出会产生 `XTBWorkerError`；超时产生
`XTBTimeoutError`。两者都派生自 `XTBRunError`。

`RunXTBBatch(returnExceptions=True)` 会按输入顺序保留这些异常；否则抛出异常并清理
其他活动的工作进程。描述符批处理生成带空值的显式错误行。反应性批处理记录工作进程
失败行，并保留检查点中的尝试历史。被强制结束的三态计算没有可恢复的部分诊断，重试时
从参考态重新开始。隔离路径中的配置、安装和存储错误会直接中止，不会被转换成看似
SCF 失败的结果。

输入分子会做快照，保留双精度坐标、属性、原子顺序和实际构象 ID。原生指针不会跨越
进程边界。工作进程会把自己持有的结果临时序列化到 Python 临时目录下私有的
`eled-worker-*` 目录。正常退出、超时和已处理的中断都会删除这些目录；如果整个父进程
被强制结束，则可能残留。这些是进程间传输文件，不是 xTB 计算日志。持久化存储仍使用
带校验和的 JSON，而不是 pickle。

## 缓存与检查点

缓存由父进程持有，父进程在各工作进程完成时写入检查点；返回的列表保持输入顺序。只有
成功的结果会进入缓存。同时发生的相同未命中可能会各算一次。缓存命中不会启动工作进程。
数值身份由原生计算设置和输入哈希决定，工作进程数和超时不影响它。缓存中的耗时记录的是
最初那次计算。

反应性检查点复用会保留失败行，除非设置 `retryFailed=True`。因此，只增大超时并不会
重试以前的失败：

```python
reports = RunReactivityBatch(
    jobs, checkpoint='dataset.sqlite', cache=cache,
    maxWorkers=4, timeout=120, retryFailed=True,
)
```

失败、重试和已完成的结果都可以审计。输入和报告仍全部保存在内存中；这是一个本地批处理
执行器，不是流式或分布式调度器，也不保证加速。
