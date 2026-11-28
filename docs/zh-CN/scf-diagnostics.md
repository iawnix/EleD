# 成功 SCC 计算的收敛诊断（0.11）

经过审查的 libxtb 扩展新增了 `xtb_rdkit_getSCFAPIVersion()` 和
`xtb_rdkit_getSCFDiagnostics(env, results, iterations, values)`。这是项目自有的
ABI 版本 1，与上游 xTB API 及基组导出 ABI 并存。所选二进制文件必须与其构建记录
及审查过的诊断补丁一致；公开计算接口会拒绝官方版本或旧版本的库。

GFN1/2 单点计算成功后，`result.metadata['scf']` 包含：

| 字段 | 含义 |
|---|---|
| `iterations` | 实际 SCC 迭代次数 |
| `energy_change_hartree` | 最后一次 SCC 迭代的能量变化绝对值 |
| `mixer_rms` | 引擎混合变量差值的最终 RMS |
| `energy_threshold_hartree` | 引擎使用的能量收敛阈值 |
| `mixer_threshold` | 引擎使用的混合变量收敛阈值 |
| `criteria_satisfied` | 两个测量值均低于对应的引擎阈值 |

GFN2 的混合变量包含壳层电荷、原子偶极和四极矩，因此其 RMS **不只是原子电荷残差**，
也不给这个混合度量指定元电荷单位。这些量描述的是数值收敛情况，不代表物理精度，
也不证明得到的是能量最低的 SCF 解。

`CalcDescriptorReport` 在 `quality.scf_details` 中给出该记录（0.17 起删除了与之
重复的 `scf_iterations` / `scf_residual` 字段）。如果观测到报告的收敛判据未满足，
会设置 error 级别的质量标记。标量描述符定义和预测特征集不受影响。

编译接口的 `XTBResults` 提供 `scfIterations`、`scfEnergyChangeHartree`、
`scfMixerRMS`、`scfEnergyThresholdHartree` 和 `scfMixerThreshold`。原始 C++ 字段
在数据不可用时为 -1；共享的 Python 提取函数则返回显式的缺失记录。GFN0 报告
`not_applicable_non_scf_method`。

C 接口只在数据可用时写入一个整数和四个双精度数。新结果、失败或重复计算、Hessian
以及 CPCM-X 入口都会清除上一次的快照；结果副本各自持有自己的诊断数据。第一版 ABI
只导出成功 SCC 的诊断。SCF 失败时只保留错误信息和请求的设置，不会编造残差，也没有
可恢复的迭代轨迹。

源码构建脚本会运行一个只依赖标准库的 C-ABI 检查，覆盖可用性、GFN0、成功的 GFN1/2、
复制的结果、失败后复用以及失效清除。集成测试比较 ctypes 与编译接口的数值。更新源码
锁定文件时，诊断补丁在已有基组补丁之后应用，并校验最终源码哈希。需要重新构建原生
引擎和编译适配器才能获得这项功能。
