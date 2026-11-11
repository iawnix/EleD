# 报告与反应性

## 单点报告

```python
from eled import XTBCalculator, CalcDescriptorReport

result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
report = CalcDescriptorReport(mol, result, record_id="molecule-1")
```

`CalcDescriptorReport()` 不会再次调用引擎。报告包含：

- `values`、`units`：固定的分子级字段和单位；
- `atoms`、`bonds`、`fragments`：按 RDKit 原始编号排列；
- `populations`、`frontier`、`density`、`shell_populations` 等嵌套结果；
- `analysis_settings`、`provenance`：阈值、几何、库和参数来源；
- `quality`、`missing`、`error`：质量标志、缺失原因和错误。

没有保留 overlap 时，密度相关区段为 `None` 并记录
`basis_not_retained`。报告不会用零填补缺失值。

## π 分析参数

`CalcDescriptorReport`、`RunDescriptorBatch` 和 `RunConformerAnalysis` 共用：

| 参数 | 默认值 | 含义 |
| --- | --- | --- |
| `pi_scheme` | `"lowdin"` | `"lowdin"` 或 `"mulliken"` 布居/投影定义 |
| `pi_plane_tol_a` | `0.1` | 局部平面最大距离容差，Å；有限且非负 |
| `pi_ref` | `None` | RDKit 原始原子编号 → 参考 π 电子数，范围 0 至 2 e |

```python
pi_settings = dict(pi_scheme="mulliken", pi_plane_tol_a=0.05, pi_ref={0: 1.0})
report = CalcDescriptorReport(mol, result, **pi_settings)
reports = RunDescriptorBatch([mol], **pi_settings)
```

这些参数写入 `analysis_settings`，失败的批处理行也保留设置。`pi_ref` 只覆盖
列出的原子，其他原子遵循[默认参考规则](density-descriptors.md#局部-π-投影)。
输入支持整数索引和标准十进制字符串索引，报告复制映射并统一保存为字符串键，
从 JSON 读回后可通过 `CalcDescriptorReport(..., **settings)` 重用。

同一批次共享参考索引映射；不同分子需要不同参考定义时，分别提取报告。
无效参数和任一有效输入分子的越界索引会在批处理/构象计算前被拒绝。
这些是分析参数，不改变底层 xTB 计算缓存键；敏感性比较要求 π 参数一致。

## 批量计算和导出

```python
from eled import RunDescriptorBatch, ExportDescriptorReports

reports = RunDescriptorBatch(
    molecules, record_ids=ids, max_workers=4, isolate=True,
)
ExportDescriptorReports(reports, "features.jsonl")
ExportDescriptorReports(reports, "features.csv")
```

批处理保持输入顺序，并保留单分子失败行。JSONL 保存完整嵌套结构；CSV 保存
固定的标量列以及 JSON 编码的来源、质量和缺失信息。不同 schema 或描述符
版本的报告不能混合导出。

## 固定几何反应性

```python
from eled import RunReactivity, XTBOptions

report = RunReactivity(
    mol,
    options=XTBOptions(method="GFN2-xTB", electronic_temperature=300.0),
    population_scheme="charges",
    return_errors=True,
)
```

反应性报告的能量定义为：

- `ip_v_ev = E(Q+1) - E(Q)`；
- `ea_v_ev = E(Q) - E(Q-1)`；
- `eta_ev = (IP - EA) / 2`；
- `s_1ev = 1 / (IP - EA)`。

三种电荷态使用相同几何、方法、参数、电子温度和引擎。300 K 是电子占据
平滑的计算设置，不是分子热力学温度。GFN0 不参与 ΔSCF 反应性分析。

## 报告版本

当前单点报告使用 schema 2、descriptor version 1.3。报告字段带有显式单位，
原子编号固定遵循输入分子。修改公开字段时需要同步更新目录、测试和版本说明。

英文详细说明见 [descriptor-reports.md](../descriptor-reports.md)。
