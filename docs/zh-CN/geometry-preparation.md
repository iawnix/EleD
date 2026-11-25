# 构象准备

`PrepareConformers()` 是 RDKit 构象生成和检查层，负责从 SMILES 或分子对象
生成显式氢、3D 坐标并记录失败原因。它不是新的嵌入算法。

```python
from eled import PrepareConformers

prepared = PrepareConformers(
    "C[C@H](O)C(=O)O",
    numConformers=4,
    method="ETKDGv3",
    optimize="MMFF94",
)
if not prepared.success:
    raise RuntimeError(prepared.report)
mol = prepared.molecule
```

`optimize="UFF"` 使用 UFF；`optimize="xTB"` 使用扩展版 xTB 梯度优化，成本
更高。几何优化的作用是为电子描述符提供可用坐标，不是全局构象搜索、过渡态
搜索或溶剂化计算。

`OptimizeXTB()` 在线搜索的试探几何上遇到 `XTBRunError` 时，将位移减半后重试。
若全部 14 次试探计算都失败，抛出最后一个引擎错误；若至少一次计算完成但没有
降低能量，则返回 `reason="line_search_failed"`，保留最后接受的几何和计算结果。
初始几何的计算失败直接抛出。

报告记录实际使用的种子、嵌入尝试、立体化学检查、优化状态、坐标哈希和
推荐构象。构象失败会保留原因，不会静默复制或用零值代替。
