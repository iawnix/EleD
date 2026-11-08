# 构象准备

`prepare_conformers()` 是 RDKit 构象生成和检查层，负责从 SMILES 或分子对象
生成显式氢、3D 坐标并记录失败原因。它不是新的嵌入算法。

```python
from eled import prepare_conformers

prepared = prepare_conformers(
    "C[C@H](O)C(=O)O",
    num_conformers=4,
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

报告记录实际使用的种子、嵌入尝试、立体化学检查、优化状态、坐标哈希和
推荐构象。构象失败会保留原因，不会静默复制或用零值代替。
