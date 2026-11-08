# 可选的 YAeHMOP/EHT 对照

`eled` 的核心协议使用项目构建的扩展版 `libxtb`，执行 GFN0、GFN1 和
GFN2 计算。RDKit 的 `rdEHTTools`/YAeHMOP 是另一套已有的电子结构模型。

`eled.descriptors` 中的辅助函数只把 YAeHMOP 的标量结果整理成便于对照
的形式：

```python
from eled.descriptors import CalcEHTDescriptors, CalcOrbitalDescriptors

# CalcEHTDescriptors 消费已有的 rdEHTTools 结果。
# CalcOrbitalDescriptors(..., backend="eht") 调用 RDKit 的 YAeHMOP 接口。
```

YAeHMOP 与 xTB 的哈密顿量、基组约定、电荷模型、能量定义和参考零点都不同。
两者的总能量和布居不能混用，也不能视为同一种描述符。EHT 字段不属于
`DescriptorCatalog()`、xTB 报告 schema、EAS 特征协议或扩展版 libxtb 能力检查。

这个适配器不修改 YAeHMOP，也不会把它作为 xTB 工作流的依赖。只使用项目
xTB 描述符的用户可以完全忽略它。
