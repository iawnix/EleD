# 密度描述符

密度描述符只消费已有的扩展版 `libxtb` 结果，不调用引擎。运行时必须使用：

```python
from eled import XTBCalculator

result = XTBCalculator(method="GFN2-xTB").calculate(mol, keep_overlap=True)
```

## 密度矩阵

在非正交 AO 基组中构造：

```text
P = C · diag(n) · Cᵀ
```

其中 `C` 是同一次计算的 AO-MO 系数，`n` 是轨道占据，`S` 是同一次计算的
overlap。电子数由 `trace(P S)` 检查。系数矩阵、overlap 和 AO/shell 映射
必须维度一致且通过正交归一检查。

## Mulliken 与 Löwdin

Mulliken 布居使用 `C` 与 `S C` 的乘积；Löwdin 布居使用 `S^(1/2) C` 的
平方。两者都是基组相关的原子分区定义。程序同时输出 s、p、d shell 电子数、
p/s 比和守恒误差。GFN0 的 EEQ 电荷与这些密度布居保持独立。

## 键和 π 特征

闭壳层可计算 Mayer/Wiberg 并列键级。HOMO/LUMO 前线键分解使用简并轨道组
投影，保留成键和反键的符号。π 特征把局部 p shell 投影到原子或环平面法向；
法向不明确时返回 `None` 和原因。

这些结果是基态密度特征。它们不是激发态跃迁密度、振子强度或电子-空穴
激发态距离。

## 验证

测试覆盖电子数守恒、shell-to-atom 守恒、AO 排列、轨道相位、简并旋转、刚体
平移/旋转，以及 ctypes 与 compiled 结果一致性。完整英文说明见
[density-descriptors.md](../density-descriptors.md)。
