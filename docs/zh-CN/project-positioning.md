# 项目定位

`rdkit-eth` 是 RDKit 的电子结构能力扩展层。RDKit 继续作为分子图、原子和
键编号、构象、立体化学以及坐标单位的唯一来源；`rdkit-eth` 在此基础上增加
电子结构计算、描述符分析和来源记录。

最接近的参考项目是 DScribe 一类的描述符库。ASE 的 calculator 设计可以借鉴，
但 `rdkit-eth` 不会发展成通用原子模拟环境。周期体系、分子动力学、轨迹管理
和完整溶剂工作流不属于核心范围。

## 扩展边界

```text
RDKit Mol / Conformer
        │
        ├── 几何准备和检查
        │
        ├── XTBCalculator → XTBResult
        │                         │
        │                         ├── 分子/原子/键描述符
        │                         ├── 密度和 shell 描述符
        │                         ├── 固定几何反应性描述符
        │                         └── RDKit 属性标注和报告
        │
        └── 可选 YAeHMOP 对照适配器（不属于核心协议）
```

公开层次保持独立：

1. **RDKit 边界层**：检查 `ROMol`、构象编号、原子顺序和来源；
2. **计算器层**：调用后端，管理数值参数和引擎身份；
3. **结果层**：保存方法相关数组、单位和同一次计算的 basis 数据；
4. **描述符层**：只消费结果，不再次调用引擎；
5. **报告层**：保存单位、缺失原因、质量标志和哈希。

YAeHMOP 适配器只整理 RDKit 已有接口，供结果对照使用；它不修改 YAeHMOP，
也不属于 xTB 描述符目录或研究协议。

## API 方向

公开计算接口使用可复用的计算器对象：

```python
from rdkit_eth import XTBCalculator, calc_xtb_descriptors

calculator = XTBCalculator(method="GFN2-xTB", electronic_temperature=300.0)
result = calculator.calculate(mol, keep_overlap=True)
features = calc_xtb_descriptors(result)
```

描述符注册表应统一提供名称、单位、所需结果能力、不变性和缺失值语义。当前
的 `descriptor_catalog()` 和 `descriptor_status()` 已经是这一注册表的基础。

## 发布方向

Python/RDKit 描述符 API 与原生 xTB runtime 应当分开：

- 核心包提供结果模型、描述符提取、报告和 RDKit 映射；
- 可选 native extra 或 conda 包提供扩展版 `libxtb` 和 `rdXTBTools`；
- 开发者仍可使用源码安装和可复现的 Linux 研究环境；
- 普通 wheel 不应携带完整编译器和 conda runtime。

## 下一阶段

下一阶段重点是计算器 facade、描述符注册表、后端能力发现、平台打包和示例。
几何优化保留为可选坐标准备工具。新增描述符必须对应明确的 RDKit 分析需求，
并通过单位、不变性、守恒和 ctypes/compiled 一致性检查后，才能进入推荐目录。

EAS 是验证这个扩展层的研究应用，应放在 research 工具中，而不是定义核心 API。
