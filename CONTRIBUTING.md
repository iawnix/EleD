# Contributing / 贡献指南

## English

Issues and pull requests are welcome. Before changing a public descriptor or
report field, open an issue describing the use case, units, invariance and
missing-value semantics.

For code changes:

1. Create a focused branch and keep unrelated history out of the change.
2. Add or update a meaningful test for behavior that can regress.
3. Run `python -m pytest -q` in the matched runtime.
4. Update the relevant English and Chinese documentation when the public API
   changes.
5. Record engine, parameter and protocol changes in the report provenance or
   release notes.

Descriptor additions should identify their physical definition, units, index
ordering, invariance and failure semantics. Do not silently convert an
unavailable value to zero. Pull requests that change calculation results should
include a small reproducible input and the before/after numerical tolerance.

## 中文

欢迎提交 Issue 和 Pull Request。修改公开描述符或报告字段前，请先说明使用
场景、单位、不变性和缺失值语义。

代码修改应当：

1. 使用专门分支，避免混入无关修改；
2. 为可能回归的行为添加或更新有意义的测试；
3. 在匹配的运行环境中执行 `python -m pytest -q`；
4. 公开 API 变化时同步更新英文和中文文档；
5. 在报告来源或 release notes 中记录引擎、参数和协议变化。

新增描述符必须说明物理定义、单位、索引顺序、不变性和失败语义。不可用值
不能静默转换为零。改变数值结果的 Pull Request 应提供最小可复现输入和前后
结果的数值容差。
