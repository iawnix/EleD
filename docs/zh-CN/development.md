# 开发 EleD

[English](../development.md) | [中文](development.md)

先完成[标准安装](install.md)。安装器会在 `.eled/runtime` 中配置 RDKit、
EleD 项目构建的 `libxtb` 和 `rdXTBTools`。

## 修改 Python 代码

在仓库根目录激活该环境，以可编辑模式安装源码：

```bash
conda activate "$PWD/.eled/runtime"
python -m pip install --no-deps --no-build-isolation -e .
python -m eled doctor --check native --method GFN2-xTB --keep-overlap
python -m pytest -q
```

这一步用于更新已配置环境中的 EleD Python 包，复用完整安装提供的计算库；
`pip` 本身不构建 `libxtb`。修改 Python 源码后，重启 Python 或重新加载模块
即可使用新代码。

## 修改原生代码

修改 `rdXTBTools` C++ 代码时，按 [C++ 接口构建指南（英文）](../rdxtb-tools.md)
重新编译。修改 xTB 源码或补丁时，按[计算库构建指南（英文）](../xtb-patched-build.md)
重新构建。替换已经加载的原生库后，需要重启 Python。
