# 构建原生 xTB 依赖（0.11）

标准安装方式是 `python tools/install_native.py`，见 [install.md](install.md)。这个更底层的
构建脚本还支持为开发单独指定引擎前缀。它使用上游的 GFN0 状态修复，并应用本项目的基组
导出、SCC 诊断、ESP、扩展性质和实空间（ABI v2）补丁。GFN0/1/2 和布居分析都使用这一个
扩展库，不涉及命令行程序或 tblite 后端。

## 复现构建

在 Linux x86-64 的源码目录中执行：

```bash
conda create --yes --prefix "$PWD/.xtb-build-env" --file locks/linux-64/xtb-toolchain.txt
conda run --no-capture-output --prefix "$PWD/.xtb-build-env" \
  python tools/build_patched_xtb.py "$PWD/.xtb-native-0.5" \
  --dependency-prefix "$PWD/.xtb-build-env" --jobs 4
source .xtb-native-0.5/activate.sh
```

编译器环境与 RDKit 运行时/SDK 相互独立。两者以及 `.xtb-native-0.5` 都需要保留：本地
二进制通过 RPATH 查找编译器和数值库。这些二进制不是可移植的发行包。构建时关闭
`WITH_TBLITE` 和 `WITH_CPCMX`；xTB 自带的 ALPB/GBSA 隐式溶剂化保留在引擎中，并通过
`XTBOptions(solvent=...)` 提供。上游可执行程序只为溯源和诊断而构建，适配器不会调用它。

[源码锁定文件](../../tools/xtb-source.json) 固定了提交
`33fe381adde775c4992610c0f90d8914b7572bb5`、上游修复后的 D4 源码、
基组/诊断/ESP/性质/实空间补丁，以及所有被修改源文件的 SHA256。这个上游 main 分支提交
自报版本为 6.7.1，但与发布标签不同，见[上游 GFN0 修复记录](../upstream/xtb-gfn0-upstream-fix.md)（英文）。
构建脚本接受干净的固定版本检出，或与补丁完全一致的已打补丁检出；遇到无关或不一致的修改
会拒绝执行，并且从不重置它们。重复运行会复用增量构建。可选的第二个位置参数用于指定已有
的检出目录。修改源码锁定文件时请使用新的前缀，以保留之前验证过的引擎。

构建会在写入能力元数据之前运行所有上游测试、12 项已安装库的 D4 方法切换检查、原生基组/
生命周期 C 诊断、SCC 诊断生命周期检查，以及 ESP C ABI/命令行对比。没有跳过测试的选项。
安装内容包括：

```text
.xtb-native-0.5/
  source/                     # 固定版本的源码及可见的本地 diff
  build/                      # 上游构建和 C 诊断
  lib/libxtb.so
  include/xtb.h
  bin/xtb
  share/xtb/                  # 匹配的 GFN 参数文件
  share/xtb/eled/             # build.json、source.diff、补丁、许可证、测试日志
  activate.sh
```

激活脚本设置 `XTB_LIBRARY` 和 `XTBPATH`，不修改 `PATH`。更换已加载的库之后请重启
Python。按 [install.md](install.md) 的说明，使用匹配的 RDKit SDK 构建适配器。

## 基组扩展 ABI v1

上游的 `xtb_getAPIVersion()` 不变。可选导出的符号有：

```c
int xtb_rdkit_getBasisAPIVersion(void);
void xtb_rdkit_setKeepOverlap(xtb_TEnvironment, xtb_TResults, bool);
void xtb_rdkit_getBasisSize(xtb_TEnvironment, xtb_TResults, int *sizes);
void xtb_rdkit_getOverlap(xtb_TEnvironment, xtb_TResults, int capacity, double *out);
void xtb_rdkit_getBasisMap(xtb_TEnvironment, xtb_TResults, int nao, int nshell,
                         int *ao_atom, int *ao_shell, int *shell_atom, int *shell_l);
```

需要在 `xtb_singlepoint` 之前对结果对象开启保留。实际的重叠矩阵和基组映射在 GFN0 或
GFN1/2 计算中采集，属于该结果的波函数。S 与 C 使用相同的原生球谐 AO 顺序。C 缓冲区按
列主序存储；两个 Python 接口都返回 `[AO, AO]` 的 S 和 `[AO, MO]` 的 C。尺寸输出恰好是
三个整数：原子数、壳层数、AO 数。映射索引从零开始；壳层角动量 s=0、p=1、d=2。AO 映射
缓冲区需提供 `nao` 个元素，壳层映射需提供 `nshell` 个。重叠矩阵容量必须恰好为
`nao * nao`；无效或空缓冲区会产生环境错误且不写入输出。与其他 C 缓冲区 API 一样，调用方
负责按声明的容量分配内存。

保留功能默认关闭，会额外占用 O(NAO²) 存储。复制结果时会深拷贝其保留的数组。更改保留
设置或开始新的计算会清除之前的快照。失败的单点计算不会暴露旧的基组数据。Hessian 和
CPCM-X 入口会使波函数快照失效并关闭它；之后如果结果的保留标志仍开启，下一次成功的单点
计算会恢复快照。Hessian 失效由 C 校验程序覆盖；CPCM-X 在本构建中关闭，未做数值测试。

只有当扩展符号/ABI、审查过的补丁哈希、清单中的能力声明以及所选二进制的哈希全部一致时，
适配器才会启用这项能力。`capabilities.orbital_basis_api` 为 1。GFN0 还单独要求经过验证
的上游修复提交/源码身份，以及 `capabilities.native_gfn0: true`（0.17 起不再接受使用
0.4 时代本地 D4 补丁构建的安装）。引擎版本字符串不能证明任何一项能力。公开适配器对每次
GFN0/1/2 计算都要求完整的扩展配置。

## ESP 扩展 ABI v1

该库可以保留单点计算的 ESP 快照，并通过 `xtb_rdkit_setKeepESP` 和 `xtb_rdkit_getESP`
在内存中计算指定空间点的静电势。单位、所有权、错误处理、示例和验证见
[ESP C API](esp-api.md)。该扩展记录为 `capabilities.esp_api = 1`。两个 Python 接口都
支持 ESP 采样；共享的 ESP/RESP 拟合函数可用于任一种结果类型。

## 验证与适用范围

原生校验程序检查 GFN0/1/2 的 CᵀSC、壳层映射、电子守恒、GFN1/2 Mulliken 电荷重建、结果
复制、关闭保留、无效缓冲区、SCF 失败、Hessian 失效与恢复。GFN0 的引擎电荷是 EEQ，因此
不要求它等于 MO 的 Mulliken 电荷。适配器测试覆盖 66 个 GFN1/2 冻结构象算例、自由基/
阳离子计算、编译接口与 ctypes 的一致性，以及简并子空间不变性。实测结果见
[当前验证记录](../validation.md)（英文）。

这并不说明 libxtb 在一般意义上是线程安全的。适配器的调用（包括两个 Python 接口）是串行
化的；独立的其他绑定不在这把锁的保护范围内。Python 适配器还基于该单点 API 提供梯度
优化器；激发态计算不在适配器范围内。

标准安装程序把这个引擎构建与锁定的 RDKit 依赖以及适配器 wheel 安装组合在一起。conda
配方只是打包脚手架，不提供可用的纯上游运行时。普通的 `pip install` 只安装适配器，不包含
Fortran 工具链；请使用标准安装程序或本构建脚本获得所需的扩展库。可移植的发行版需要单独
的打包和干净环境检查。xTB 采用 LGPL-3.0-or-later 许可；两份许可证文本、确切的补丁、源码
diff 和提交信息都随安装一起保留。各依赖有各自的许可证。
