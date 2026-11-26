# 离线 Linux 运行时安装包

[English](../runtime-bundle.md) | [中文](runtime-bundle.md)

离线运行时包包含匹配的 RDKit/Python 运行时、已编译的 `rdXTBTools` 扩展、项目构建的
扩展版 `libxtb`、xTB 参数文件以及独立的原生依赖前缀（`toolchain`）。安装过程不需要
conda、编译器、网络连接或单独安装 xTB 命令行程序。

该产物面向 Linux x86-64 环境分发，不是通用的 manylinux wheel。目标机器的
CPU/glibc 兼容性仍需在对应平台上单独验证；Windows/macOS 与公共包索引发布仍不在
当前支持范围内。

## 离线安装

拿到包含 `install.py` 和 `bundle.json` 的可信离线包目录后，将其安装到一个新的目标目录：

```bash
python3 /path/to/bundle/install.py --directory /absolute/new/eled
/absolute/new/eled/runtime/bin/python -m eled doctor --check compiled --keep-overlap
```

安装脚本需要 Python 3.10+（带 `tarfile` `data` 过滤器的版本，或 Python 3.12+）。
安装器会执行以下步骤：

1. 校验 `bundle.json` 中记录的全部载荷 SHA256；
2. 解压 `runtime` 与 `toolchain` 两个前缀并运行 `conda-unpack`；
3. 使用随包携带的 `patchelf` 修复跨前缀的 ELF RPATH，并将重定位前后的哈希写入
   `share/xtb/eled/build.json`；
4. 在源码目录之外且清除环境变量覆盖的隔离 Python 进程中，验证 `native` 与
   `compiled` 两种接口下的 GFN0/1/2、overlap/布居与 GFN1/2 反应性计算；
5. 检查 `ldd` 输出，拒绝存在缺失符号或残留原始构建路径的动态库依赖。

安装完成后的验证记录写入目标目录的 `installation.json`。使用期间需将 `runtime` 与
`toolchain` 一起保留在安装路径；如需安装到其他路径，请从原始离线包重新运行安装器。
失败的安装目录不会被自动删除，便于排查问题。

Python 的 `compiled` 与 `native` 路径均自动发现当前运行时前缀中的 `libxtb.so`；
显式传入的 `xtbLibrary` 或 `XTB_LIBRARY` 仍保持最高优先级。独立的 C++ 消费程序在
重定位后应通过 `XTBOptions.libraryPath` 指定库路径。

## 构建离线包

先完成并验证标准的 wheel 安装（打包脚本会在打包前自动检查 `runtime` 中安装的
`eled` 位于该前缀内部而非可编辑源码目录，且实际加载的是前缀内的 `libxtb.so`）。
然后在独立的打包虚拟环境中安装构建依赖并打包：

```bash
python -m pip install --require-hashes -r tools/distribution-requirements.txt
python tools/build_runtime_bundle.py \
  --runtime /path/to/installation/runtime \
  --toolchain /path/to/installation/toolchain \
  --output /new/bundle-directory
```

打包产物保留 conda 元数据与许可证、对应版本的补丁后 xTB 源码及许可证，以及 EleD
的 BSD-3-Clause 许可证。因完整保留了锁定 SDK 与依赖环境，包体积较大，并非最小化
镜像。构建与安装工具均不会执行任何网络上传。校验和用于检测文件意外损坏并关联来源，
不替代发布者签名，请仅使用来自可信来源的离线包。
