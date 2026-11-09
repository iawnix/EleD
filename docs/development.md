# Developing EleD

[English](development.md) | [中文](zh-CN/development.md)

First complete the [standard installation](install.md). It creates
`.eled/runtime` with RDKit, EleD's project-built `libxtb` and `rdXTBTools`.

## Edit Python code

From the repository root, activate that environment and install the source in
editable mode:

```bash
conda activate "$PWD/.eled/runtime"
python -m pip install --no-deps --no-build-isolation -e .
python -m eled doctor --check native --method GFN2-xTB --keep-overlap
python -m pytest -q
```

This command updates the EleD Python package in an already configured environment.
It reuses the installed calculation library; `pip` does not build `libxtb`.
Python source edits take effect when you restart Python or reload the module.

## Edit native code

For `rdXTBTools` C++ changes, follow the [compiled interface build guide](rdxtb-tools.md).
For xTB source or patch changes, follow the [calculation library build guide](xtb-patched-build.md).
Restart Python after rebuilding a loaded native library.
