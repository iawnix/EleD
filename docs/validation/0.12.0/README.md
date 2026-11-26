# 0.12.0 release validation

Validated locally on 2026-11-26 on Linux x86-64 / Python 3.11 / RDKit 2026.09.1 with `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.

| Check | Result |
| --- | --- |
| Full extended-engine test suite (`.eled/runtime`, `ELED_REQUIRE_RDXTB=1 ELED_EXPECT_BASIS=1 ELED_EXPECT_NATIVE_GFN0=1`) | 618 passed, 12 skipped, 21 benchmark deselected |
| Pure-Python test suite (`.runtime-env`, `-m 'not integration and not benchmark'`) | 410 passed, 11 skipped |
| Standalone C++ core CTest (`rdXTBTools/core`, including native `calcCM5Charges`) | 1/1 passed |
| All 10 runnable scripts in `examples/` (including `reactivity_batch.py` checkpoint resume) | 10/10 passed |
| In-process vs CLI xTB geometry optimization (`tools/benchmark_xtb_optimization.py --check`, water & ethanol) | Passed |
| Native vs CLI `--vfukui` comparison (`tools/compare_fukui_cli.py`, `tests/data/reactivity_cases.json`) | Passed |
| Clean sdist + pure wheel (`eled-0.12.0-py3-none-any.whl`) outside-checkout venv verification | Passed |
| Clean sdist + compiled wheel (`eled-0.12.0-cp311-cp311-linux_x86_64.whl`) outside-checkout venv verification (`GFN2-xTB` native/compiled + overlap + reactivity; `GFN0-xTB` native/compiled + overlap) | Passed |
| Inherited `--system-site-packages` venv without `XTB_LIBRARY`/`XTBPATH` overrides (`native` and `compiled`) | Passed (`-5.0703716808 Hartree` water parity) |
| Offline Linux runtime bundle build (`tools/build_runtime_bundle.py`) and relocated offline install (`tools/install_bundle.py`, GFN0/1/2 × native/compiled + `ldd` + isolated `parallel_descriptors.py` & 35 worker/SCC/RDKit/visualization tests) | Passed (1.1 GB bundle, 35 passed in relocated runtime) |

## Scope of 0.12.0 changes verified

- Full RDKit `PascalCase` / `camelCase` public Python API alignment (`XTBResults`, `XTBOptions`, `XTBCalculator.Calculate`, `XTBCache`, descriptor/report/batch/conformer/geometry/visualization entry points).
- Single-step `CalcMolXTBDescriptors`, `RegisterRDKitProperties` (`rdMolDescriptors.Properties` integration via `Descriptors.PropertyFunctor` with per-conformer single-point deduplication), and `AnnotateMol(..., privateProps=True)`.
- Native C++ `RDKit::XTBTools::calcCM5Charges` + `XTBResults::cm5Charges` and C++/NumPy option and ESP coordinate validation in `Wrap.cpp`.
- Upgraded RDKit-native `DrawDescriptorMap` (fragment-free reports, direct `XTBResults` input, bond descriptors, continuous diverging/sequential colormaps, implicit-H-preserving `RemoveHs` mapping, and `provenance` validation).
- Protocol and storage version bumps: `DESCRIPTOR_VERSION = "1.4"`, `REACTIVITY_PROTOCOL_VERSION = "2.1"`, `STORAGE_VERSION = 2`, `XTBCache` schema `6`.

Remote GitHub Actions status must be checked separately on GitHub after pushing the release commit and tag.
