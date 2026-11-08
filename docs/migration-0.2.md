# Migrating from 0.1 to 0.2

Version 0.2 removes the external xTB CLI calculation backend. The default
`run_xtb(mol)` and `calc_orbital_descriptors(mol)` now use libxtb directly.
The optional `rdXTBTools.RunMol()` C++ interface and its units are unchanged.

| Previous API | 0.2 behavior / replacement |
| --- | --- |
| `backend="cli"` or `backend="auto"` | Rejected. Use `backend="native"` (default) or explicit `backend="tblite"`. |
| `xtb_path` / `XTB_PATH` | Executable discovery removed. Use `xtb_library` / `XTB_LIBRARY` to select a shared library. |
| `keep_files`, `timeout` | Removed; calculations run in process and do not create work directories. No wall-clock timeout is enforced. |
| `optimize`, `solvent` | Removed. The current API supports gas-phase single points only. |
| `add_optimized_conformer` | Removed with the geometry optimization API. |
| `result.command`, `result.workdir`, `result.optimized_coordinates_angstrom` | Removed. Numerical results and provenance remain on the result object. |
| `XTBNotFoundError` | Removed with executable lookup. Catch `XTBNativeNotFoundError` for library-loading failures. |
| `run_xtb_batch(..., max_workers>1)` | Rejected; batch execution is sequential. Cache and per-molecule error handling remain available. |
| `doctor --check cli`, `doctor --executable` | Removed. Use `doctor --check native --library /path/libxtb.so`. |

Removed keyword arguments raise `TypeError`; removed backend values raise
`ValueError`. No unsupported option is silently ignored. `XTBPATH` still
selects parameter files and is distinct from the removed `XTB_PATH` variable.

A missing native library raises an error without switching engines. Native
GFN0 requires the verified patched xTB build; an official/unverified library
is rejected before GFN0 initialization. Install using [these instructions](install.md).

```python
from eled import run_xtb
result = run_xtb(mol, xtb_library="/absolute/prefix/lib/libxtb.so")
print(result.total_energy_hartree)
print(result.metadata["engine"])
```

The conda xTB dependency and the upstream source builder still contain an
executable as part of that distribution. The adapter does not launch it for
calculations. Upstream diagnostic scripts and archived 0.1 evidence may still
refer to executable runs; those records describe historical experiments.
The `eled doctor` command remains an installation diagnostic utility.

Existing `dist/0.1.0` files are historical artifacts and still contain the old
backend. Use the 0.2 source or newly built 0.2 wheels. Compiled wheels continue
to require the matching local SDK and xTB prefixes.
