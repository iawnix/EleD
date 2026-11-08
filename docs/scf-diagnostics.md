# Successful SCC diagnostics (0.11)

The reviewed libxtb extension adds `xtb_rdkit_getSCFAPIVersion()` and
`xtb_rdkit_getSCFDiagnostics(env, results, iterations, values)`. This is a
project ABI version 1 alongside the upstream xTB API and basis-export ABI. The
selected binary must match its build record and reviewed diagnostic patch;
official or older libraries are rejected by the public calculation APIs.

After a successful GFN1/2 single point, `result.metadata['scf']` contains:

| Field | Meaning |
|---|---|
| `iterations` | Actual SCC iteration count |
| `energy_change_hartree` | Absolute energy change in the final SCC iteration |
| `mixer_rms` | Final RMS of the engine's mixed-variable differences |
| `energy_threshold_hartree` | Energy convergence threshold used by the engine |
| `mixer_threshold` | Mixed-variable convergence threshold used by the engine |
| `criteria_satisfied` | Both measured values are below their engine thresholds |

GFN2 mixing includes shell charges, atomic dipoles and quadrupoles. Its RMS is
**not simply an atomic-charge residual**, and no elementary-charge unit is
assigned to this mixed metric. These quantities describe numerical convergence,
not physical accuracy or proof of the lowest SCF solution.

`CalcDescriptorReport` exposes the record under `quality.scf_details`, with
`scf_iterations` and `scf_residual` as compatibility fields. An observed failure
of the reported criteria sets an error quality flag. The scalar descriptor
definitions and prediction feature sets do not change.

Compiled `XTBResults` exposes `scfIterations`, `scfEnergyChangeHartree`,
`scfMixerRMS`, `scfEnergyThresholdHartree` and `scfMixerThreshold`. Raw C++ fields
use -1 for unavailable data; shared Python extractors return explicit missing
records. GFN0 reports `not_applicable_non_scf_method`.

The C getter writes one integer and four doubles only when data are available.
New results, failed/repeated calculations, Hessian and CPCM-X entry points clear
the previous snapshot; result copies own their diagnostics. This first ABI
exports successful SCC diagnostics only. A failed SCF still has an error and
requested settings, not a fabricated residual or a recoverable iteration trace.

The source builder runs a standard-library C-ABI check covering availability,
GFN0, successful GFN1/2, copied results, failed reuse and invalidation. Integration
tests compare ctypes and compiled values. Updating the source lock adds the
diagnostic patch after the existing basis patch and verifies final source hashes.
Rebuild the native engine and compiled adapter to obtain the new capability.
