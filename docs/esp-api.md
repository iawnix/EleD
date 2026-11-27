# In-memory xTB electrostatic potential (ESP)

The two fitting functions validate the same sampled points and potentials as `CalcESP`. Fitting also checks the indexed source geometry against the atomic numbers and the recorded atom count when available.

The extended `libxtb` exports an ESP C ABI for nonperiodic GFN0-, GFN1- and
GFN2-xTB single points. Upstream already compiles the ESP integral routines
into the shared library, but exposes no corresponding public C API. Its
Fortran plotting entry point reads/writes files and uses a fixed temporary
filename; it is not a stable interface for external callers.

EleD adds an opt-in result snapshot and evaluates arbitrary supplied points
entirely in memory, using the same primitive integrals and screening as native
`xtb --esp`. No additional SCF or command-line process is needed. ESP evaluation
is separate from the optional ESP/RESP charge fitting described below.
The EleD and rdXTBTools interfaces sample ESP during a singlepoint and return
owned, serializable numerical arrays. ESP/RESP fits use those same arrays.

## Python workflow

```python
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import (XTBCalculator, GenerateESPGrid, CalcESP,
                  FitESPCharges, FitRESPCharges)

mol = Chem.AddHs(Chem.MolFromSmiles("O"))
AllChem.EmbedMolecule(mol, randomSeed=42)
points = GenerateESPGrid(mol, scales=(1.4, 1.6, 1.8, 2.0), pointsPerAtom=100)
result = XTBCalculator(method="GFN2-xTB").Calculate(mol, espPointsAngstrom=points)
potential = CalcESP(result)  # tuple, Hartree/e
esp = FitESPCharges(result)
resp = FitRESPCharges(result, equivalenceGroups=[(1, 2)])
print(esp.charges, resp.GetCharges())  # original RDKit atom order, e
record = resp.ToDict()  # JSON-compatible fit, diagnostics, settings and provenance
```

`espPointsAngstrom` accepts a finite real `[npoints, 3]` array or sequence;
`None`/empty disables sampling. `XTBOptions` freezes a private copy. The sampled
coordinates and values are returned as `result.espPointsAngstrom` and
`result.espValuesAU` (or `result.GetESPPoints()` / `result.GetESPValues()`). Sampling uses the same singlepoint as the other descriptors.
`CalcESP` reads these values; it does not evaluate new points on an old result.
Supply the desired points when calculating. All native handles are released
before the Python result is returned; disk caches and worker processes need no
live engine objects. Grid coordinates participate in cache and checkpoint identity.
Batch options apply the same supplied grid to each molecule; generate a separate
molecular grid for each input when appropriate.

`CalcESP` raises `DescriptorUnavailableError` for an unsampled result and
`ValueError` for malformed stored data. Sampling requires a verified ESP-capable
library; no substitution of point-charge potentials is made.

The compiled equivalent is:

```python
from eled import rdXTBTools
ok, compiled = rdXTBTools.RunMol(mol, method="GFN2-xTB", espPointsAngstrom=points)
if not ok:
    raise RuntimeError(compiled.errorMessage)
xyz = compiled.GetESPPoints()  # independent NumPy copy, Angstrom [point, xyz]
v = compiled.GetESPValues()   # independent NumPy copy, Hartree/e [point]
resp = FitRESPCharges(compiled, equivalenceGroups=[(1, 2)])
```

`compiled.espSampled` reports availability. The common `CalcESP`, `FitESPCharges`
and `FitRESPCharges` functions accept either result type. Standalone C++ uses
`XTBOptions::espPointsAngstrom` and `XTBResults::espPointsAngstrom/espValuesAU`
(flat point-major vectors). Existing Mulliken/EEQ/CM5 fields retain their meanings;
fitted charges are returned in a separate `ESPChargeFit` with `chargeModel`.

## Sampling and charge-fitting definitions

`GenerateESPGrid(mol, confId=-1, scales=(1.4,1.6,1.8,2.0), pointsPerAtom=100)`
uses RDKit vdW radii and fixed lab-frame Fibonacci directions on each atom/shell.
Points buried inside another atom's sphere on that shell are removed. Explicit
hydrogen atoms and a finite 3D conformer are required. This is a documented
sampling rule, not CHELPG or Merz–Kollman. Store the returned coordinates and
chosen grid settings. A newly generated grid changes under molecular rotation;
for rigid-motion comparisons transform the same grid with the molecule.

Both fits use `A[k,i] = 1 / distance_bohr(point_k, atom_i)`. The target is the
sampled molecular potential, and total charge is read from the calculation's
provenance. The fitted charge unit is e. Interfaces are:

```python
FitESPCharges(result, *, weights=None, equivalenceGroups=(), fixedCharges=None,
              rcond=1e-12)
FitRESPCharges(result, *, restraintStrength=0.0005, restraintWidth=0.1,
               restrainHydrogens=False, weights=None, equivalenceGroups=(),
               fixedCharges=None, maxIterations=1000, tolerance=1e-8, rcond=1e-12)
```

ESP minimizes `0.5 * sum(w * (Aq - V)**2)`. Single-stage RESP additionally uses
`sum(a_i * (sqrt(q_i**2 + b**2) - b))`, where `b` is `restraintWidth` in e,
and `a_i` is `restraintStrength` (zero on H by default). Strength is the
coefficient in this numerical atomic-unit objective. Weights multiply squared
residuals, default to one, may be zero, and must be nonnegative with a positive
total. They are **not normalized**: changing their scale or sampling density
changes the relative effect of the restraint.

Total charge, user-specified equalities (`equivalenceGroups`) and fixed values
(`fixedCharges={atom_index: charge_e}`) are enforced exactly through a null-space
solve. Indices follow the input molecule. Equivalent atoms are not inferred from
the graph. Inconsistent constraints fail. ESP requires full rank in the remaining
free parameters. RESP can resolve deficient data only if its restraints remove
the remaining ambiguity; otherwise it also fails. `rcond` controls rank decisions.
RESP iterations stop when the maximum charge change is at most `tolerance` e;
nonconvergence raises `RuntimeError` and does not return a successful fit.

`ESPChargeFit` includes `charges` (`GetCharges()`), predicted potentials, residuals (fit minus
reference), weighted RMS, maximum residual, constraint error, data rank, number
of free parameters, condition number, iteration count and objective. A deficient
data matrix has `conditionNumber=None`, including when a RESP fit is determined
by restraints. Settings, source engine/parameter/calculation provenance and a
SHA256 of numerical inputs plus fitting settings accompany the result. Protocol
IDs are `eled-esp-v1` and `eled-resp-single-stage-v1`.

This RESP interface implements one hyperbolically restrained stage. It does not
claim the complete Amber two-stage protocol, prescribed QM level/grid, or a
validated force-field charge set. Multiple conformers can be analyzed separately;
joint multi-conformer fitting and automatic chemical equivalences are outside
this interface. Finite-basis density, geometry, sampling and restraints affect
the fitted charges; inspect residuals and conditioning before interpreting them.

## Politzer–Murray surface ESP descriptors

`GenerateESPSurfaceGrid` and `CalcESPSurfaceDescriptors` compute molecular and
per-atom Politzer–Murray electrostatic-surface statistics from sampled ESP
results without invoking xTB again:

```python
from eled import XTBCalculator, GenerateESPSurfaceGrid, CalcESPSurfaceDescriptors

grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=100)
result = XTBCalculator(method="GFN2-xTB").Calculate(
    mol, espPointsAngstrom=grid.pointsAngstrom
)
surface = CalcESPSurfaceDescriptors(mol, result, surfaceGrid=grid)
print(surface["v_min_au"], surface["v_max_au"], surface["pi_esp_au"], surface["nu"])
print(surface["atoms"][0]["v_mean_au"], surface["atoms"][0]["sigma2_au2"])
```

- `GenerateESPSurfaceGrid(mol, confId=-1, *, scales=(1.4,), pointsPerAtom=100)`
  builds exposed Fibonacci spheres around RDKit periodic-table van der Waals
  radii scaled by `scales`, returning an immutable `ESPSurfaceGrid` with
  `pointsAngstrom` (Å), `atomIndices` (original RDKit atom order), and
  `areaWeightsAngstrom2` (Å²). This is an explicit scaled vdW shell surface,
  **not** a 0.001 a.u. quantum-mechanical electron-density isosurface.
- `CalcESPSurfaceDescriptors(mol, result, *, atomIndices=None, weights=None, surfaceGrid=None, includeKcalMol=True)`
  reports molecular `v_min_au`, `v_max_au`, `v_mean_au`, `v_pos_mean_au`,
  `v_neg_mean_au`, `v_pos_frac`, `v_neg_frac`, `pi_esp_au`, `sigma2_pos_au2`,
  `sigma2_neg_au2`, `sigma2_tot_au2`, `nu`, `nu_sigma2_au2` and per-atom
  `v_min_au`, `v_max_au`, `v_mean_au`, `v_pos_mean_au`, `v_neg_mean_au`,
  `sigma2_au2` (plus kcal/mol and (kcal/mol)² counterparts when
  `includeKcalMol=True`). Single-sign surfaces or buried atoms return `None`
  with explicit `missing` reasons rather than zero.

## Exact electrostatic potential at nuclei (EPN)

The extended properties C ABI v1 also exports `xtb_rdkit_getNuclearESP`,
evaluating the exact valence electronic plus other-nuclei electrostatic
potential $V_N$ at each nuclear position while omitting the singular
self-nuclear term $Z_A / r_{iA}$:

```python
from eled import XTBCalculator, CalcNuclearESP

result = XTBCalculator(method="GFN2-xTB").Calculate(mol)
epn = CalcNuclearESP(result, mol=mol)
print(epn["epn_min_au"], epn["atoms"][0]["epn_au"], epn["atoms"][0]["epn_kcalmol"])
```

`CalcNuclearESP(result, *, mol=None)` accepts native `XTBResults` and compiled
`rdXTBTools.XTBResults` (`result.GetNuclearESP()`), preserves original RDKit
atom order, and reports both Hartree/e (`epn_au`) and kcal/(mol·e)
(`epn_kcalmol`). Values are rigid-motion and atom-permutation invariant and
match $\lim_{\mathbf{r}\to\mathbf{R}_A}\bigl(V_{\mathrm{ESP}}(\mathbf{r}) - Z_A / |\mathbf{r}-\mathbf{R}_A|\bigr)$.

## C interface

```c
int xtb_rdkit_getESPAPIVersion(void);  /* returns 1 */
void xtb_rdkit_setKeepESP(xtb_TEnvironment env, xtb_TResults res, bool keep);
int xtb_rdkit_getESP(xtb_TEnvironment env, xtb_TResults res,
                    int npoints, const double *points,
                    int capacity, double *output);
```

With an initialized environment, molecule, loaded calculator and result:

```c
xtb_rdkit_setKeepESP(env, res, true);
xtb_singlepoint(env, mol, calc, res);
/* Check xtb_checkEnvironment(env) before reading results. */
double points[2][3] = {{4.1, 2.3, 3.7}, {-3.2, 4.5, 2.1}};
double potential[2];
int ok = xtb_rdkit_getESP(env, res, 2, &points[0][0], 2, potential);
/* Read potential only if ok == 1; otherwise inspect xtb_getError. */
```

Points are contiguous `double[npoints][3]` in **Bohr**, in the calculation's
coordinate frame. Output is `double[npoints]` in **atomic units of electric
potential (Hartree/e)**. Its sign is the potential energy sign for a positive
unit test charge. The caller allocates both buffers; `capacity` must equal
`npoints > 0`. The output order exactly follows the supplied point order.

The potential combines the molecular effective nuclear cores and electronic
density in xTB's valence basis. It excludes solvent reaction fields and external
point-charge contributions, even if those polarized the underlying wavefunction.
It is not an all-electron potential near atomic cores or a periodic potential.

The function returns 1 on success. Missing snapshots, unsupported calculations,
null pointers, invalid counts/capacities, nonfinite coordinates, points within
1e-10 Bohr of a nucleus, and nonfinite computed values return 0 with an
`xtb_getError` diagnostic. Output buffers remain untouched on error; unavailable
values are never encoded as zero. A null environment returns 0 without a diagnostic.

## Ownership and numerical scope

Retention defaults to false. Enable it **before** `xtb_singlepoint`; enabling it
afterwards does not reconstruct a prior result. Snapshots own geometry, basis,
occupations and Cartesian orbital coefficients, independently of molecule and
calculator lifetimes. `xtb_copyResults` deep-copies the snapshot. Any retention
setter, new singlepoint attempt (including failed setup or SCF), Hessian or
CPCM-X call invalidates the previous snapshot. A subsequent supported successful
singlepoint restores it when retention remains enabled. Storage is O(NAO²).
GFN-FF and periodic calculations do not produce an ESP snapshot.

Screening matches CLI defaults: occupation cutoff 1e-6, packed-density cutoff
1e-4, primitive exponent cutoff 10, primitive density/contraction cutoff 1e-5,
and weighted primitive cutoff 5e-7. Summation order differs from the CLI.
Cartesian screening can limit strict rotational invariance; this is a numerical
approximation, not a changed charge partition. The API does not establish general
libxtb thread safety; callers should serialize engine operations as EleD does.

Validation covers points at molecular-surface distances. At extreme distances
(for example 10,000 Bohr), the upstream integral routine can print convergence
warnings while returning finite values. This extension retains that behavior;
finite output alone does not certify accuracy in that regime.

## Build and verification

Reuse the fixed runtime and toolchain from the standard installation:

```bash
conda run --no-capture-output -p "$PWD/.eled/toolchain" \
  python tools/build_patched_xtb.py "$PWD/.eled/runtime" \
  --dependency-prefix "$PWD/.eled/toolchain" --jobs 6
.eled/runtime/bin/python tools/diagnostics/verify_esp.py .eled/runtime
```

The builder verifies the ESP patch and resulting source hashes, runs upstream
CTest plus existing basis/SCC checks, and compares the installed ESP ABI with
CLI values before writing `capabilities.esp_api = 1` and `esp_checks` in the
build manifest. Restart processes that already loaded the old library.

The verifier covers water, hydrogen sulfide (Cartesian d conversion), ammonium,
and the hydroxyl radical under all three methods; result copying, disabled
retention, invalid inputs, failed setup/SCF, absence of ESP files, and water
rigid-motion/atom/point permutation checks. It also runs as an integration test
in `tests/test_esp.py`. CPCM-X is disabled in the standard build.
