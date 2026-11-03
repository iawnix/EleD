**Describe the bug**

Calling `xtb_loadGFN0xTB` changes the energies and orbital eigenvalues of subsequently created GFN2 calculators in the same process. Loading GFN0 parameters alone is sufficient; no GFN0 single-point calculation is needed.

This reproduces with sequential, single-threaded calls and fresh environment, molecule, calculator, and results handles for every calculation, with all handles deleted between stages. The example below uses only the public C API.

Confirmed on Linux x86-64 with:

- conda-forge xTB 6.7.1 (`gfortran_ha02e54d_5`);
- an independent source build of tag `v6.7.1`;
- an independent source build of `main` at [`a9040512dd7c15793472aa20f9647d36a266f7cd`](https://github.com/grimme-lab/xtb/commit/a9040512dd7c15793472aa20f9647d36a266f7cd).

**To Reproduce**

1. Calculate GFN2 for the fixed ethanol geometry below, then delete all handles.
2. Create fresh handles, load GFN0 parameters, then delete all handles.
3. Create fresh handles and repeat GFN2 on exactly the same geometry.

Save the C program below as `repro.c`. Set `xtb_prefix` to an installation containing `include/xtb.h`, `lib/libxtb.so`, and the parameter files. For a source build, use its header/build-library paths and the source root as `xtb_param_dir`.

```bash
xtb_prefix=/path/to/xtb-installation
xtb_param_dir="$xtb_prefix/share/xtb"

cc -std=c11 -Wall -Wextra -Werror repro.c \
  -I"$xtb_prefix/include" -L"$xtb_prefix/lib" \
  -Wl,-rpath,"$xtb_prefix/lib" -lxtb -lm -o repro

env -u LD_LIBRARY_PATH -u LD_PRELOAD -u XTBPATH -u XTBHOME \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./repro "$xtb_param_dir" 0 --expect-clean
```

On the unpatched 6.7.1 builds, this prints:

```json
{"stage":"gfn2_before","energy_Eh":-11.38860400183772,"delta_energy_Eh":0,"max_orbital_delta_eV":0}
{"stage":"gfn2_after_parameter_load","energy_Eh":-11.38871517219716,"delta_energy_Eh":-0.0001111703594478541,"max_orbital_delta_eV":0.0007635684228234822}
```

The exit status is **3**, indicating a change above `1e-10 Hartree` in total energy or `1e-8 eV` in any orbital energy. API errors instead exit with status 1.

Controls, each in a fresh process with the same thread/environment settings:

- Replace `0` with `1` or `2`: both GFN1/GFN2 intermediate loads leave subsequent GFN2 unchanged, exit 0.
- Add `--run-intermediate`: performing the intermediate GFN0 single-point calculation produces the same subsequent GFN2 differences as merely loading parameters.
- Add `--water`: the independent fixed water geometry also shows the problem.

<details>
<summary>Complete standalone C reproducer (fixed coordinates in Bohr; public API only)</summary>

```c
/* Standalone xTB 6.7.1 C API reproducer; no RDKit or Python adapter required.
 * Usage: repro PARAMETER_DIRECTORY [0|1|2]
 *              [--run-intermediate] [--water] [--expect-clean]
 * Every calculation owns fresh handles and deletes them before returning.
 */
#include <xtb.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const int ethanol_numbers[9] = {6, 6, 8, 1, 1, 1, 1, 1, 1};
/* Ethanol, RDKit EmbedMolecule seed 1, coordinates in Bohr. */
static const double ethanol_xyz[27] = {
    1.6945481234523154, -0.14513687028427885, -0.028442010502842605,
    -1.0893542800733351, -0.66384098408490566, 0.11328282457504339,
    -2.4873192390810432, 1.5334078933434283, 0.38677566079066467,
    2.4044388440269864, 0.2745547719270563, 1.9143719419814178,
    1.9799246375778272, 1.5544839955851846, -1.1819828050167818,
    2.7358373903821924, -1.8244437346627014, -0.65185897389749514,
    -1.510886620841805, -1.8074743990489552, 1.8412705812374754,
    -1.6242186256041504, -1.7819806495787187, -1.5540216971434768,
    -2.1029702298389861, 2.8604299768038937, -0.8393955220240058
};
static const int water_numbers[3] = {8, 1, 1};
static const double water_xyz[9] = {0., 0., 0., 1.43, 1.11, 0., -1.43, 1.11, 0.};
static int atom_count = 9;
static const int *numbers = ethanol_numbers;
static const double *xyz = ethanol_xyz;

struct result { double energy; int nao; double *orbitals; };

static void check(xtb_TEnvironment env, const char *operation) {
    if (xtb_checkEnvironment(env)) {
        char message[4096] = {0};
        const int size = sizeof(message);
        xtb_getError(env, message, &size);
        fprintf(stderr, "%s: %s\n", operation, message);
        exit(1);
    }
}

static struct result calculate(const char *directory, int method, int load_only) {
    struct result value = {0, 0, NULL};
    const int n = atom_count, uhf = 0;
    const double charge = 0.;
    char path[4096];
    int length = snprintf(path, sizeof(path), "%s/param_gfn%d-xtb.txt", directory, method);
    if (length < 0 || (size_t)length >= sizeof(path)) exit(1);
    xtb_TEnvironment env = xtb_newEnvironment();
    if (!env) exit(1);
    xtb_setVerbosity(env, 0);
    xtb_TMolecule mol = xtb_newMolecule(env, &n, numbers, xyz, &charge, &uhf, NULL, NULL);
    check(env, "molecule");
    xtb_TCalculator calc = xtb_newCalculator();
    if (!mol || !calc) exit(1);
    switch (method) {
        case 0: xtb_loadGFN0xTB(env, mol, calc, path); break;
        case 1: xtb_loadGFN1xTB(env, mol, calc, path); break;
        case 2: xtb_loadGFN2xTB(env, mol, calc, path); break;
        default: exit(2);
    }
    check(env, "load");
    if (!load_only) {
        xtb_TResults res = xtb_newResults();
        if (!res) exit(1);
        xtb_singlepoint(env, mol, calc, res);
        check(env, "singlepoint");
        xtb_getEnergy(env, res, &value.energy);
        xtb_getNao(env, res, &value.nao);
        check(env, "results");
        if (value.nao < 1) exit(1);
        value.orbitals = calloc((size_t)value.nao, sizeof(double));
        if (!value.orbitals) exit(1);
        xtb_getOrbitalEigenvalues(env, res, value.orbitals);
        check(env, "orbitals");
        xtb_delResults(&res);
    }
    xtb_delCalculator(&calc);
    xtb_delMolecule(&mol);
    xtb_delEnvironment(&env);
    return value;
}

static double maxdiff(const double *a, const double *b, int n) {
    double difference = 0.;
    for (int i = 0; i < n; ++i) difference = fmax(difference, fabs(a[i] - b[i]));
    return difference;
}

static void report(const char *stage, struct result baseline, struct result value) {
    if (baseline.nao != value.nao) exit(1);
    printf("{\"stage\":\"%s\",\"energy_Eh\":%.16g,\"delta_energy_Eh\":%.16g,"
           "\"max_orbital_delta_eV\":%.16g}\n", stage, value.energy,
           value.energy - baseline.energy,
           maxdiff(baseline.orbitals, value.orbitals, baseline.nao) * 27.211386245988);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: %s PARAMETER_DIRECTORY [0|1|2] "
                "[--run-intermediate] [--water] [--expect-clean]\n", argv[0]);
        return 2;
    }
    int method = 0, run_intermediate = 0, expect_clean = 0;
    int first_flag = 2;
    if (argc >= 3 && argv[2][0] != '-') {
        if (strlen(argv[2]) != 1 || argv[2][0] < '0' || argv[2][0] > '2') return 2;
        method = argv[2][0] - '0';
        first_flag = 3;
    }
    for (int i = first_flag; i < argc; ++i) {
        if (strcmp(argv[i], "--run-intermediate") == 0) run_intermediate = 1;
        else if (strcmp(argv[i], "--expect-clean") == 0) expect_clean = 1;
        else if (strcmp(argv[i], "--water") == 0) {
            atom_count = 3;
            numbers = water_numbers;
            xyz = water_xyz;
        } else return 2;
    }
    struct result before = calculate(argv[1], 2, 0);
    report("gfn2_before", before, before);
    /* Loading GFN0 parameters alone is sufficient; no GFN0 SCF/singlepoint. */
    struct result intermediate = calculate(argv[1], method, !run_intermediate);
    if (run_intermediate) {
        printf("{\"stage\":\"intermediate_calculation\",\"method\":%d,\"energy_Eh\":%.16g}\n",
               method, intermediate.energy);
        free(intermediate.orbitals);
    }
    struct result after = calculate(argv[1], 2, 0);
    report("gfn2_after_parameter_load", before, after);
    int contaminated = fabs(after.energy - before.energy) > 1e-10 ||
        maxdiff(after.orbitals, before.orbitals, before.nao) * 27.211386245988 > 1e-8;
    free(before.orbitals);
    free(after.orbitals);
    return expect_clean && contaminated ? 3 : 0;
}
```

</details>

**Expected behaviour**

For an identical molecule, geometry, charge, spin, and GFN2 settings, a newly created GFN2 calculator should give the same result regardless of whether a previously deleted GFN0 calculator was initialized earlier in the process.

**Additional context**

The differences below compare GFN2 before/after GFN0 within each build. Completing a GFN0 single-point calculation gives the same differences as loading its parameters alone.

| Build | Molecule | Delta total energy (Hartree) | Maximum absolute orbital-energy change (eV) |
|---|---|---:|---:|
| Source-built v6.7.1 | Ethanol | -1.111703594478541e-4 | 7.635684228234822e-4 |
| Source-built v6.7.1 | Water | +4.229524108723126e-7 | 1.084464535799272e-5 |
| Source-built main `a9040512` | Ethanol | -1.111703626488492e-4 | 7.635450311748853e-4 |
| Source-built main `a9040512` | Water | +4.229524188659184e-7 | 1.084465219467583e-5 |
| Either source build with the candidate patch below | Both | 0 | 0 |

For each source build, I checked two molecules × three intermediate methods × two actions (parameter load only or single point), giving 12 sequences before and 12 after the patch. All four GFN0 sequences show the issue before patching; all eight GFN1/GFN2 controls are unchanged. After patching, all 12 sequences have exactly zero computed differences in this environment. Fresh GFN2 baseline energies and the recorded intermediate GFN0/GFN1/GFN2 single-point energies also remain unchanged by the patch. The conda-forge ethanol load-only test was repeated three times in fresh processes with identical results.

The source builds used GNU C/Fortran 16.2.0, CMake 4.4.5, mctc-lib 0.7.0, test-drive 0.6.1, and OpenBLAS 0.3.34. They were Release builds with shared libraries enabled and optional tblite/CPCM-X disabled. The v6.7.1 C sources were compiled with `-std=c17`.

**Root cause and candidate repair**

In [`newD4Model`](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/disp/dftd4.F90#L106-L181), the Goedecker mode used by GFN0 assigns:

```fortran
refq = clsq
refh = clsh
secq = clss
```

`refq` and `refh` are [module-level arrays](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/disp/dftd4_parameters.f90#L53-L71), initialized from `param_ref.fh`. The next GFN2 initialization resets `secq`, but does not restore `refq`/`refh`. Its default branch copies those modified tables into `dispm%q` and `tmp_hq`. The GFN0/GFN2 mode choices are in [`read_gfn_param.f90`](https://github.com/grimme-lab/xtb/blob/a9040512dd7c15793472aa20f9647d36a266f7cd/src/read_gfn_param.f90#L250-L306).

The candidate patch below keeps the original GFN2 reference tables unchanged and selects each method's tables into model/local storage, including local secondary reference charges. It applies to both tested source versions and eliminates the observed contamination without accessing internal symbols from the reproducer.

<details>
<summary>Candidate patch to src/disp/dftd4.F90</summary>

```diff
--- a/src/disp/dftd4.F90
+++ b/src/disp/dftd4.F90
@@ -113,52 +113,32 @@
    integer  :: i,ia,is,icn,j,ii,jj
    integer  :: cncount(0:18)
    real(wp) :: sec_al(23),iz,c6,alpha(23)
-   real(wp) :: tmp_hq(7,118)
+   real(wp) :: tmp_hq(7,118), secondary_q(17)
 
    intrinsic :: nint
 
    call init(dispm)
 
-   secq = 0.0_wp
-   select case(mode)
-   case(p_refq_hirshfeld,p_refq_periodic)
-!     print'(1x,''* using PBE0/def2-TZVP Hirshfeld charges'')'
-      refq = dftq
-      refh = dfth
-      secq = dfts
-!  case(2)
-!     refq = pbcq
-!     refh = pbch
-!     secq = pbcs
-   case(p_refq_gasteiger)
-!     print'(1x,''* using classical Gasteiger charges'')'
-      refq = gffq
-      refh = gffh
-      secq = gffs
-   case(p_refq_goedecker)
-      refq = clsq
-      refh = clsh
-      secq = clss
-   case(p_refq_gfn2xtb_gbsa_h2o)
-!     print'(1x,''* using GFN2-xTB//GBSA(H2O) charges'')'
-      refq = solq
-      refh = solh
-      secq = sols
-   end select
-
+   ! Keep method-specific reference charges local to this model. The module
+   ! refq/refh tables contain the original GFN2 references and must not change.
+   secondary_q = 0.0_wp
    select case(mode)
    case(p_refq_hirshfeld,p_refq_periodic)
       dispm%q = dftq
       tmp_hq = dfth
+      secondary_q = dfts
    case(p_refq_gasteiger)
       dispm%q = gffq
       tmp_hq = gffh
+      secondary_q = gffs
    case(p_refq_goedecker)
       dispm%q = clsq
       tmp_hq = clsh
+      secondary_q = clss
    case(p_refq_gfn2xtb_gbsa_h2o)
       dispm%q = solq
       tmp_hq = solh
+      secondary_q = sols
    case default
       dispm%q = refq
       tmp_hq = refh
@@ -175,7 +155,7 @@
          is = refsys(j,ia)
          iz = zeff(is)
          sec_al = sscale(is)*secaiw(:,is) &
-            &  * zeta(g_a,gam(is)*g_c,secq(is)+iz,tmp_hq(j,ia)+iz)
+            &  * zeta(g_a,gam(is)*g_c,secondary_q(is)+iz,tmp_hq(j,ia)+iz)
          dispm%cn(j,ia) = refcovcn(j,ia)
          icn =nint(refcn(j,ia))
          cncount(icn) = cncount(icn) + 1
```

</details>

Upstream CTest results before/after the patch:

- `main`: **45/45 passed** both before and after.
- `v6.7.1`: **34/38 passed** both before and after. Under this toolchain, `xtb/gfn1`, `xtb/gfn2`, `xtb/hessian`, and `xtb/oniom` already fail in the unpatched build, with identical assertion values after the patch. I have not diagnosed those pre-existing failures.

This report concerns sequential initialization of new calculators; it does not establish general thread safety or cover interleavings with already-created calculators. A practical workaround for callers is to isolate GFN0 in a separate process, since deleting all C API handles does not reset these module arrays.
