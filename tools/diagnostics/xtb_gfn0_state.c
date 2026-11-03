/* Standalone xTB 6.7.1 C API reproducer; no RDKit or Python adapter required.
 * Usage: repro PARAMETER_DIRECTORY [0|1|2] [--inspect-restore]
 *              [--run-intermediate] [--water] [--expect-clean]
 * Every calculation owns fresh handles and deletes them before returning.
 * The optional diagnostic accesses GNU Fortran internals to test causality.
 * It is NOT a supported API or a production workaround.
 */
#include <xtb.h>
#include <dlfcn.h>
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
        fprintf(stderr, "usage: %s PARAMETER_DIRECTORY [0|1|2] [--inspect-restore] "
                "[--run-intermediate] [--water] [--expect-clean]\n", argv[0]);
        return 2;
    }
    int method = 0, inspect = 0, run_intermediate = 0, expect_clean = 0;
    int first_flag = 2;
    if (argc >= 3 && argv[2][0] != '-') {
        if (strlen(argv[2]) != 1 || argv[2][0] < '0' || argv[2][0] > '2') return 2;
        method = argv[2][0] - '0';
        first_flag = 3;
    }
    for (int i = first_flag; i < argc; ++i) {
        if (strcmp(argv[i], "--inspect-restore") == 0) inspect = 1;
        else if (strcmp(argv[i], "--run-intermediate") == 0) run_intermediate = 1;
        else if (strcmp(argv[i], "--expect-clean") == 0) expect_clean = 1;
        else if (strcmp(argv[i], "--water") == 0) {
            atom_count = 3;
            numbers = water_numbers;
            xyz = water_xyz;
        } else return 2;
    }
    /* Source declarations in xtb_disp_dftd4param: refq(7,118), refh(7,118). */
    double *refq = NULL, *refh = NULL;
    double initial_q[7 * 118], initial_h[7 * 118];
    if (inspect) {
        refq = dlsym(RTLD_DEFAULT, "__xtb_disp_dftd4param_MOD_refq");
        refh = dlsym(RTLD_DEFAULT, "__xtb_disp_dftd4param_MOD_refh");
        if (!refq || !refh) {
            fprintf(stderr, "diagnostic requires GNU Fortran xTB 6.7.1 symbols\n");
            return 2;
        }
        memcpy(initial_q, refq, sizeof(initial_q));
        memcpy(initial_h, refh, sizeof(initial_h));
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
    if (inspect) {
        printf("{\"stage\":\"after_parameter_load\",\"refq_max_change\":%.16g,"
               "\"refh_max_change\":%.16g}\n",
               maxdiff(refq, initial_q, 7 * 118), maxdiff(refh, initial_h, 7 * 118));
    }
    struct result after = calculate(argv[1], 2, 0);
    report("gfn2_after_parameter_load", before, after);
    int contaminated = fabs(after.energy - before.energy) > 1e-10 ||
        maxdiff(after.orbitals, before.orbitals, before.nao) * 27.211386245988 > 1e-8;
    if (inspect) {
        memcpy(refq, initial_q, sizeof(initial_q));
        memcpy(refh, initial_h, sizeof(initial_h));
        struct result restored = calculate(argv[1], 2, 0);
        report("gfn2_after_restoring_only_refq_refh", before, restored);
        free(restored.orbitals);
    }
    free(before.orbitals);
    free(after.orbitals);
    return expect_clean && contaminated ? 3 : 0;
}
