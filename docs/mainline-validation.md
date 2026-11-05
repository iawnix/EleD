# Descriptor validation update (0.11)

These studies complete the two pending QM8 evaluations and add a second physical
target on nonoverlapping QM9 molecules. All failures and unfavorable comparisons
are retained. They do not introduce an embedding algorithm or an excited-state
xTB solver. Evidence is in [0.11 validation](validation/0.11.0/README.md).

## Correcting the historical QM8 study

The original 1,024-molecule cohort had 41 embedding failures. Independent stereo
checks reject nine of its 983 stored geometries, leaving **974/1,024** usable
rows. We reuse the original accepted coordinates and features, retain original
train/test memberships, and **refit all models**, including the source-geometry
models and 2D baseline, on this identical retained cohort. Removing only test
predictions would have left rejected training structures in fitted models.

The original evidence remains unchanged. Source/geometry checksums, exact
rejections, filtered splits, refitted predictions and conditional paired
bootstrap intervals are saved under `qm8-stereo-refit`.

Mean test MAE in eV across the three frozen seeds (the seeds are not independent
datasets; per-seed metrics and uncertainty are in the evidence):

| Features / geometry | Random | Scaffold |
|---|---:|---:|
| RDKit 2D + Morgan | 0.35244 | 0.34877 |
| Source + basic xTB | 0.21306 | 0.20769 |
| Source + basic xTB + orbital distributions | 0.20473 | 0.20725 |
| Generated + basic xTB | 0.23872 | 0.22990 |
| Generated + basic xTB + orbital distributions | 0.23076 | 0.22836 |

For generated-geometry scaffold splits, orbital-minus-basic MAE differences are
+0.00282, +0.00218 and -0.00961 eV. The first two conditional 95% intervals cross
zero; the third is [-0.02561, -0.00099] eV. One favorable split does not establish
a stable general benefit, especially without multiplicity adjustment.

## Current checked preparation workflow

`qm8-checked-preparation-v1` is a distinct protocol: eight requested conformers,
ETKDGv3 → ETKDGv4 → srETKDGv3, two seeds per initialization strategy, 300 embedding
iterations, MMFF94 up to 500 iterations, and stereo checks before/after relaxation.
Seeds are hashes of protocol and source ID. Stop at the first accepted pool and
choose its lowest-MMFF-energy conformer, with conformer ID breaking ties.
Source coordinates never initialize generation. The full protocol is frozen
before generation; model/feature definitions and original split memberships are
unchanged. Both preparation protocols use unmodified RDKit.

**999/1,024** molecules produce accepted geometries and finite xTB features;
**25** have no accepted conformer. All selected geometries, rejected attempts,
engine/parameter fingerprints and failure records are archived. Models for source
and generated geometries are refitted on the same 999 rows.

| Features / geometry | Random | Scaffold |
|---|---:|---:|
| RDKit 2D + Morgan | 0.35132 | 0.34951 |
| Source + basic xTB | 0.21436 | 0.20799 |
| Source + basic xTB + orbital distributions | 0.20557 | 0.20586 |
| Generated + basic xTB | 0.24435 | 0.23260 |
| Generated + basic xTB + orbital distributions | 0.23402 | 0.22945 |

The three generated-geometry scaffold orbital-minus-basic differences are
-0.00094, -0.00031 and -0.00819 eV; **all three conditional 95% intervals cross
zero**. Basic xTB features help in this pilot; additional orbital-distribution
features still lack a stable demonstrated benefit. Different cohort sizes and
retry budgets prevent interpreting the two tables as a controlled comparison
of embedding algorithms or equal-budget speed.

Both evaluations use the original 0.10 descriptor implementation from commit
`eb7eef9`, verified against the recorded Python-source hashes, and the original
engine/parameter hashes. The 0.11 SCC extension is not mixed into the frozen
source-feature comparison. Some jobs overlapped on the host; recorded wall times
are reproducibility metadata, not a performance benchmark.

## A second target: QM9 dipole moments

From the DeepChem QM9 archive (SHA256 recorded), require exactly nine heavy atoms,
neutral closed-shell input with explicit hydrogen coordinates, and unique
nonstereochemical identities. This leaves **86,278** eligible unique identities.
Select **256** by fixed SHA256 order of identity without using target values.
All selected molecules are outside QM8's at-most-eight-heavy-atom domain.

Compare GFN1/2 gas-phase dipole norms on the supplied SDF geometries directly
against QM9's B3LYP/6-31G(2df,p) dipole magnitude in Debye. There is no regression
model, parameter fitting, geometry optimization or experimental reference here.
SDF records join CSV labels by molecule ID, not row number.

| Method | Completed | MAE (Debye) | RMSE (Debye) | Mean signed error (Debye) |
|---|---:|---:|---:|---:|
| GFN1-xTB | 256/256 | 0.92144 | 1.07292 | +0.85395 |
| GFN2-xTB | 256/256 | 0.42859 | 0.49933 | +0.35288 |

GFN2-minus-GFN1 MAE is -0.49285 Debye; the paired molecule-bootstrap 95% interval
is [-0.53956, -0.44653]. This supports better GFN2 agreement **in this selected
cohort**. Both show positive bias. It does not validate every descriptor or
imply that 0.43 Debye is acceptable for every application. The two datasets
share the GDB family: this is nonoverlapping compound/target coverage, not a
broad external chemical-space benchmark. IP/EA calibration, experimental
reactivity-site validation, larger molecules and non-GDB datasets remain open.

QM9 reference: Ramakrishnan et al., *Quantum chemistry structures and properties
of 134 kilo molecules*, [DOI:10.1038/sdata.2014.22](https://doi.org/10.1038/sdata.2014.22).
The original QM8 protocol and its limitations are described in [the pilot](qm8-pilot.md).

## Reproduce

For QM8, use the original verified runtime/engine and separate sklearn environment
described by the frozen 0.9/0.10 evidence; mismatched engine hashes are rejected:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  build/benchmark-env/bin/python tools/benchmark_qm8_checked.py \
  --mode historical --output build/my-stereo-refit
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  build/benchmark-env/bin/python tools/benchmark_qm8_checked.py \
  --mode checked --output build/my-checked-preparation
```

QM9 accepts the downloaded, checksum-verified `gdb9.tar.gz` archive:

```bash
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  /path/to/0.11/runtime/bin/python tools/validate_qm9_dipoles.py \
  --archive build/qm9/gdb9.tar.gz --output build/my-qm9-dipoles
```

Each output directory must be new/empty. Existing evidence is never overwritten.
