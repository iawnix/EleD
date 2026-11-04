# QM8 first-excitation prediction pilot (0.9)

This experiment asks whether **ground-state xTB descriptors** improve a small
supervised prediction task. It does not calculate excited states with xTB and
does not identify the xTB HOMO–LUMO gap with an optical excitation energy.
The library remains an RDKit/libxtb adapter; sklearn/scipy are separate research
tools, not new calculation backends or runtime package dependencies.

## Data and fixed protocol

- Source: DeepChem's [gdb8 archive](https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/gdb8.tar.gz),
  SHA256 `e668b58f8dfba329bc20cd7e0eb3944e8a73a5d9db37160f686a8b31b20548a0`.
  Original study: Ramakrishnan et al., *Electronic spectra from TDDFT and machine
  learning in chemical space*, J. Chem. Phys. 143, 084111 (2015),
  [DOI:10.1063/1.4928757](https://doi.org/10.1063/1.4928757).
  See also the [DeepChem QM8 loader](https://github.com/deepchem/deepchem/blob/master/deepchem/molnet/load_function/qm8_datasets.py).
- Primary target: first vertical singlet excitation energy `E1-CC2`,
  RI-CC2/def2-TZVP. Source Hartree values are multiplied by 27.211386245988 for
  eV. The two sets of identically named PBE0 columns are different basis sets;
  they are not used. Property rows are joined to SDF titles by **gdb9 ID**,
  never by position after sanitization or filtering.
- Geometries: provided ground-state SDF coordinates, rounded to four decimal
  places, without reoptimization. Embedded source-header properties are never
  prediction features. This benchmark benefits from provided geometries and
  does not measure the cost or accuracy of generating them from SMILES.
- Of 21,786 source records, 39 fail RDKit sanitization, 1,413 fall outside our
  neutral/closed-shell protocol as represented by the parsed SDF, and 25 are
  repeated non-stereochemical canonical identities. These are representation
  and protocol exclusions, not claims that the original dataset physically
  contains 1,413 radicals/ions. No bond/charge/radical repair is attempted.
- From 20,309 eligible identities, select 1,024 by sorted SHA256 of a fixed
  prefix plus canonical non-stereochemical SMILES. No target-based selection.
  Deduplication conservatively retains only the first source stereoisomer for
  a given non-stereochemical identity. The full exclusions and selected IDs,
  targets and geometry hashes are frozen in the manifest.
- GFN2-xTB, 300 K electronic temperature, accuracy 1, independent single-point
  initial guesses; same extended libxtb as the package. All 1,024 rows complete
  with finite features. Any calculation/quality/missing-feature failure would
  be retained explicitly and excluded from **all three** feature groups, with
  no replacement. Software/environment errors abort the run.
- Same fixed ExtraTrees model for each group: 256 trees, minimum leaf size 2,
  all features considered at splits, seeds 17/29/43; no hyperparameter search.
  Median imputation is fitted only on training data. No target normalization,
  feature selection, or test-guided tuning.
- Random splits use unique chemical identities; scaffold splits use RDKit
  Bemis–Murcko scaffolds, with **all acyclic molecules in one empty-scaffold
  group**. GroupShuffleSplit holds out 20% of groups, not necessarily 20% of
  molecules. Each comparison uses identical training/test rows. Exact IDs and
  every test prediction are retained. The three seeds reuse molecules and are
  not three independent datasets.

## Ablation and observed results

The baseline includes 24 explicitly listed RDKit 2D descriptors and a 1,024-bit
Morgan fingerprint (radius 2, no chirality), rather than an artificially small
baseline. The second group adds nine basic xTB scalars: total energy, HOMO,
LUMO, gap, dipole norm, charge min/max/std and Wiberg pair sum. The final group
adds eight orbital distribution scalars: HOMO/LUMO IPR, participation ratio,
entropy, atom overlap and atom-centroid separation. These are frame-invariant
features; raw coordinate components are excluded.

Test MAE in eV, **mean ± sample standard deviation across three fixed splits**:

| Features | Random | Scaffold |
| --- | ---: | ---: |
| RDKit 2D + Morgan | 0.3469 ± 0.0121 | 0.3473 ± 0.0237 |
| + basic xTB | 0.2085 ± 0.0125 | 0.2105 ± 0.0217 |
| + orbital distribution | 0.2028 ± 0.0127 | 0.2093 ± 0.0246 |

Basic xTB features are useful in this particular pilot. Additional orbital
distribution features reduce mean MAE by only about 0.0057 eV under random
splits and 0.0012 eV under scaffold splits; one scaffold split becomes worse.
This does **not** establish a robust additional benefit or statistical
significance. IPR and participation ratio are reciprocals, so the eight added
columns do not represent eight independent pieces of information.

Recorded featurization cost on this local CPU run is about ten seconds for
1,024 molecules with one OpenMP/BLAS thread; each fixed forest fits/predicts in
roughly two to three seconds using two sklearn workers. Exact component times,
software versions and engine/parameter fingerprints are recorded with the
results. These timings exclude dataset download, sanitization, source geometry
generation, native compilation and dependency installation. They are not a
hardware-independent performance guarantee.

This is a small source-geometry pilot, not the full QM8 benchmark, a comparison
to state-of-the-art models, or a publication-ready generalization claim. Before
making a stronger scientific claim, freeze a larger protocol and evaluate
additional target families, models and independent datasets; test generated
geometries separately. The conformer example is a separate sensitivity study;
we have not shown that conformer averaging improves these predictions.

## Reproduce

Install the standard package/runtime first. Keep research dependencies separate:

```bash
.rdkit-eth/runtime/bin/python -m venv --system-site-packages build/benchmark-env
build/benchmark-env/bin/python -m pip install -r tools/benchmark-requirements.txt
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
  build/benchmark-env/bin/python tools/benchmark_qm8.py \
  --data build/qm8 --output build/qm8-new-run --size 1024
```

Use a new/empty output directory; the script refuses to overwrite prior results.
Downloaded data remain under `build/`, outside Git. Archived derived evidence is
under [`validation/0.9.0`](validation/0.9.0/README.md): the manifest and split IDs,
features/provenance, failures, all predictions and metrics. Small protocol tests
check ID joins with reordered properties and duplicate column names, checksums,
deduplication, reproducible splits and empty-scaffold isolation without any
external download. No benchmark dependency is added to `pyproject.toml`.
