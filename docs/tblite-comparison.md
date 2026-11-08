# Optional tblite backend

> Historical 0.1–0.3 validation evidence. Version 0.4 has no tblite backend.
> Reproduce the commands below with the preserved `dist/0.3.0` source archive;
> comparison sources are archived in `tools/diagnostics/archive/tblite-0.3/`.
> Current population analysis uses the [native basis extension](xtb-patched-build.md).

`eled` can use [tblite](https://github.com/tblite/tblite) as an explicit,
optional comparison backend:

```bash
python -m pip install -e '.[tblite]'  # from the source checkout
```

This installs the tblite Python wheel and its two runtime Python dependencies,
NumPy and CFFI. It does **not** install ASE. The base `eled` dependency
set is unchanged, and the default backend remains native libxtb. A
current Linux CPython 3.11 tblite wheel is about 7.2 MB compressed and about
18 MB unpacked including its bundled native library; actual size varies by
platform and Python version. This is an optional runtime cost, not a cost paid
by users of the native backend.

The adapter is deliberately a thin RDKit boundary layer. It converts atom
numbers and RDKit Angstrom coordinates into tblite's NumPy structure arrays in
Bohr, sets the same numerical accuracy, 250 iteration limit and 300 K
electronic temperature used by our native test path, and normalizes the
results into the existing `XTBResult` fields:

- total energy: Hartree;
- orbital energies: eV;
- occupations: electrons;
- charges: elementary charge;
- dipole: e·Bohr;
- bond orders: dimensionless atom-by-atom matrix;
- orbital coefficients: `[AO, MO]`.

tblite's result interface also exposes density, overlap and Hamiltonian
matrices, localized orbitals, orbital-to-shell maps and shell-to-atom maps.
The public `keep_overlap=True` option now retains overlap in
`XTBResult.overlap_matrix`, with a typed `ao_atom_indices` map. Other integral
matrices remain outside the public result. The diagnostic map arrays are
also retained in `result.raw_json["orbital_map"]`, `shell_map`, and
`angular_momenta`. This makes the comparison useful for deciding whether a
future descriptor layer should use tblite's richer basis metadata.

## What is compared

The optional tests run the same RDKit molecule and 3D conformer through:

```text
RDKit → tblite Calculator
RDKit → libxtb C API
```

They compare GFN1 and GFN2 on water, ethanol, benzene and ammonium. The tests
also check atom ordering, charge/spin defaults, units, coefficient axes,
mapping metadata, cache identity and explicit backend selection. The compiled
module is also compared against the ctypes adapter, and a frozen water energy
from the earlier executable comparison remains a numerical regression check.

These original tests compare the ctypes adapter against tblite. Direct
compiled `rdXTBTools` comparisons now use 33 frozen conformers for each method,
with per-case reports, source-verified s/p/d row conventions, d-row phases,
and overlap-metric sign/subspace diagnostics. See
[orbital alignment and reproducible results](orbital-alignment.md).

The earlier exploratory GFN2 maxima came from water, ethanol, benzene, acetic
acid and ammonium; they were not maxima over a broad GFN1/GFN2 benchmark.
The larger comparison finds some differences outside the original tolerances
at default accuracy 1. Tightening both engines to accuracy 0.01 brings all 66
cases from the local source build within the recorded limits. Official 6.7.1
retains two GFN1 exceptions. The reports preserve both outcomes.

Coefficient rows differ between engines, including a d-basis phase. Merely
comparing absolute coefficients is insufficient. Neither raw nor diagnostically
aligned coefficient squares are automatically atom populations in a
nonorthogonal basis.

## Current scope

The optional backend supports gas-phase, closed-shell GFN1-xTB and GFN2-xTB
single points. tblite 0.8.0 does not expose GFN0-xTB, so GFN0 remains on the
verified native libxtb path. The current adapter rejects nonzero `uhf` because `XTBResult` does not
represent separate spin channels. This is a conservative adapter restriction,
not a statement that tblite cannot calculate radicals; use native libxtb for
radicals in this package. The
backend does not currently add tblite's ALPB/GBSA, gradient, periodic or
density-matrix functionality to the public `eled` API. ASE is not part
of this path.

This separation is intentional: tblite is useful for validating numerical
behavior and for a future richer basis-aware descriptor layer, but it is not a
mandatory replacement for the native libxtb adapter.

Version 0.3 uses the same-calculation overlap for [Löwdin/Mulliken and frontier
subspace descriptors](descriptors.md). This is an explicit optional feature;
the native backend and its dependency set are unchanged.
