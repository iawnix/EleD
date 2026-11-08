# Density-derived descriptors

`eled` can analyze the density returned by the extended `libxtb` build
without a second engine call.  Request `keep_overlap=True`; the result then
contains the same calculation's overlap matrix, orbital coefficients and
AO/shell maps.

The density is the spin-summed matrix `P = C diag(n) C^T`; its matrix entries
are basis-dependent and have no standalone physical unit. Every report carries
`trace(P S)` in electrons and the modeled electron count. Mulliken populations use the
diagonal of `P S`; Löwdin populations use the symmetric orthogonalized density.
Shell rows are reported in electrons for `s`, `p`, and `d`, with a dimensionless
`p_to_s_ratio`.  For GFN1/GFN2, a diagnostic charge can be reconstructed as
reference valence electrons minus the density population.  GFN0 shell
populations are density populations and must not be interpreted as its EEQ
charges. GFN1/GFN2 rows contain separate Mulliken and Löwdin reconstructed
charges (reference valence minus population), plus a Löwdin-compatible alias.

Mayer and Wiberg bond orders are listed for every RDKit bond. Mayer order is
withheld for open-shell results because this API has no spin-resolved density;
the row contains `spin_density_unavailable_for_open_shell`. HOMO/LUMO bond
overlap populations use the projector of the complete near-degenerate frontier
group, so orbital phase and degenerate rotations do not change them.

The π descriptors apply to aromatic or `sp2` atoms with a defined local plane.
The projected p-shell population is in electrons and π bond order is
dimensionless. Missing normals, incomplete p-shell maps, unsupported atoms and
open-shell Mayer values are explicit reasons, never zero-filled values. All
atom and bond values are invariant to rigid translation/rotation; normal sign
does not affect a π population or the reported bond-order magnitude.

The current validation target is atom-level EAS regioselectivity from the
public RegioML/RegioSQM data. QM8 molecular scalars are retained as smoke tests
but are not treated as evidence for atom-level utility. The baseline uses xTB
engine partial charges; this package does not label them CM5 charges. A frozen specification
is in [`protocol.json`](../protocol.json).

The implementation follows the standard density analyses of Mulliken (1955),
Löwdin (1950), Mayer (1983), and conceptual-DFT condensed Fukui indices
(Parr and Yang, 1989; Yang and Mortier, 1986). Frontier overlap populations
are an orbital-density feature, not an excitation or transition-density
calculation. The xTB/ML literature provides useful baselines for charges and
shell features, while shell-resolved π and frontier-bond decompositions are
treated here as explicitly testable extensions rather than established
experimental observables.

`DescriptorSensitivity(reference, gfn1=..., gfn2=..., temperature_300=...,
temperature_1000=...)` compares reports that already exist. It records
conformer standard deviation, the absolute GFN1/GFN2 difference, and the
absolute 300 K/1000 K difference independently. Omitted comparison reports
produce `None` with a missing flag and never trigger a new calculation.

The frozen dataset stage is reproducible from local public-data exports:

```bash
python tools/PrepareEASDataset.py regio-ml.jsonl \
  --source RegioML --output build/eas/manifest.json
python tools/run_eas_features.py build/eas/manifest.json \
  --xtb-library /path/to/libxtb.so --output build/eas/features.json
```

The preparation command validates atom-level labels and records excluded rows;
the feature command generates ETKDGv3 plus GFN2 geometry/descriptor reports and
retains every calculation failure. Neither command downloads data.
