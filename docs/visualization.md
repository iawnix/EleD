# Visualization

EleD's visualization layer consumes an existing RDKit molecule and either a
descriptor report (`CalcDescriptorReport`, `RunReactivity`,
`CalcDensityDescriptors`) or an `XTBResults` directly. It never starts an xTB
calculation. RDKit provides the 2D layout, atom/bond indexing and molecular
rendering; EleD supplies the electronic values, continuous colormaps and
fragment summaries.

## Molecular and bond maps

`DrawDescriptorMap` works with or without fragment partitions. When a molecule
carries 3D coordinates, EleD generates clean 2D drawing coordinates on an
internal copy without modifying the input molecule. By default, explicit
hydrogens are omitted while preserving exact original atom/bond indices (even
when RDKit retains isotopic, charged or stereochemical hydrogens), and numeric
values are rendered as RDKit `atomNote` / `bondNote` annotations so element
symbols remain readable.

```python
from eled import DrawDescriptorMap

# Single-molecule charge heatmap directly from XTBResults or a report (no fragments required)
charge_svg = DrawDescriptorMap(mol, result, atomField="q_e")

# Bond-order or density-derived fields (e.g. wbo, pi_bo, coop_e, frontier_bonds.homo.mul_ovlp_e)
bond_svg = DrawDescriptorMap(mol, report, bondField="wbo", colorBy="value")

# Custom RDKit MolDrawOptions callback
styled_svg = DrawDescriptorMap(
    mol,
    report,
    atomField="cm5_q_e",
    colorBy="value",
    labelStyle="note",      # "note" (default), "replace", or "none"
    addAtomIndices=True,
    includeHydrogens=False,
    showColorbar=True,
    drawOptions=lambda opts: setattr(opts, "bondLineWidth", 2.0),
)
```

Supported field names include top-level atomic fields (`q_e`, `cm5_q_e`,
`wbo_sum`, `f_plus`, `f_minus`, `f_radical`, `dual`, `s_loc_1ev`,
`omega_loc_ev`, `n_loc_1ev`), nested frontier/shell/π atomic fields
(`homo_population`, `lumo_population`, `mul_p_e`, `low_p_e`, `low_ps`,
`pi_n_e`, `pi_q_e`), and bond fields (`wbo`, `pi_bo`,
`frontier_bonds.homo.mul_ovlp_e`, `frontier_bonds.lumo.mul_ovlp_e`). Dotted
section paths (such as `pi.pi_q_e` or `pi.pi_bo`) are also accepted.
When `report` includes provenance metadata, `BuildFragmentVisualData` and
`DrawDescriptorMap` verify atom count and ordered structure hash against `mol`.

`colorBy` controls highlighting:

- `"auto"` (default): colors by fragments when fragments are present and no
  `bondField` is requested; otherwise uses continuous value coloring (`"value"`).
- `"value"`: uses a zero-symmetric diverging red/blue scale for signed fields
  (`q_e`, `cm5_q_e`, `dual`, `pi_q_e`, `mul_ovlp_e`, etc.) and a sequential blue
  scale for nonnegative fields, with an optional bottom colorbar (`showColorbar`).
- `"fragment"`: colors atoms by fragment assignment.
- `"none"`: disables highlight fills.

## Fragment summaries and combined HTML view

```python
from eled import (
    CalcDescriptorReport,
    DrawDescriptorMap,
    PlotFragmentSummary,
    PlotFragmentInteractions,
    RenderFragmentAnalysis,
)

fragments = {
    "acetyl": [0, 1, 2],
    "phenoxy": [3, 4, 5, 6, 7, 8, 9],
}
report = CalcDescriptorReport(mol, result, fragments=fragments)

molecule_svg = DrawDescriptorMap(mol, report, atomField="q_e")
summary_svg = PlotFragmentSummary(report, field="q_e")
network_svg = PlotFragmentInteractions(report)
html_fragment = RenderFragmentAnalysis(mol, report, atomField="q_e")
```

- `PlotFragmentSummary` renders signed horizontal bars for any numeric fragment
  field (such as `q_e`, `cm5_q_e`, `f_plus`, `dual` or `s_loc_1ev`), sizing the
  label column dynamically and marking unavailable values explicitly as
  `missing`.
- `PlotFragmentInteractions` renders fragment-pair Wiberg sums from
  `report["fragments"]["pairs"]` using dynamically sized pill nodes without
  double-counting directional sums.
- `RenderFragmentAnalysis` combines the molecular map with fragment views when
  fragments are present, links `summaryField` to `atomField` by default, and
  includes light/dark-mode CSS matching the embedded SVGs.

The core package does not add matplotlib, Plotly or NetworkX dependencies; all
functions return self-contained SVG/HTML strings suitable for notebooks,
reports or web views.

Fragment assignment is explicit. Unassigned atoms remain visible in the
molecule view and are reported by `BuildFragmentVisualData`. The current layer
does not draw real-space electron-density or orbital isosurfaces because an
xTB result contains AO matrices, not a volumetric grid.
