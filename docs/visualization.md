# Visualization

EleD's visualization layer consumes an existing RDKit molecule and descriptor
report. It never starts an xTB calculation. RDKit provides the molecular graph,
atom order and 2D drawing; EleD supplies the electronic values and fragment
assignment.

```python
from eled import (
    CalcDescriptorReport,
    DrawDescriptorMap, PlotFragmentSummary, PlotFragmentInteractions,
    RenderFragmentAnalysis,
)

fragments = {
    "acetyl": [0, 1, 2],
    "phenoxy": [3, 4, 5, 6, 7, 8, 9],
}
report = CalcDescriptorReport(mol, result, fragments=fragments)

molecule_svg = DrawDescriptorMap(mol, report, atom_field="q_e")
summary_svg = PlotFragmentSummary(report, field="q_e")
network_svg = PlotFragmentInteractions(report)
html_fragment = RenderFragmentAnalysis(mol, report)
```

`DrawDescriptorMap` uses RDKit's `rdMolDraw2D` and keeps the original RDKit atom
indices. Explicit hydrogens are omitted from the map by default so their labels
do not crowd the heavy-atom values; the report still contains their descriptors.
Use `include_hydrogens=True` for a hydrogen-level map. Set
`add_atom_indices=True` to include the original index in each value label.
`PlotFragmentSummary` renders signed fragment values such as charge, Fukui
indices or local softness. `PlotFragmentInteractions` uses the explicit
fragment-pair Wiberg values from `report["fragments"]["pairs"]`; it does not add
the two directional external sums a second time.

The CamelCase functions are the documented RDKit-style public entry points. The
core package does not add
matplotlib, Plotly or NetworkX dependencies; the functions return SVG/HTML
strings that can be embedded in notebooks, reports or a web view.

Fragment assignment is explicit. Unassigned atoms remain visible in the
molecule view and are reported by `BuildFragmentVisualData`. The current layer
does not draw real-space electron-density or orbital isosurfaces because an
xTB result contains AO matrices, not a volumetric grid.
