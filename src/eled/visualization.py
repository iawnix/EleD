"""RDKit-native visualizations for EleD descriptor reports.

The functions in this module only consume an RDKit molecule and an existing
report. They never call xTB. RDKit draws the molecular graph; the summary and
interaction views are small dependency-free SVGs so the core package does not
need matplotlib, Plotly or NetworkX.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from html import escape
import math


_COLORS = (
    (0.20, 0.48, 0.79),
    (0.84, 0.39, 0.15),
    (0.28, 0.58, 0.40),
    (0.61, 0.39, 0.73),
    (0.75, 0.58, 0.16),
    (0.25, 0.62, 0.64),
)


def _finite_number(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _svg_color(rgb):
    return '#{:02x}{:02x}{:02x}'.format(*(round(255 * value) for value in rgb))


def _unit(field):
    if field.endswith('_1ev'):
        return '1/eV'
    if field.endswith('_ebohr'):
        return 'e*Bohr'
    if field.endswith('_ev'):
        return 'eV'
    if field.endswith('_e') or field in {'charge', 'q_e'}:
        return 'e'
    return 'dimensionless'


def _field_label(field):
    return {
        'q_e': 'charge',
        'f_plus': 'Fukui f+',
        'f_minus': 'Fukui f−',
        's_loc_1ev': 'local softness',
        'omega_loc_ev': 'local electrophilicity',
        'n_loc_1ev': 'local nucleophilicity',
        'homo_population': 'HOMO population',
        'lumo_population': 'LUMO population',
        'wbo_sum': 'Wiberg sum',
    }.get(field, field.replace('_', ' '))


def _report_section(report):
    if not isinstance(report, Mapping):
        raise TypeError('report must be a mapping returned by an EleD analysis')
    section = report.get('fragments', report)
    if not isinstance(section, Mapping):
        raise ValueError('report does not contain a fragment section')
    raw = section.get('fragments', section)
    if not isinstance(raw, Mapping):
        raise ValueError('fragment section must map names to records')
    records = []
    for name, record in raw.items():
        if not isinstance(record, Mapping):
            raise ValueError(f'fragment {name!r} is not a mapping')
        atoms = record.get('atom_indices', ())
        if isinstance(atoms, (str, bytes)):
            raise ValueError(f'fragment {name!r} atom_indices must be a sequence')
        try:
            atom_indices = [int(index) for index in atoms]
        except (TypeError, ValueError):
            raise ValueError(f'fragment {name!r} atom_indices must be integer indices') from None
        if len(set(atom_indices)) != len(atom_indices) or any(index < 0 for index in atom_indices):
            raise ValueError(f'fragment {name!r} atom_indices must be unique nonnegative indices')
        row = dict(record)
        row['name'] = str(name)
        row['atom_indices'] = atom_indices
        records.append(row)
    pairs = section.get('pairs', ())
    if not isinstance(pairs, Sequence) or isinstance(pairs, (str, bytes)):
        raise ValueError('fragment pairs must be a sequence')
    return records, [dict(pair) for pair in pairs], section


def _atom_rows(report):
    rows = report.get('atoms', ()) if isinstance(report, Mapping) else ()
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        return {}
    return {int(row['atom_index']): row for row in rows
            if isinstance(row, Mapping) and 'atom_index' in row}


def BuildFragmentVisualData(mol, report, *, atom_field=None):
    """Build a JSON-compatible mapping for molecular and fragment views.

    The input fragment assignment is explicit and disjoint, as required by the
    descriptor API. Atom indices remain RDKit indices. No calculation is run.
    """
    if not hasattr(mol, 'GetNumAtoms'):
        raise TypeError('mol must be an RDKit molecule')
    records, pairs, section = _report_section(report)
    natoms = int(mol.GetNumAtoms())
    assignments = {}
    for record in records:
        for index in record['atom_indices']:
            if index >= natoms:
                raise ValueError(f'fragment atom index {index} is outside the molecule')
            if index in assignments:
                raise ValueError(f'atom {index} belongs to more than one fragment')
            assignments[index] = record['name']
    atom_rows = _atom_rows(report)
    if atom_field is None:
        for candidate in ('q_e', 'f_plus', 's_loc_1ev', 'homo_population'):
            if any(candidate in row for row in atom_rows.values()):
                atom_field = candidate
                break
    atoms = []
    for index in range(natoms):
        row = atom_rows.get(index, {})
        atoms.append({'atom_index': index, 'fragment': assignments.get(index),
                      'value': _finite_number(row.get(atom_field)) if atom_field else None})
    names = {record['name'] for record in records}
    normalized_pairs = []
    for pair in pairs:
        first, second = pair.get('fragment_a'), pair.get('fragment_b')
        if first not in names or second not in names:
            continue
        value = pair.get('wbo', pair.get('value'))
        normalized_pairs.append({'fragment_a': str(first), 'fragment_b': str(second),
                                 'value': _finite_number(value),
                                 'wbo': _finite_number(pair.get('wbo'))})
    return {
        'atom_field': atom_field,
        'atom_unit': _unit(atom_field or ''),
        'atoms': atoms,
        'fragments': records,
        'pairs': normalized_pairs,
        'unassigned_atom_indices': sorted(set(range(natoms)) - set(assignments)),
        'units': dict(report.get('units', {})) if isinstance(report, Mapping) else {},
        'source': section.get('source', 'EleD descriptor report'),
    }


def _strip_xml_header(svg):
    if svg.lstrip().startswith('<?xml') and '?>' in svg:
        return svg.split('?>', 1)[1].lstrip('\n')
    return svg


def DrawDescriptorMap(mol, report, *, atom_field=None, width=520, height=320,
                     add_atom_indices=True):
    """Return an SVG RDKit drawing with fragments and atom values highlighted."""
    from rdkit.Chem.Draw import rdMolDraw2D

    data = BuildFragmentVisualData(mol, report, atom_field=atom_field)
    records = data['fragments']
    assignments = {row['atom_index']: row['fragment'] for row in data['atoms']}
    colors = {record['name']: _COLORS[index % len(_COLORS)]
              for index, record in enumerate(records)}
    highlight_atoms = [index for index, name in assignments.items() if name is not None]
    highlight_colors = {index: colors[name] for index, name in assignments.items() if name is not None}
    drawer = rdMolDraw2D.MolDraw2DSVG(int(width), int(height))
    options = drawer.drawOptions()
    options.addAtomIndices = bool(add_atom_indices)
    if data['atom_field']:
        for row in data['atoms']:
            if row['value'] is not None:
                options.atomLabels[int(row['atom_index'])] = (
                    f"{row['atom_index']}:{row['value']:+.2f}"
                    if not add_atom_indices else f"{row['value']:+.2f}"
                )
    drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms,
                        highlightAtomColors=highlight_colors)
    drawer.FinishDrawing()
    return _strip_xml_header(drawer.GetDrawingText())


def PlotFragmentSummary(report, *, field='q_e', width=620, row_height=38):
    """Return an SVG signed bar chart for one fragment field."""
    data = BuildFragmentVisualData(_report_atom_count(report), report)
    records = data['fragments']
    values = [_finite_number(row.get(field)) for row in records]
    finite = [value for value in values if value is not None]
    height = max(130, 70 + len(records) * row_height)
    title = escape(f"{_field_label(field)} ({_unit(field)})")
    if not finite:
        return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
                f'role="img" aria-label="No {title} values available"><title>No {title} values available</title>'
                f'<text x="20" y="48" fill="currentColor">No values available</text></svg>')
    maximum = max(abs(value) for value in finite) or 1.0
    left, chart_width, center = 160, width - 210, 160 + (width - 210) / 2
    body = [f'<line x1="{center:.1f}" y1="55" x2="{center:.1f}" y2="{height - 18}" stroke="currentColor" opacity=".35"/>']
    for index, (record, value) in enumerate(zip(records, values)):
        y = 72 + index * row_height
        name = escape(str(record['name']))
        body.append(f'<text x="12" y="{y + 5}" fill="currentColor">{name}</text>')
        if value is not None:
            span = abs(value) / maximum * chart_width / 2
            x = center if value >= 0 else center - span
            fill = 'var(--viz-series-1)' if value >= 0 else 'var(--viz-series-2)'
            body.append(f'<rect x="{x:.1f}" y="{y - 8}" width="{span:.1f}" height="16" fill="{fill}" rx="3"/>')
            body.append(f'<text x="{width - 12}" y="{y + 5}" text-anchor="end" fill="currentColor">{value:+.3f}</text>')
        else:
            body.append(f'<text x="{width - 12}" y="{y + 5}" text-anchor="end" fill="currentColor">missing</text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'role="img" aria-label="Fragment {title} summary"><title>Fragment {title} summary</title>'
            f'<text x="12" y="25" fill="currentColor" font-size="16">{title}</text>{"".join(body)}</svg>')


def PlotFragmentInteractions(report, *, width=620, height=280):
    """Return an SVG network of fragment-pair Wiberg couplings."""
    data = BuildFragmentVisualData(_report_atom_count(report), report)
    records, pairs = data['fragments'], data['pairs']
    names = [record['name'] for record in records]
    if not names:
        return f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img"><text x="20" y="40" fill="currentColor">No fragments available</text></svg>'
    positions = {}
    if len(names) == 1:
        positions[names[0]] = (width / 2, height / 2)
    else:
        for index, name in enumerate(names):
            positions[name] = (width / 2 + (width / 2 - 90) * math.cos(2 * math.pi * index / len(names)),
                               height / 2 + (height / 2 - 70) * math.sin(2 * math.pi * index / len(names)))
    finite = [pair['value'] for pair in pairs if pair['value'] is not None]
    maximum = max((abs(value) for value in finite), default=1.0)
    body = []
    for pair in pairs:
        if pair['value'] is None:
            continue
        x1, y1 = positions[pair['fragment_a']]
        x2, y2 = positions[pair['fragment_b']]
        stroke_width = 2 + 10 * abs(pair['value']) / maximum
        body.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="var(--viz-series-3)" stroke-width="{stroke_width:.1f}" stroke-linecap="round" opacity=".75"/>')
        body.append(f'<text x="{(x1+x2)/2:.1f}" y="{(y1+y2)/2-10:.1f}" text-anchor="middle" fill="currentColor" font-size="12">{pair["value"]:.2f}</text>')
    for index, name in enumerate(names):
        x, y = positions[name]
        fill = _svg_color(_COLORS[index % len(_COLORS)])
        body.append(f'<circle cx="{x:.1f}" cy="{y:.1f}" r="38" fill="{fill}" opacity=".88"/>')
        body.append(f'<text x="{x:.1f}" y="{y+5:.1f}" text-anchor="middle" fill="white" font-size="13">{escape(name)}</text>')
    return (f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'role="img" aria-label="Fragment interaction network"><title>Fragment interaction network</title>'
            f'{"".join(body)}</svg>')


def _report_atom_count(report):
    """Make a tiny molecule-like object for report-only SVG functions."""
    class _Count:
        def __init__(self, count): self._count = count
        def GetNumAtoms(self): return self._count
    rows = report.get('atoms', ()) if isinstance(report, Mapping) else ()
    indices = [int(row['atom_index']) for row in rows if isinstance(row, Mapping) and 'atom_index' in row]
    if not indices:
        records, _, _ = _report_section(report)
        indices = [index for record in records for index in record['atom_indices']]
    return _Count(max(indices, default=-1) + 1)


def RenderFragmentAnalysis(mol, report, *, atom_field=None, width=620):
    """Return a self-contained HTML fragment with the three EleD views."""
    data = BuildFragmentVisualData(mol, report, atom_field=atom_field)
    title = escape(_field_label(data['atom_field'] or 'descriptor'))
    molecule = DrawDescriptorMap(mol, report, atom_field=atom_field, width=width, height=320)
    summary = PlotFragmentSummary(report, field='q_e', width=width)
    interactions = PlotFragmentInteractions(report, width=width)
    return f'''<div class="eled-fragment-analysis">
<style>
.eled-fragment-analysis {{ font-family: sans-serif; color: #20242a; display: grid; gap: 12px; max-width: {width}px; }}
.eled-fragment-analysis section {{ border: 1px solid #d8dde5; border-radius: 8px; padding: 12px; }}
.eled-fragment-analysis h2 {{ margin: 0 0 8px; font-size: 16px; }}
.eled-fragment-analysis svg {{ max-width: 100%; height: auto; }}
.eled-fragment-analysis .note {{ color: #59636f; font-size: 12px; }}
</style>
<section><h2>RDKit molecule with EleD {_field_label(data['atom_field'] or 'descriptor')}</h2>{molecule}</section>
<section><h2>Fragment summary</h2>{summary}</section>
<section><h2>Fragment interaction network</h2>{interactions}<div class="note">Edges are fragment-pair Wiberg sums.</div></section>
</div>'''


__all__ = [
    'BuildFragmentVisualData', 'DrawDescriptorMap', 'PlotFragmentSummary',
    'PlotFragmentInteractions', 'RenderFragmentAnalysis',
]
