"""RDKit-native visualizations for EleD descriptor reports and results.

The functions in this module only consume an RDKit molecule and an existing
report or ``XTBResults``. They never invoke xTB.
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

_SIGNED_FIELDS = frozenset({
    'q_e',
    'cm5_q_e',
    'dual',
    'pi_q_e',
    'mul_qrec_e',
    'low_qrec_e',
    'mul_ovlp_e',
    'homo_mul_ovlp_e',
    'lumo_mul_ovlp_e',
    'frontier_bonds.homo.mul_ovlp_e',
    'frontier_bonds.lumo.mul_ovlp_e',
})

_FIELD_LABELS = {
    'q_e': 'charge',
    'cm5_q_e': 'CM5 charge',
    'wbo_sum': 'Wiberg sum',
    'wbo': 'Wiberg bond order',
    'f_plus': 'Fukui f+',
    'f_minus': 'Fukui f−',
    'f_radical': 'Fukui f0',
    'dual': 'dual descriptor',
    's_loc_1ev': 'local softness',
    'omega_loc_ev': 'local electrophilicity',
    'n_loc_1ev': 'local nucleophilicity',
    'homo_population': 'HOMO population',
    'lumo_population': 'LUMO population',
    'frontier.homo.atom_pop': 'HOMO population',
    'frontier.lumo.atom_pop': 'LUMO population',
    'pi_n_e': 'π electron population',
    'pi_q_e': 'π charge',
    'pi_bo': 'π bond order',
    'mul_s_e': 'Mulliken s population',
    'mul_p_e': 'Mulliken p population',
    'mul_d_e': 'Mulliken d population',
    'low_s_e': 'Löwdin s population',
    'low_p_e': 'Löwdin p population',
    'low_d_e': 'Löwdin d population',
    'mul_ps': 'Mulliken p/s ratio',
    'low_ps': 'Löwdin p/s ratio',
    'mul_qrec_e': 'Mulliken reconstructed charge',
    'low_qrec_e': 'Löwdin reconstructed charge',
    'homo_mul_ovlp_e': 'HOMO overlap population',
    'lumo_mul_ovlp_e': 'LUMO overlap population',
    'frontier_bonds.homo.mul_ovlp_e': 'HOMO overlap population',
    'frontier_bonds.lumo.mul_ovlp_e': 'LUMO overlap population',
}


def _finite_number(value):
    if value is None or isinstance(value, bool):
        return None
    try:
        value = float(value)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def _svg_color(rgb):
    return '#{:02x}{:02x}{:02x}'.format(*(round(255 * max(0.0, min(1.0, value))) for value in rgb))


def _lerp_rgb(start, end, t):
    t = max(0.0, min(1.0, float(t)))
    return tuple(a + (b - a) * t for a, b in zip(start, end))


_DIV_NEG = (0.89, 0.44, 0.39)
_DIV_MID = (0.97, 0.98, 0.99)
_DIV_POS = (0.36, 0.62, 0.87)
_SEQ_LOW = (0.93, 0.96, 0.99)
_SEQ_HIGH = (0.28, 0.58, 0.85)


def _value_rgb(value, vmin, vmax, *, diverging):
    if diverging:
        span = max(abs(vmin), abs(vmax)) or 1.0
        ratio = max(-1.0, min(1.0, value / span))
        return _lerp_rgb(_DIV_MID, _DIV_POS, ratio) if ratio >= 0 else _lerp_rgb(_DIV_MID, _DIV_NEG, -ratio)
    span = (vmax - vmin) if vmax > vmin else (abs(vmax) or 1.0)
    ratio = (value - vmin) / span if vmax > vmin else (value / span)
    return _lerp_rgb(_SEQ_LOW, _SEQ_HIGH, ratio)


def _pastel_rgb(rgb, factor=0.55):
    return _lerp_rgb(rgb, (1.0, 1.0, 1.0), factor)


_SVG_SERIES = tuple(_svg_color(rgb) for rgb in _COLORS[:3])
_SVG_STYLE = '''<style>
.eled-viz-bg { fill: #ffffff; }
.eled-viz-text { fill: #172033; font-family: system-ui, -apple-system, sans-serif; }
.eled-viz-sub { fill: #526174; font-family: system-ui, -apple-system, sans-serif; }
.eled-viz-muted { stroke: #94a3b8; }
@media (prefers-color-scheme: dark) {
  .eled-viz-bg { fill: #0f172a; }
  .eled-viz-text { fill: #f8fafc; }
  .eled-viz-sub { fill: #94a3b8; }
  .eled-viz-muted { stroke: #64748b; }
}
</style>'''


def _unit_from_name(field):
    if not field:
        return 'dimensionless'
    leaf = field.rsplit('.', 1)[-1]
    if leaf.endswith('_1ev'):
        return '1/eV'
    if leaf.endswith('_ebohr'):
        return 'e*Bohr'
    if leaf.endswith('_debye'):
        return 'Debye'
    if leaf.endswith('_ev'):
        return 'eV'
    if leaf.endswith('_ha'):
        return 'Hartree'
    if leaf.endswith('_a'):
        return 'Angstrom'
    if leaf.endswith('_e') or leaf in {'charge', 'q_e'}:
        return 'e'
    return 'dimensionless'


def _resolve_unit(report, field):
    if not field:
        return 'dimensionless'
    leaf = field.rsplit('.', 1)[-1]
    if isinstance(report, Mapping):
        unit_maps = [
            report.get('atom_units'),
            report.get('bond_units'),
            report.get('units'),
        ]
        for section_name in ('shell_populations', 'frontier_bonds', 'pi'):
            section = report.get(section_name)
            if isinstance(section, Mapping):
                unit_maps.append(section.get('units'))
        for mapping in unit_maps:
            if isinstance(mapping, Mapping):
                if field in mapping and isinstance(mapping[field], str):
                    return mapping[field]
                if leaf in mapping and isinstance(mapping[leaf], str):
                    return mapping[leaf]
    return _unit_from_name(field)


def _field_label(field):
    if not field:
        return 'descriptor'
    if field in _FIELD_LABELS:
        return _FIELD_LABELS[field]
    leaf = field.rsplit('.', 1)[-1]
    if leaf in _FIELD_LABELS:
        return _FIELD_LABELS[leaf]
    return leaf.replace('_', ' ')


def _is_xtb_result(obj):
    from .result import XTBResults
    return isinstance(obj, XTBResults) or (
        hasattr(obj, 'GetAtomicCharges') and hasattr(obj, 'totalEnergyHartree')
    )


def _normalize_report(mol, report):
    if _is_xtb_result(report):
        from .reports import CalcDescriptorReport
        return CalcDescriptorReport(mol, report)
    if not isinstance(report, Mapping):
        raise TypeError('report must be a mapping returned by an EleD analysis or an XTBResults')
    return report


def _validate_provenance(mol, report):
    if not isinstance(report, Mapping):
        return
    provenance = report.get('provenance')
    if not isinstance(provenance, Mapping):
        return
    natoms = int(mol.GetNumAtoms())
    atom_count = provenance.get('atom_count')
    if atom_count is not None and int(atom_count) != natoms:
        raise ValueError('report atom count disagrees with the indexed molecule')
    expected_hash = provenance.get('ordered_structure_sha256')
    if expected_hash:
        from .input import _ordered_structure_hash
        if _ordered_structure_hash(mol) != expected_hash:
            raise ValueError('molecule structure or atom order does not match report provenance')


def _extract_fragments(report, *, required=False):
    if not isinstance(report, Mapping):
        raise TypeError('report must be a mapping returned by an EleD analysis')
    if 'fragments' in report:
        section = report['fragments']
        if section is None:
            if required:
                raise ValueError('report does not contain a fragment section')
            return [], [], {}
    else:
        section = report
        # A plain report with atoms/bonds/values and no fragment records has no fragments.
        if any(key in report for key in ('atoms', 'bonds', 'values', 'report_type', 'provenance')):
            if required:
                raise ValueError('report does not contain a fragment section')
            return [], [], {}
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
        if isinstance(atoms, (str, bytes)) or not isinstance(atoms, Sequence):
            raise ValueError(f'fragment {name!r} atom_indices must be a sequence')
        atom_indices = []
        for index in atoms:
            if isinstance(index, bool):
                raise ValueError(f'fragment {name!r} atom_indices must be integer indices')
            try:
                atom_indices.append(int(index))
            except (TypeError, ValueError):
                raise ValueError(f'fragment {name!r} atom_indices must be integer indices') from None
        if len(set(atom_indices)) != len(atom_indices) or any(index < 0 for index in atom_indices):
            raise ValueError(f'fragment {name!r} atom_indices must be unique nonnegative indices')
        row = dict(record)
        row['name'] = str(name)
        row['atom_indices'] = atom_indices
        records.append(row)
    seen = {}
    for record in records:
        for index in record['atom_indices']:
            if index in seen:
                raise ValueError(f'atom {index} belongs to more than one fragment')
            seen[index] = record['name']
    pairs = section.get('pairs', ())
    if not isinstance(pairs, Sequence) or isinstance(pairs, (str, bytes)):
        raise ValueError('fragment pairs must be a sequence')
    names = {record['name'] for record in records}
    normalized_pairs = []
    for pair in pairs:
        if not isinstance(pair, Mapping):
            raise ValueError('each fragment pair must be a mapping')
        first, second = pair.get('fragment_a'), pair.get('fragment_b')
        if first not in names or second not in names:
            continue
        value = pair.get('wbo', pair.get('value'))
        normalized_pairs.append({
            'fragment_a': str(first),
            'fragment_b': str(second),
            'value': _finite_number(value),
            'wbo': _finite_number(pair.get('wbo')),
        })
    return records, normalized_pairs, section


def _index_rows(rows, key):
    if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
        return {}
    indexed = {}
    for row in rows:
        if isinstance(row, Mapping) and key in row and not isinstance(row[key], bool):
            indexed[int(row[key])] = row
    return indexed


def _resolve_atom_values(report, field, natoms):
    if not field or not isinstance(report, Mapping):
        return {}
    # 1. Explicit frontier population aliases / dotted paths
    if field in {'homo_population', 'frontier.homo.atom_pop', 'lumo_population', 'frontier.lumo.atom_pop'}:
        group_name = 'homo' if 'homo' in field else 'lumo'
        frontier = report.get('frontier')
        if isinstance(frontier, Mapping):
            group = frontier.get(group_name)
            if isinstance(group, Mapping) and isinstance(group.get('atom_pop'), Sequence):
                return {
                    idx: _finite_number(val)
                    for idx, val in enumerate(group['atom_pop'])
                    if idx < natoms
                }
        # Fallback to top-level atoms rows if caller supplied custom column
        rows = _index_rows(report.get('atoms'), 'atom_index')
        return {
            idx: _finite_number(row.get(field))
            for idx, row in rows.items()
            if idx < natoms and field in row
        }

    section_name, leaf = field.split('.', 1) if '.' in field else (None, field)
    candidate_maps = []
    if section_name in (None, 'atoms'):
        candidate_maps.append(_index_rows(report.get('atoms'), 'atom_index'))
    if section_name in (None, 'shell_populations'):
        shells = report.get('shell_populations')
        if isinstance(shells, Mapping):
            candidate_maps.append(_index_rows(shells.get('atoms'), 'atom_index'))
    if section_name in (None, 'pi'):
        pi = report.get('pi')
        if isinstance(pi, Mapping):
            candidate_maps.append(_index_rows(pi.get('atoms'), 'atom_index'))

    for row_map in candidate_maps:
        if any(leaf in row for row in row_map.values()):
            return {
                idx: _finite_number(row.get(leaf))
                for idx, row in row_map.items()
                if idx < natoms
            }
    return {}


def _resolve_bond_values(report, field, nbonds):
    if not field or not isinstance(report, Mapping):
        return {}
    if field in {'homo_mul_ovlp_e', 'frontier_bonds.homo.mul_ovlp_e', 'lumo_mul_ovlp_e', 'frontier_bonds.lumo.mul_ovlp_e'}:
        group_name = 'homo' if 'homo' in field else 'lumo'
        fb = report.get('frontier_bonds')
        rows = fb.get('bonds') if isinstance(fb, Mapping) else None
        if isinstance(rows, Sequence) and not isinstance(rows, (str, bytes)):
            out = {}
            for row in rows:
                if (isinstance(row, Mapping)
                        and row.get('frontier_group') == group_name
                        and 'bond_index' in row):
                    idx = int(row['bond_index'])
                    if idx < nbonds:
                        out[idx] = _finite_number(row.get('mul_ovlp_e'))
            return out
        return {}

    section_name, leaf = field.split('.', 1) if '.' in field else (None, field)
    candidate_maps = []
    if section_name in (None, 'bonds'):
        candidate_maps.append(_index_rows(report.get('bonds'), 'bond_index'))
    if section_name in (None, 'pi'):
        pi = report.get('pi')
        if isinstance(pi, Mapping):
            candidate_maps.append(_index_rows(pi.get('bonds'), 'bond_index'))

    for row_map in candidate_maps:
        if any(leaf in row for row in row_map.values()):
            return {
                idx: _finite_number(row.get(leaf))
                for idx, row in row_map.items()
                if idx < nbonds
            }
    return {}


def _default_atom_field(report, natoms):
    for candidate in ('q_e', 'cm5_q_e', 'f_plus', 's_loc_1ev', 'homo_population', 'wbo_sum'):
        values = _resolve_atom_values(report, candidate, natoms)
        if any(val is not None for val in values.values()):
            return candidate
    return None


def BuildFragmentVisualData(mol, report, *, atomField=None, bondField=None):
    """Build a JSON-compatible mapping for molecular, bond and fragment views."""
    if not hasattr(mol, 'GetNumAtoms'):
        raise TypeError('mol must be an RDKit molecule')
    report = _normalize_report(mol, report)
    _validate_provenance(mol, report)
    records, normalized_pairs, section = _extract_fragments(report, required=False)
    natoms = int(mol.GetNumAtoms())
    nbonds = int(mol.GetNumBonds()) if hasattr(mol, 'GetNumBonds') else 0
    assignments = {}
    for record in records:
        for index in record['atom_indices']:
            if index >= natoms:
                raise ValueError(f'fragment atom index {index} is outside the molecule')
            assignments[index] = record['name']
    if atomField is None and bondField is None:
        atomField = _default_atom_field(report, natoms)
    atom_values = _resolve_atom_values(report, atomField, natoms) if atomField else {}
    atoms = [
        {
            'atom_index': index,
            'fragment': assignments.get(index),
            'value': atom_values.get(index),
        }
        for index in range(natoms)
    ]
    bond_values = _resolve_bond_values(report, bondField, nbonds) if bondField else {}
    bonds = [
        {
            'bond_index': bond.GetIdx(),
            'begin_atom_index': bond.GetBeginAtomIdx(),
            'end_atom_index': bond.GetEndAtomIdx(),
            'value': bond_values.get(bond.GetIdx()),
        }
        for bond in (mol.GetBonds() if hasattr(mol, 'GetBonds') else ())
    ]
    return {
        'atom_field': atomField,
        'atom_unit': _resolve_unit(report, atomField),
        'bond_field': bondField,
        'bond_unit': _resolve_unit(report, bondField),
        'atoms': atoms,
        'bonds': bonds,
        'fragments': records,
        'pairs': normalized_pairs,
        'unassigned_atom_indices': sorted(set(range(natoms)) - set(assignments)),
        'units': dict(report.get('units', {})) if isinstance(report, Mapping) else {},
        'source': section.get('source', 'EleD descriptor report') if isinstance(section, Mapping) else 'EleD descriptor report',
    }


def _strip_xml_header(svg):
    if svg.lstrip().startswith('<?xml') and '?>' in svg:
        return svg.split('?>', 1)[1].lstrip('\n')
    return svg


def _prepare_draw_mol(mol, *, includeHydrogens):
    from rdkit import Chem
    from rdkit.Chem import rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    work = Chem.Mol(mol)
    for atom in work.GetAtoms():
        atom.SetIntProp('_eled_orig_idx', atom.GetIdx())
    for bond in work.GetBonds():
        bond.SetIntProp('_eled_orig_bond_idx', bond.GetIdx())

    if includeHydrogens:
        draw_mol = work
    else:
        try:
            draw_mol = Chem.RemoveHs(work, sanitize=True)
        except Exception:
            draw_mol = Chem.RemoveHs(work, sanitize=False)
            draw_mol.UpdatePropertyCache(strict=False)

    original_to_draw = {}
    for atom in draw_mol.GetAtoms():
        orig_idx = atom.GetIntProp('_eled_orig_idx')
        original_to_draw[orig_idx] = atom.GetIdx()
        atom.ClearProp('_eled_orig_idx')

    orig_bond_to_draw = {}
    for bond in draw_mol.GetBonds():
        orig_bidx = bond.GetIntProp('_eled_orig_bond_idx')
        orig_bond_to_draw[orig_bidx] = bond.GetIdx()
        bond.ClearProp('_eled_orig_bond_idx')

    if draw_mol.GetNumConformers() == 0 or draw_mol.GetConformer().Is3D():
        rdDepictor.Compute2DCoords(draw_mol)
    try:
        draw_mol = rdMolDraw2D.PrepareMolForDrawing(draw_mol, kekulize=True, addChiralHs=False)
    except Exception:
        draw_mol = rdMolDraw2D.PrepareMolForDrawing(draw_mol, kekulize=False, addChiralHs=False)
    return draw_mol, original_to_draw, orig_bond_to_draw


def _is_diverging_field(field, values):
    if field in _SIGNED_FIELDS or (field and field.rsplit('.', 1)[-1] in _SIGNED_FIELDS):
        return True
    return any(v < 0.0 for v in values)


def _format_value(value, field, valueFormat):
    if valueFormat is not None:
        formatted = format(value, valueFormat)
    elif _is_diverging_field(field, (value,)):
        formatted = f'{value:+.2f}'
    else:
        formatted = f'{value:.2f}'
    if formatted.startswith('-') and set(formatted[1:]) <= {'0', '.'}:
        return ('+' + formatted[1:]) if (valueFormat is None and _is_diverging_field(field, (value,))) else formatted[1:]
    return formatted


def _append_colorbar(svg, *, width, height, label, unit, vmin, vmax, diverging, fmt):
    bar_h = 42
    total_h = height + bar_h
    grad_id = 'eled-colorbar-grad'
    if diverging:
        bound = max(abs(vmin), abs(vmax)) or 1.0
        left_val, mid_val, right_val = -bound, 0.0, bound
        stops = (
            f'<stop offset="0%" stop-color="{_svg_color(_DIV_NEG)}"/>'
            f'<stop offset="50%" stop-color="{_svg_color(_DIV_MID)}"/>'
            f'<stop offset="100%" stop-color="{_svg_color(_DIV_POS)}"/>'
        )
    else:
        left_val, right_val = vmin, vmax
        mid_val = (vmin + vmax) / 2.0
        stops = (
            f'<stop offset="0%" stop-color="{_svg_color(_SEQ_LOW)}"/>'
            f'<stop offset="100%" stop-color="{_svg_color(_SEQ_HIGH)}"/>'
        )
    bar_w = max(130, min(200, int(width * 0.44)))
    bar_x = width - bar_w - 16
    bar_y = height + 10
    caption = escape(f'{label} ({unit})' if unit and unit != 'dimensionless' else label)
    left_txt = escape(format(left_val, fmt))
    mid_txt = escape(format(mid_val, fmt))
    right_txt = escape(format(right_val, fmt))
    legend = (
        f'<defs><linearGradient id="{grad_id}" x1="0%" y1="0%" x2="100%" y2="0%">{stops}</linearGradient></defs>'
        f'<rect x="0" y="{height}" width="{width}" height="{bar_h}" fill="#ffffff"/>'
        f'<line x1="16" y1="{height}" x2="{width - 16}" y2="{height}" stroke="#e2e8f0" stroke-width="1"/>'
        f'<text x="16" y="{bar_y + 13}" font-family="system-ui, -apple-system, sans-serif" font-size="11.5" font-weight="500" fill="#334155">{caption}</text>'
        f'<rect x="{bar_x:.1f}" y="{bar_y}" width="{bar_w}" height="8" rx="4" fill="url(#{grad_id})" stroke="#cbd5e1" stroke-width="0.6"/>'
        f'<text x="{bar_x:.1f}" y="{bar_y + 22}" font-family="system-ui, -apple-system, sans-serif" font-size="10.5" fill="#64748b">{left_txt}</text>'
        f'<text x="{bar_x + bar_w / 2:.1f}" y="{bar_y + 22}" text-anchor="middle" font-family="system-ui, -apple-system, sans-serif" font-size="10.5" fill="#64748b">{mid_txt}</text>'
        f'<text x="{bar_x + bar_w:.1f}" y="{bar_y + 22}" text-anchor="end" font-family="system-ui, -apple-system, sans-serif" font-size="10.5" fill="#64748b">{right_txt}</text>'
    )
    svg = svg.replace(f"height='{height}px'", f"height='{total_h}px'", 1)
    svg = svg.replace(f"viewBox='0 0 {width} {height}'", f"viewBox='0 0 {width} {total_h}'", 1)
    idx = svg.rfind('</svg>')
    return f'{svg[:idx]}{legend}</svg>' if idx != -1 else svg


def DrawDescriptorMap(
    mol,
    report,
    *,
    atomField=None,
    bondField=None,
    colorBy='auto',
    labelStyle='note',
    valueFormat=None,
    width=520,
    height=240,
    addAtomIndices=False,
    includeHydrogens=False,
    showColorbar=None,
    drawOptions=None,
):
    """Return an SVG RDKit drawing with fragments, atom values and/or bond values."""
    from rdkit.Chem.Draw import rdMolDraw2D

    if colorBy not in {'auto', 'fragment', 'value', 'none'}:
        raise ValueError("colorBy must be 'auto', 'fragment', 'value', or 'none'")
    if labelStyle not in {'note', 'replace', 'none'}:
        raise ValueError("labelStyle must be 'note', 'replace', or 'none'")
    if int(width) <= 0 or int(height) <= 0:
        raise ValueError('width and height must be positive integers')

    width, height = int(width), int(height)
    data = BuildFragmentVisualData(mol, report, atomField=atomField, bondField=bondField)
    records = data['fragments']
    assignments = {row['atom_index']: row['fragment'] for row in data['atoms']}

    draw_mol, original_to_draw, orig_bond_to_draw = _prepare_draw_mol(
        mol, includeHydrogens=bool(includeHydrogens)
    )

    atom_finite = [
        row['value']
        for row in data['atoms']
        if row['value'] is not None and row['atom_index'] in original_to_draw
    ]
    bond_finite = [
        row['value']
        for row in data['bonds']
        if row['value'] is not None and row['bond_index'] in orig_bond_to_draw
    ]

    if colorBy == 'auto':
        if records and bondField is None:
            resolved_color = 'fragment'
        elif atom_finite or bond_finite:
            resolved_color = 'value'
        elif records:
            resolved_color = 'fragment'
        else:
            resolved_color = 'none'
    else:
        resolved_color = colorBy

    highlight_atoms = []
    highlight_atom_colors = {}
    highlight_bonds = []
    highlight_bond_colors = {}
    colorbar_spec = None

    if resolved_color == 'fragment':
        colors = {record['name']: _pastel_rgb(_COLORS[index % len(_COLORS)]) for index, record in enumerate(records)}
        for orig_idx, name in assignments.items():
            draw_idx = original_to_draw.get(orig_idx)
            if name is not None and draw_idx is not None:
                highlight_atoms.append(draw_idx)
                highlight_atom_colors[draw_idx] = colors[name]
    elif resolved_color == 'value':
        if atom_finite:
            div = _is_diverging_field(data['atom_field'], atom_finite)
            vmin, vmax = min(atom_finite), max(atom_finite)
            for row in data['atoms']:
                draw_idx = original_to_draw.get(row['atom_index'])
                if draw_idx is not None and row['value'] is not None:
                    highlight_atoms.append(draw_idx)
                    highlight_atom_colors[draw_idx] = _value_rgb(row['value'], vmin, vmax, diverging=div)
            colorbar_spec = (
                _field_label(data['atom_field']),
                data['atom_unit'],
                vmin,
                vmax,
                div,
                valueFormat or ('+.2f' if div else '.2f'),
            )
        if bond_finite:
            bdiv = _is_diverging_field(data['bond_field'], bond_finite)
            bmin, bmax = min(bond_finite), max(bond_finite)
            for row in data['bonds']:
                draw_bidx = orig_bond_to_draw.get(row['bond_index'])
                if draw_bidx is not None and row['value'] is not None:
                    highlight_bonds.append(draw_bidx)
                    highlight_bond_colors[draw_bidx] = _value_rgb(row['value'], bmin, bmax, diverging=bdiv)
            if colorbar_spec is None:
                colorbar_spec = (
                    _field_label(data['bond_field']),
                    data['bond_unit'],
                    bmin,
                    bmax,
                    bdiv,
                    valueFormat or ('+.2f' if bdiv else '.2f'),
                )

    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    options = drawer.drawOptions()
    font_size = max(11, min(16, int(min(width, height) / 17)))
    options.fixedFontSize = font_size
    options.minFontSize = max(8, font_size - 3)
    options.maxFontSize = font_size
    options.bondLineWidth = 1.8
    options.padding = 0.14
    options.additionalAtomLabelPadding = 0.08
    options.atomHighlightsAreCircles = True
    options.highlightBondWidthMultiplier = 7
    options.addAtomIndices = False
    options.annotationFontScale = 0.78
    options.setAnnotationColour((0.12, 0.16, 0.23))
    if resolved_color == 'value':
        options.useBWAtomPalette()
    highlight_radii = {idx: 0.28 for idx in highlight_atoms}

    if labelStyle != 'none' and data['atom_field']:
        for row in data['atoms']:
            orig_idx = int(row['atom_index'])
            draw_idx = original_to_draw.get(orig_idx)
            if draw_idx is not None and row['value'] is not None:
                formatted = _format_value(row['value'], data['atom_field'], valueFormat)
                text = f'{orig_idx}:{formatted}' if addAtomIndices else formatted
                if labelStyle == 'note':
                    draw_mol.GetAtomWithIdx(draw_idx).SetProp('atomNote', text)
                else:
                    options.atomLabels[draw_idx] = text
    elif addAtomIndices:
        for orig_idx, draw_idx in original_to_draw.items():
            draw_mol.GetAtomWithIdx(draw_idx).SetProp('atomNote', str(orig_idx))

    if labelStyle != 'none' and data['bond_field']:
        for row in data['bonds']:
            draw_bidx = orig_bond_to_draw.get(int(row['bond_index']))
            if draw_bidx is not None and row['value'] is not None:
                formatted = _format_value(row['value'], data['bond_field'], valueFormat)
                draw_mol.GetBondWithIdx(draw_bidx).SetProp('bondNote', formatted)

    if callable(drawOptions):
        drawOptions(options)
    elif drawOptions is not None:
        raise TypeError('drawOptions must be a callable accepting rdMolDraw2D.MolDrawOptions')

    drawer.DrawMolecule(
        draw_mol,
        highlightAtoms=highlight_atoms,
        highlightAtomColors=highlight_atom_colors,
        highlightBonds=highlight_bonds,
        highlightBondColors=highlight_bond_colors,
        highlightAtomRadii=highlight_radii,
    )
    drawer.FinishDrawing()
    svg = _strip_xml_header(drawer.GetDrawingText())
    want_colorbar = (resolved_color == 'value') if showColorbar is None else bool(showColorbar)
    if want_colorbar and colorbar_spec is not None:
        label, unit, vmin, vmax, div, fmt = colorbar_spec
        svg = _append_colorbar(
            svg,
            width=width,
            height=height,
            label=label,
            unit=unit,
            vmin=vmin,
            vmax=vmax,
            diverging=div,
            fmt=fmt,
        )
    return svg


def _truncate_label(text, max_chars):
    if len(text) <= max_chars:
        return text
    return f'{text[:max(1, max_chars - 1)]}…'


def PlotFragmentSummary(report, *, field='q_e', width=620, rowHeight=38):
    """Return an SVG signed bar chart for one fragment field."""
    width, rowHeight = int(width), int(rowHeight)
    if width <= 0 or rowHeight <= 0:
        raise ValueError('width and rowHeight must be positive integers')
    records, _, _ = _extract_fragments(report, required=True)
    values = [_finite_number(row.get(field)) for row in records]
    finite = [value for value in values if value is not None]
    height = max(130, 82 + len(records) * rowHeight)
    unit = _resolve_unit(report, field)
    title_text = f'{_field_label(field)} ({unit})' if unit != 'dimensionless' else _field_label(field)
    title = escape(title_text)
    if not finite:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
            f'role="img" aria-label="No {title} values available">{_SVG_STYLE}'
            f'<rect class="eled-viz-bg" width="{width}" height="{height}"/>'
            f'<title>No {title} values available</title>'
            f'<text class="eled-viz-text" x="20" y="48" font-size="14">No values available</text></svg>'
        )
    maximum = max(abs(value) for value in finite) or 1.0
    longest = max((len(str(record['name'])) for record in records), default=8)
    right_margin = 76
    left_margin = min(max(96, 20 + longest * 8), max(96, int(width * 0.36)))
    chart_width = max(60, width - left_margin - right_margin)
    center = left_margin + chart_width / 2.0
    max_label_chars = max(6, int((left_margin - 20) / 8))

    baseline_top = 50
    baseline_bottom = height - 24
    body = [
        f'<line class="eled-viz-muted" x1="{center:.1f}" y1="{baseline_top}" '
        f'x2="{center:.1f}" y2="{baseline_bottom}" opacity=".55"/>',
        f'<text class="eled-viz-sub" x="{center:.1f}" y="{height - 8}" '
        f'text-anchor="middle" font-size="11">0</text>',
    ]
    for index, (record, value) in enumerate(zip(records, values)):
        y = 70 + index * rowHeight
        raw_name = str(record['name'])
        display_name = escape(_truncate_label(raw_name, max_label_chars))
        body.append(
            f'<text class="eled-viz-text" x="12" y="{y + 5}" font-size="13">'
            f'<title>{escape(raw_name)}</title>{display_name}</text>'
        )
        if value is not None:
            span = abs(value) / maximum * (chart_width / 2.0 - 4.0)
            x = center if value >= 0 else center - span
            fill = _SVG_SERIES[0] if value >= 0 else _SVG_SERIES[1]
            body.append(
                f'<rect x="{x:.1f}" y="{y - 8}" width="{max(1.0, span):.1f}" '
                f'height="16" fill="{fill}" rx="3"/>'
            )
            body.append(
                f'<text class="eled-viz-text" x="{width - 12}" y="{y + 5}" '
                f'text-anchor="end" font-size="13">{value:+.3f}</text>'
            )
        else:
            body.append(
                f'<text class="eled-viz-sub" x="{width - 12}" y="{y + 5}" '
                f'text-anchor="end" font-size="12" font-style="italic">missing</text>'
            )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'role="img" aria-label="Fragment {title} summary">{_SVG_STYLE}'
        f'<rect class="eled-viz-bg" width="{width}" height="{height}"/>'
        f'<title>Fragment {title} summary</title>'
        f'<text class="eled-viz-text" x="12" y="26" font-size="15" font-weight="600">{title}</text>'
        f'{"".join(body)}</svg>'
    )


def PlotFragmentInteractions(report, *, width=620, height=280):
    """Return an SVG network of fragment-pair Wiberg couplings."""
    width, height = int(width), int(height)
    if width <= 0 or height <= 0:
        raise ValueError('width and height must be positive integers')
    records, pairs, _ = _extract_fragments(report, required=True)
    names = [record['name'] for record in records]
    if not names:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" role="img">'
            f'{_SVG_STYLE}<rect class="eled-viz-bg" width="{width}" height="{height}"/>'
            f'<text class="eled-viz-text" x="20" y="40" font-size="14">No fragments available</text></svg>'
        )
    node_dims = {}
    for name in names:
        display = _truncate_label(name, 16)
        w = max(76, min(148, 24 + len(display) * 8))
        node_dims[name] = (w, 38, display)

    max_half_w = max(dims[0] / 2.0 for dims in node_dims.values())
    positions = {}
    if len(names) == 1:
        positions[names[0]] = (width / 2.0, height / 2.0)
    else:
        rx = max(40.0, width / 2.0 - max_half_w - 20.0)
        ry = max(40.0, height / 2.0 - 44.0)
        for index, name in enumerate(names):
            angle = 2.0 * math.pi * index / len(names)
            positions[name] = (
                width / 2.0 + rx * math.cos(angle),
                height / 2.0 + ry * math.sin(angle),
            )

    finite = [pair['value'] for pair in pairs if pair['value'] is not None]
    maximum = max((abs(value) for value in finite), default=1.0) or 1.0
    name_order = {name: idx for idx, name in enumerate(names)}
    body = []
    for pair in pairs:
        if pair['value'] is None:
            continue
        x1, y1 = positions[pair['fragment_a']]
        x2, y2 = positions[pair['fragment_b']]
        stroke_width = 1.5 + 9.5 * abs(pair['value']) / maximum
        body.append(
            f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" '
            f'stroke="{_SVG_SERIES[2]}" stroke-width="{stroke_width:.1f}" stroke-linecap="round" opacity=".78"/>'
        )
        # Offset non-adjacent chords away from the center to avoid label collisions.
        n = len(names)
        delta = abs(name_order[pair['fragment_a']] - name_order[pair['fragment_b']])
        t = 0.36 if (n >= 4 and min(delta, n - delta) > 1) else 0.5
        lx = x1 + (x2 - x1) * t
        ly = y1 + (y2 - y1) * t
        dx, dy = x2 - x1, y2 - y1
        length = math.hypot(dx, dy) or 1.0
        nx, ny = -dy / length, dx / length
        lx += nx * 12.0
        ly += ny * 12.0
        body.append(
            f'<text class="eled-viz-text" x="{lx:.1f}" y="{ly + 4:.1f}" '
            f'text-anchor="middle" font-size="13">{pair["value"]:.2f}</text>'
        )
    for index, name in enumerate(names):
        x, y = positions[name]
        w, h, display = node_dims[name]
        fill = _svg_color(_COLORS[index % len(_COLORS)])
        body.append(
            f'<rect x="{x - w / 2:.1f}" y="{y - h / 2:.1f}" width="{w:.1f}" height="{h:.1f}" '
            f'rx="19" fill="{fill}" opacity=".92"><title>{escape(name)}</title></rect>'
        )
        body.append(
            f'<text x="{x:.1f}" y="{y + 5:.1f}" text-anchor="middle" '
            f'fill="#ffffff" font-family="system-ui, sans-serif" font-size="13">{escape(display)}</text>'
        )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" '
        f'role="img" aria-label="Fragment interaction network">{_SVG_STYLE}'
        f'<rect class="eled-viz-bg" width="{width}" height="{height}"/>'
        f'<title>Fragment interaction network</title>'
        f'{"".join(body)}</svg>'
    )


def _select_summary_field(records, summaryField, atomField):
    if summaryField is not None:
        return summaryField
    for candidate in (atomField, 'q_e', 'cm5_q_e', 'f_plus', 's_loc_1ev', 'dual'):
        if candidate and any(_finite_number(row.get(candidate)) is not None for row in records):
            return candidate
    return atomField or 'q_e'


def RenderFragmentAnalysis(
    mol,
    report,
    *,
    atomField=None,
    bondField=None,
    summaryField=None,
    colorBy='auto',
    labelStyle='note',
    includeHydrogens=False,
    addAtomIndices=False,
    width=620,
    height=240,
):
    """Return a self-contained HTML fragment combining EleD molecular and fragment views."""
    report = _normalize_report(mol, report)
    data = BuildFragmentVisualData(mol, report, atomField=atomField, bondField=bondField)
    molecule = DrawDescriptorMap(
        mol,
        report,
        atomField=atomField,
        bondField=bondField,
        colorBy=colorBy,
        labelStyle=labelStyle,
        width=width,
        height=height,
        includeHydrogens=includeHydrogens,
        addAtomIndices=addAtomIndices,
    )
    active_label = _field_label(data['atom_field'] or data['bond_field'] or 'descriptor')
    sections = [
        f'<section><h2>RDKit molecule with EleD {escape(active_label)}</h2>{molecule}</section>'
    ]
    if data['fragments']:
        chosen_summary = _select_summary_field(data['fragments'], summaryField, data['atom_field'])
        summary = PlotFragmentSummary(report, field=chosen_summary, width=width)
        sections.append(f'<section><h2>Fragment summary</h2>{summary}</section>')
        interactions = PlotFragmentInteractions(report, width=width)
        sections.append(
            f'<section><h2>Fragment interaction network</h2>{interactions}'
            f'<div class="note">Edges are fragment-pair Wiberg sums.</div></section>'
        )
    return f'''<div class="eled-fragment-analysis">
<style>
.eled-fragment-analysis {{ font-family: system-ui, -apple-system, sans-serif; color: #172033; background: transparent; display: grid; gap: 12px; max-width: {int(width)}px; }}
.eled-fragment-analysis section {{ background: #ffffff; border: 1px solid #d8dde5; border-radius: 8px; padding: 12px; }}
.eled-fragment-analysis h2 {{ margin: 0 0 8px; font-size: 16px; }}
.eled-fragment-analysis svg {{ max-width: 100%; height: auto; }}
.eled-fragment-analysis .note {{ color: #59636f; font-size: 12px; margin-top: 6px; }}
@media (prefers-color-scheme: dark) {{
  .eled-fragment-analysis {{ color: #f8fafc; }}
  .eled-fragment-analysis section {{ background: #111827; border-color: #334155; }}
  .eled-fragment-analysis .note {{ color: #94a3b8; }}
}}
</style>
{"".join(sections)}
</div>'''


__all__ = [
    'BuildFragmentVisualData',
    'DrawDescriptorMap',
    'PlotFragmentSummary',
    'PlotFragmentInteractions',
    'RenderFragmentAnalysis',
]
