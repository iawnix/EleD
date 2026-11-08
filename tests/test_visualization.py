from rdkit import Chem
import pytest

from eled import (
    BuildFragmentVisualData, DrawDescriptorMap, PlotFragmentInteractions,
    PlotFragmentSummary, RenderFragmentAnalysis,
    build_fragment_visual_data, draw_descriptor_map,
)


def _report():
    mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1')
    return mol, {
        'atoms': [
            {'atom_index': i, 'partial_charge_e': (i - 4) / 100.,
             'f_plus': (i + 1) / 55.}
            for i in range(mol.GetNumAtoms())
        ],
        'fragments': {
            'fragments': {
                'acetyl': {'atom_indices': [0, 1, 2], 'charge_e': -0.10, 'f_plus': 0.31},
                'phenoxy': {'atom_indices': [3, 4, 5, 6, 7, 8, 9], 'charge_e': 0.10, 'f_plus': 0.69},
            },
            'pairs': [{'fragment_a': 'acetyl', 'fragment_b': 'phenoxy', 'wiberg_sum': 1.2}],
        },
    }


def test_visual_data_preserves_rdkit_indices_and_pair_values():
    mol, report = _report()
    data = BuildFragmentVisualData(mol, report, atom_field='f_plus')
    assert data['atom_field'] == 'f_plus'
    assert data['unassigned_atom_indices'] == []
    assert data['atoms'][0] == {'atom_index': 0, 'fragment': 'acetyl', 'value': pytest.approx(1 / 55)}
    assert data['pairs'][0]['value'] == pytest.approx(1.2)
    assert build_fragment_visual_data is BuildFragmentVisualData


def test_rdkit_map_and_svg_views_are_self_contained():
    mol, report = _report()
    molecule = DrawDescriptorMap(mol, report, atom_field='partial_charge_e')
    summary = PlotFragmentSummary(report, field='charge_e')
    network = PlotFragmentInteractions(report)
    html = RenderFragmentAnalysis(mol, report)
    assert '<svg' in molecule and '<?xml' not in molecule
    assert 'acetyl' in summary and 'phenoxy' in summary
    assert '1.20' in network and 'Fragment interaction network' in network
    assert '<section>' in html and '<svg' in html
    assert draw_descriptor_map is DrawDescriptorMap


def test_visual_data_rejects_overlapping_fragments():
    mol, report = _report()
    report['fragments']['fragments']['phenoxy']['atom_indices'] = [2, 3, 4, 5, 6, 7, 8, 9]
    with pytest.raises(ValueError, match='more than one fragment'):
        BuildFragmentVisualData(mol, report)
