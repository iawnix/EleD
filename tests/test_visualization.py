from rdkit import Chem
from rdkit.Chem import AllChem
import pytest

from eled import (
    BuildFragmentVisualData,
    DrawDescriptorMap,
    PlotFragmentInteractions,
    PlotFragmentSummary,
    RenderFragmentAnalysis,
)
from eled.input import _geometry_hash, _ordered_structure_hash, _xyz_from_mol
from eled.result import XTBResults


def _report():
    mol = Chem.MolFromSmiles('CC(=O)Oc1ccccc1')
    return mol, {
        'atoms': [
            {'atom_index': i, 'q_e': (i - 4) / 100., 'cm5_q_e': (i - 4) / 90.,
             'f_plus': (i + 1) / 55.}
            for i in range(mol.GetNumAtoms())
        ],
        'bonds': [
            {'bond_index': b.GetIdx(), 'begin_atom_index': b.GetBeginAtomIdx(),
             'end_atom_index': b.GetEndAtomIdx(), 'wbo': 1.0 + 0.05 * b.GetIdx()}
            for b in mol.GetBonds()
        ],
        'fragments': {
            'fragments': {
                'acetyl': {'atom_indices': [0, 1, 2], 'q_e': -0.10, 'f_plus': 0.31},
                'phenoxy': {'atom_indices': [3, 4, 5, 6, 7, 8, 9], 'q_e': 0.10, 'f_plus': 0.69},
            },
            'pairs': [{'fragment_a': 'acetyl', 'fragment_b': 'phenoxy', 'wbo': 1.2}],
        },
    }


def _water_result():
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    xyz = _xyz_from_mol(mol, -1)
    result = XTBResults(
        method='GFN2-xTB',
        totalEnergyHartree=-5.07,
        orbitalEnergiesEV=(-15.0, -10.0, 2.0),
        orbitalOccupancies=(2.0, 2.0, 0.0),
        atomicCharges=(-0.60, 0.30, 0.30),
        cm5Charges=(-0.55, 0.275, 0.275),
        dipoleAtomicUnits=(0.0, 0.0, 0.72),
        wibergBondOrders=(
            (0.0, 0.92, 0.92),
            (0.92, 0.0, 0.01),
            (0.92, 0.01, 0.0),
        ),
        numElectrons=8,
        numUnpairedElectrons=0,
        chargeModel='xtb_scc_mulliken',
        metadata={
            'atom_count': mol.GetNumAtoms(),
            'charge': 0,
            'conformer_id': 0,
            'ordered_structure_sha256': _ordered_structure_hash(mol),
            'geometry_sha256': _geometry_hash(xyz),
        },
    )
    return mol, result


def test_visual_data_preserves_rdkit_indices_and_pair_values():
    mol, report = _report()
    data = BuildFragmentVisualData(mol, report, atomField='f_plus')
    assert data['atom_field'] == 'f_plus'
    assert data['unassigned_atom_indices'] == []
    assert data['atoms'][0] == {'atom_index': 0, 'fragment': 'acetyl', 'value': pytest.approx(1 / 55)}
    assert data['pairs'][0]['value'] == pytest.approx(1.2)


def test_rdkit_map_and_svg_views_are_self_contained():
    mol, report = _report()
    molecule = DrawDescriptorMap(mol, report, atomField='q_e')
    summary = PlotFragmentSummary(report, field='q_e')
    network = PlotFragmentInteractions(report)
    html = RenderFragmentAnalysis(mol, report)
    assert '<svg' in molecule and '<?xml' not in molecule
    assert 'acetyl' in summary and 'phenoxy' in summary
    assert '1.20' in network and 'Fragment interaction network' in network
    assert 'eled-viz-bg' in summary and 'currentColor' not in summary
    assert 'eled-viz-bg' in network and 'currentColor' not in network
    assert '<section>' in html and '<svg' in html


def test_explicit_hydrogens_are_hidden_from_map_by_default():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    report = {
        'atoms': [{'atom_index': index, 'q_e': 0.01 * index}
                  for index in range(mol.GetNumAtoms())],
        'fragments': {
            'fragments': {'all': {'atom_indices': list(range(mol.GetNumAtoms()))}},
            'pairs': [],
        },
    }
    compact = DrawDescriptorMap(mol, report, atomField='q_e')
    complete = DrawDescriptorMap(mol, report, atomField='q_e',
                                 includeHydrogens=True)
    assert "class='atom-8'" not in compact
    assert "class='atom-8'" in complete
    # Removing explicit hydrogens must restore implicit H on heteroatoms (OH, not bare O)
    assert compact.count("class='atom-2'") >= 3


def test_visual_data_rejects_overlapping_fragments():
    mol, report = _report()
    report['fragments']['fragments']['phenoxy']['atom_indices'] = [2, 3, 4, 5, 6, 7, 8, 9]
    with pytest.raises(ValueError, match='more than one fragment'):
        BuildFragmentVisualData(mol, report)


def test_draw_descriptor_map_without_fragments_uses_value_colormap_and_colorbar():
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    report = {
        'atoms': [{'atom_index': i, 'q_e': -0.44 if i == 2 else 0.08 * i}
                  for i in range(mol.GetNumAtoms())],
        'fragments': None,
    }
    svg = DrawDescriptorMap(mol, report, atomField='q_e')
    assert '<svg' in svg
    assert 'eled-colorbar-grad' in svg
    assert 'charge (e)' in svg
    # Input 3D conformer must remain untouched
    assert mol.GetConformer().Is3D()
    # Colorbar can be disabled explicitly
    no_bar = DrawDescriptorMap(mol, report, atomField='q_e', showColorbar=False)
    assert 'eled-colorbar-grad' not in no_bar


def test_remove_hs_preserves_indices_with_isotopic_hydrogen():
    mol = Chem.AddHs(Chem.MolFromSmiles('[2H]C(O)C'))
    report = {
        'atoms': [{'atom_index': i, 'q_e': -0.50 if i == 2 else 0.05 * i}
                  for i in range(mol.GetNumAtoms())],
        'fragments': None,
    }
    # Deuterium is retained by RDKit RemoveHs; mapping must not fail or misalign.
    svg = DrawDescriptorMap(mol, report, atomField='q_e', addAtomIndices=True)
    assert '<svg' in svg
    assert "class='atom-3'" in svg


def test_bond_field_and_nested_report_sections_resolution():
    mol, report = _report()
    report['frontier'] = {
        'homo': {'atom_pop': [0.1] * mol.GetNumAtoms()},
        'lumo': {'atom_pop': [0.05] * mol.GetNumAtoms()},
    }
    report['shell_populations'] = {
        'atoms': [{'atom_index': i, 'mul_p_e': 1.2 + 0.1 * i} for i in range(mol.GetNumAtoms())],
        'units': {'mul_p_e': 'e'},
    }
    report['pi'] = {
        'atoms': [{'atom_index': i, 'pi_q_e': -0.02 if i >= 4 else None} for i in range(mol.GetNumAtoms())],
        'bonds': [{'bond_index': b.GetIdx(), 'pi_bo': 0.35 if b.GetIdx() >= 4 else None} for b in mol.GetBonds()],
        'units': {'pi_q_e': 'e', 'pi_bo': 'dimensionless'},
    }
    report['frontier_bonds'] = {
        'bonds': [
            {'frontier_group': 'homo', 'bond_index': b.GetIdx(), 'mul_ovlp_e': -0.08 if b.GetIdx() % 2 else 0.12}
            for b in mol.GetBonds()
        ],
        'units': {'mul_ovlp_e': 'e'},
    }

    homo_data = BuildFragmentVisualData(mol, report, atomField='homo_population')
    assert homo_data['atoms'][0]['value'] == pytest.approx(0.1)
    assert homo_data['atom_unit'] == 'dimensionless'

    shell_data = BuildFragmentVisualData(mol, report, atomField='mul_p_e')
    assert shell_data['atoms'][1]['value'] == pytest.approx(1.3)
    assert shell_data['atom_unit'] == 'e'

    pi_data = BuildFragmentVisualData(mol, report, atomField='pi_q_e', bondField='pi_bo')
    assert pi_data['atoms'][0]['value'] is None
    assert pi_data['atoms'][4]['value'] == pytest.approx(-0.02)
    assert pi_data['bonds'][4]['value'] == pytest.approx(0.35)

    fb_data = BuildFragmentVisualData(mol, report, bondField='frontier_bonds.homo.mul_ovlp_e')
    assert fb_data['atom_field'] is None
    assert fb_data['bond_unit'] == 'e'
    assert fb_data['bonds'][0]['value'] == pytest.approx(0.12)
    assert fb_data['bonds'][1]['value'] == pytest.approx(-0.08)

    bond_svg = DrawDescriptorMap(mol, report, bondField='pi_bo', colorBy='value')
    assert '<svg' in bond_svg
    assert 'π bond order' in bond_svg


def test_draw_descriptor_map_accepts_xtb_results_and_validates_provenance():
    mol, result = _water_result()
    svg = DrawDescriptorMap(mol, result, atomField='cm5_q_e', bondField='wbo')
    assert '<svg' in svg
    assert 'CM5 charge (e)' in svg

    wrong_mol = Chem.AddHs(Chem.MolFromSmiles('N'))
    with pytest.raises(ValueError, match='atom count|structure or atom order'):
        DrawDescriptorMap(wrong_mol, result)

    # Same atom count (3 atoms), different structure
    other_three = Chem.MolFromSmiles('CCC')
    bad_report = {
        'atoms': [{'atom_index': i, 'q_e': 0.0} for i in range(3)],
        'provenance': {
            'atom_count': 3,
            'ordered_structure_sha256': _ordered_structure_hash(mol),
        },
    }
    with pytest.raises(ValueError, match='structure or atom order'):
        BuildFragmentVisualData(other_three, bad_report)


def test_label_styles_and_draw_options_callback():
    mol, report = _report()
    called = []
    svg_note = DrawDescriptorMap(
        mol, report, atomField='q_e', labelStyle='note',
        drawOptions=lambda opts: called.append(opts.annotationFontScale),
    )
    svg_replace = DrawDescriptorMap(mol, report, atomField='q_e', labelStyle='replace')
    svg_none = DrawDescriptorMap(mol, report, atomField='q_e', labelStyle='none')
    assert len(called) == 1
    assert '<svg' in svg_note and '<svg' in svg_replace and '<svg' in svg_none
    assert svg_note != svg_replace

    with pytest.raises(ValueError, match='colorBy'):
        DrawDescriptorMap(mol, report, colorBy='invalid')
    with pytest.raises(ValueError, match='labelStyle'):
        DrawDescriptorMap(mol, report, labelStyle='invalid')
    with pytest.raises(TypeError, match='drawOptions'):
        DrawDescriptorMap(mol, report, drawOptions='not-callable')


def test_plot_fragment_interactions_handles_zero_wbo_and_four_fragments():
    report = {
        'fragments': {
            'fragments': {
                'frag_alpha_long_name': {'atom_indices': [0]},
                'frag_beta': {'atom_indices': [1]},
                'frag_gamma': {'atom_indices': [2]},
                'frag_delta': {'atom_indices': [3]},
            },
            'pairs': [
                {'fragment_a': 'frag_alpha_long_name', 'fragment_b': 'frag_gamma', 'wbo': 0.0},
                {'fragment_a': 'frag_beta', 'fragment_b': 'frag_delta', 'wbo': 0.0},
            ],
        }
    }
    svg = PlotFragmentInteractions(report)
    assert '<svg' in svg
    assert '0.00' in svg
    assert '<title>frag_alpha_long_name</title>' in svg


def test_plot_fragment_summary_narrow_width_and_long_names():
    report = {
        'fragments': {
            'fragments': {
                'very_long_substituent_fragment_name': {'atom_indices': [0, 1], 'f_plus': 0.42},
                'missing_group': {'atom_indices': [2], 'f_plus': None},
            },
            'pairs': [],
        }
    }
    svg = PlotFragmentSummary(report, field='f_plus', width=220)
    assert '<svg' in svg
    assert 'width="-' not in svg
    assert 'Fukui f+' in svg
    assert 'missing' in svg
    assert '<title>very_long_substituent_fragment_name</title>' in svg


def test_render_fragment_analysis_links_summary_field_and_handles_no_fragments():
    mol, report = _report()
    html_fplus = RenderFragmentAnalysis(mol, report, atomField='f_plus', includeHydrogens=False)
    assert 'Fukui f+' in html_fplus
    assert 'prefers-color-scheme: dark' in html_fplus
    assert html_fplus.count('<section>') == 3

    no_frag = {'atoms': report['atoms'], 'fragments': None}
    html_single = RenderFragmentAnalysis(mol, no_frag, atomField='q_e')
    assert html_single.count('<section>') == 1
    assert 'eled-colorbar-grad' in html_single
