"""Analytical geometry references and invariance of public local descriptors."""
import json

import numpy as np
import pytest
from rdkit import Chem

from eled import CalcBuriedVolume, CalcHOMA, CalcPlaneDescriptors, CalcRingPuckering, CalcSterimol


def molecule(smiles, coordinates):
    mol = Chem.MolFromSmiles(smiles)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for i, xyz in enumerate(coordinates):
        conf.SetAtomPosition(i, xyz)
    conf.SetId(7)
    mol.AddConformer(conf, assignId=False)
    return Chem.AddHs(mol, addCoords=True)


def ring(n=6, radius=1.388, heights=None):
    angle = np.arange(n)*2*np.pi/n
    coords = np.column_stack((radius*np.cos(angle), radius*np.sin(angle), np.zeros(n)))
    if heights is not None:
        coords[:, 2] = heights
    return molecule('C1' + 'C'*(n-1) + '1', coords)


def moved(mol):
    copy = Chem.Mol(mol)
    # Proper rotation with no axis aligned with its original direction.
    axis = np.array([1., 2., 3.]); axis /= np.linalg.norm(axis)
    cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    rotation = np.eye(3)*np.cos(.7) + (1-np.cos(.7))*np.outer(axis, axis) + np.sin(.7)*cross
    for i, p in enumerate(mol.GetConformer().GetPositions() @ rotation.T + [4., -2., 3.]):
        copy.GetConformer().SetAtomPosition(i, p)
    return copy


def test_homa_analytic_lengths_and_reordering():
    for length in (1.388, 1.5):
        mol = ring(radius=length)
        before = mol.ToBinary()
        for candidate in (mol, moved(mol)):
            for order in (list(range(6)), [2, 3, 4, 5, 0, 1], list(reversed(range(6)))):
                report = CalcHOMA(candidate, order)
                assert report['values']['homa'] == pytest.approx(1-257.7*(length-1.388)**2, abs=1e-12)
                assert report['provenance']['conformer_id'] == 7
                json.dumps(report, allow_nan=False)
        assert mol.ToBinary() == before
    order = list(reversed(range(mol.GetNumAtoms())))
    permuted = Chem.RenumberAtoms(mol, order)
    assert CalcHOMA(permuted, [order.index(i) for i in range(6)])['values'] == pytest.approx(CalcHOMA(mol, range(6))['values'])


def test_homa_heteroatom_parameters_and_ring_validation():
    mol = ring()
    mol.GetAtomWithIdx(0).SetAtomicNum(7)
    report = CalcHOMA(mol, range(6))
    assert report['values']['homa'] is None
    assert report['missing']['homa'] == 'unsupported_element_pair'
    params = {(6, 6): (1.388, 257.7), (6, 7): (1.388, 100.)}
    assert CalcHOMA(mol, range(6), parameters=params, parameter_source='analytic test')['values']['homa'] == pytest.approx(1.)
    with pytest.raises(ValueError, match='source'):
        CalcHOMA(mol, range(6), parameters=params)
    for indices in ([0, 2, 1, 3, 4, 5], [0, 1, 2], [0, 1, 1], [True, 1, 2]):
        with pytest.raises(ValueError):
            CalcHOMA(mol, indices)


def test_sterimol_axial_sphere_and_transverse_sphere_reference():
    mol = molecule('CCC', [(0, 0, 0), (1.5, 0, 0), (1.5, 2., 0)])
    radii = np.ones(mol.GetNumAtoms())
    for candidate in (mol, moved(mol)):
        axial = CalcSterimol(candidate, [1], axis=(0, 1), radii_angstrom=radii)
        assert axial['values'] == pytest.approx(dict(l_a=2.9, b1_a=1., b5_a=1.))
        full = CalcSterimol(candidate, [1, 2], axis=(0, 1), radii_angstrom=radii)
        assert full['values'] == pytest.approx(dict(l_a=2.9, b1_a=1., b5_a=3.))
        assert full['analysis_settings']['b1_error_bound_angstrom'] < .001
    for selection, axis in [([0, 1], [0, 1]), ([2], [0, 1]), ([2], [0, 2])]:
        with pytest.raises(ValueError):
            CalcSterimol(mol, selection, axis=axis)


def test_buried_volume_sphere_fraction_and_union():
    mol = molecule('[He].[He].[He]', [(0, 0, 0)]*3)
    report = CalcBuriedVolume(mol, [1], center=0, radii_angstrom=[1]*3,
                             sphere_radius_angstrom=2., spacing_angstrom=.05)
    assert report['values']['vbur_percent'] == pytest.approx(100/8, abs=.12)
    duplicate = CalcBuriedVolume(mol, [1, 2], center=0, radii_angstrom=[1]*3,
                                 sphere_radius_angstrom=2., spacing_angstrom=.05)
    assert duplicate['values'] == report['values']
    assert CalcBuriedVolume(mol, [1], center=0, radii_angstrom=[3]*3,
                            sphere_radius_angstrom=2.)['values']['vbur_percent'] == 100.
    with pytest.raises(ValueError, match='exclude'):
        CalcBuriedVolume(mol, [0, 1], center=0)
    mol.GetConformer().SetAtomPosition(1, (8., 0, 0))
    assert CalcBuriedVolume(mol, [1], center=0)['values']['vbur_percent'] == 0.


def test_buried_volume_partial_overlap_convergence_and_motion():
    mol = molecule('[He].[He]', [(0, 0, 0), (1, 0, 0)])
    # Equal spheres R=1 separated by d=1: intersection=5*pi/12,
    # divided by sphere volume 4*pi/3 gives 5/16.
    for spacing, tol in ((.1, .8), (.04, .3)):
        for candidate in (mol, moved(mol)):
            report = CalcBuriedVolume(candidate, [1], center=0, radii_angstrom=[1, 1],
                                       sphere_radius_angstrom=1., spacing_angstrom=spacing)
            assert report['values']['vbur_percent'] == pytest.approx(31.25, abs=tol)


def test_planes_angle_distance_and_degenerate_missing():
    mol = molecule('[He].[He].[He].[He].[He].[He]',
                   [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (0, 1, 1), (0, 0, 2)])
    for candidate in (mol, moved(mol)):
        report = CalcPlaneDescriptors(candidate, [0, 1, 2], other_atom_indices=[3, 4, 5])
        assert report['values']['plane_angle_deg'] == pytest.approx(90.)
        assert report['values']['plane_rms_a'] == pytest.approx(0., abs=1e-12)
        assert report['values']['centroid_distance_a'] == pytest.approx(np.sqrt(17)/3)
    mol.GetConformer().SetAtomPosition(2, (2, 0, 0))
    report = CalcPlaneDescriptors(mol, [0, 1, 2])
    assert report['values']['plane_rms_a'] is None
    assert report['missing']['plane_rms_a'] == 'plane_normal_not_unique'


def test_cremer_pople_planar_chair_and_fourier_mode():
    assert CalcRingPuckering(ring(), range(6))['values']['q_total_a'] == pytest.approx(0.)
    assert CalcRingPuckering(ring(), range(6))['values']['theta_deg'] is None
    for n in (4, 5, 6, 7, 8):
        heights = .2*np.cos(4*np.pi*np.arange(n)/n)
        for mol in (ring(n, heights=heights), moved(ring(n, heights=heights))):
            report = CalcRingPuckering(mol, range(n))
            assert report['values']['q_total_a'] == pytest.approx(np.linalg.norm(heights), abs=1e-12)
            mode_sum = sum(m['q_a']**2 for m in report['modes'])
            mode_sum += report['values'].get('q_alternating_a', 0)**2
            assert mode_sum == pytest.approx(np.dot(heights, heights), abs=1e-12)
    chair = ring(heights=.3*(-1.)**np.arange(6))
    report = CalcRingPuckering(chair, range(6))
    assert report['values']['q_total_a'] == pytest.approx(np.sqrt(6)*.3)
    assert report['modes'][0]['phi_deg'] is None
    assert report['values']['theta_deg'] == pytest.approx(180.)


def test_conformer_selection_input_validation_and_json():
    mol = ring()
    for fn, kwargs in [(CalcHOMA, {}), (CalcRingPuckering, {}), (CalcPlaneDescriptors, {})]:
        with pytest.raises(ValueError):
            fn(mol, range(6), conf_id=True, **kwargs)
        report = fn(mol, range(6), conf_id=7, **kwargs)
        assert json.loads(json.dumps(report, allow_nan=False)) == report
    with pytest.raises(ValueError, match='radii'):
        CalcBuriedVolume(mol, [1], center=0, radii_angstrom=[1.])
    bad = Chem.Mol(mol)
    bad.GetConformer().SetAtomPosition(0, (float('nan'), 0., 0.))
    with pytest.raises(ValueError, match='finite'):
        CalcHOMA(bad, range(6))


def test_fused_rings_conformer_ids_and_missing_puckering_schema():
    mol = molecule('c1ccc2ccccc2c1', [(i, i % 3, 0.) for i in range(10)])
    rings = mol.GetRingInfo().AtomRings()
    assert set(rings[0]) & set(rings[1])
    for indices in rings:
        assert CalcHOMA(mol, indices)['values']['homa'] is not None
    mol = ring()
    conf = Chem.Conformer(mol.GetConformer())
    conf.SetId(12)
    for i, p in enumerate(conf.GetPositions()):
        conf.SetAtomPosition(i, p * 1.1)
    mol.AddConformer(conf, assignId=False)
    assert CalcHOMA(mol, range(6), conf_id=7)['values']['homa'] == pytest.approx(1.)
    assert CalcHOMA(mol, range(6), conf_id=12)['values']['homa'] < 1.
    good = CalcRingPuckering(mol, range(6))
    for i in range(mol.GetNumAtoms()):
        mol.GetConformer().SetAtomPosition(i, (i, 0, 0))
    degenerate = CalcRingPuckering(mol, range(6))
    assert degenerate['values'].keys() == good['values'].keys()
    assert all(v is None for v in degenerate['values'].values())
    assert degenerate['modes'][0]['q_a'] is None
    json.dumps(degenerate, allow_nan=False)
