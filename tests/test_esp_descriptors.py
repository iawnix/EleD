"""Analytical, invariance, guard, and real-xTB tests for Politzer-Murray ESP surface descriptors."""
from __future__ import annotations

from dataclasses import replace
import json
import os

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from eled import DescriptorUnavailableError, GenerateESPGrid, XTBCalculator, XTBResults
from eled.errors import XTBNativeNotFoundError
from eled.esp_descriptors import (
    HARTREE_TO_KCAL_MOL,
    ESPSurfaceGrid,
    CalcESPSurfaceDescriptors,
    GenerateESPSurfaceGrid,
)
from eled.input import _input_metadata, _xyz_from_mol
from eled.runtime import _find_xtb_library, _native_esp_build


def _make_mol(smiles: str = 'O', seed: int = 42) -> Chem.Mol:
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _make_two_atom_mol() -> Chem.Mol:
    mol = Chem.MolFromSmiles('[O]=[C]')
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, Point3D(2.0, 0.0, 0.0))
    conf.Set3D(True)
    mol.AddConformer(conf, assignId=True)
    return mol


def _make_buried_mol() -> Chem.Mol:
    """Carbon tetrachloride: central C (index 0) is completely buried inside Cl vdW shells."""
    return _make_mol('C(Cl)(Cl)(Cl)Cl', seed=42)


def _synthetic_result(
    mol: Chem.Mol,
    points: np.ndarray | tuple[tuple[float, float, float], ...],
    values: np.ndarray | tuple[float, ...],
    *,
    confId: int = -1,
) -> XTBResults:
    xyz = _xyz_from_mol(mol, confId)
    charge = Chem.GetFormalCharge(mol)
    uhf = sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    meta = _input_metadata(mol, confId, charge, uhf, xyz)
    meta['calculation_status'] = 'success'
    natoms = mol.GetNumAtoms()
    return XTBResults(
        method='GFN2-xTB',
        version='synthetic',
        totalEnergyHartree=-1.0,
        homoLumoGapEV=5.0,
        numOrbitals=max(2, natoms),
        numElectrons=sum(a.GetAtomicNum() for a in mol.GetAtoms()) - charge,
        numUnpairedElectrons=uhf,
        atomicCharges=tuple(0.0 for _ in range(natoms)),
        espPointsAngstrom=tuple(tuple(map(float, row)) for row in points),
        espValuesAU=tuple(float(v) for v in values),
        metadata=meta,
    )


def _transform_mol(mol: Chem.Mol, rotation: np.ndarray, translation: np.ndarray, confId: int = -1) -> Chem.Mol:
    copy = Chem.Mol(mol)
    conf = copy.GetConformer(confId)
    coords = np.asarray(conf.GetPositions(), dtype=float)
    moved = coords @ rotation.T + translation
    for i, row in enumerate(moved):
        conf.SetAtomPosition(i, Point3D(float(row[0]), float(row[1]), float(row[2])))
    return copy


def _require_esp_libxtb():
    try:
        library = _find_xtb_library(None)
    except XTBNativeNotFoundError:
        if os.environ.get('XTB_LIBRARY'):
            raise
        pytest.skip('extended libxtb is not installed')
    if not _native_esp_build(library):
        pytest.skip('verified libxtb ESP extension is not installed')
    return library


def test_surface_grid_matches_generate_esp_grid_area_and_ownership():
    mol = _make_mol('O')
    grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=60)
    plain = GenerateESPGrid(mol, scales=(1.4,), pointsPerAtom=60)

    assert isinstance(grid, ESPSurfaceGrid)
    assert grid.pointsAngstrom == plain
    assert len(grid.atomIndices) == len(plain)
    assert len(grid.areaWeightsAngstrom2) == len(plain)
    assert len(grid.shellScales) == len(plain)
    assert grid.surfaceModel == 'rdkit_vdw_fibonacci_shell'
    assert '0.001 a.u.' in grid.settings['surface_definition']
    assert set(grid.atomIndices).issubset(set(range(mol.GetNumAtoms())))
    assert all(w > 0 for w in grid.areaWeightsAngstrom2)
    assert all(s == pytest.approx(1.4) for s in grid.shellScales)

    # Owned copies do not mutate grid state.
    pts = grid.GetPoints()
    pts[:] = 999.0
    assert grid.pointsAngstrom == plain
    owners = grid.GetAtomIndices()
    owners[:] = -1
    assert all(0 <= idx < mol.GetNumAtoms() for idx in grid.atomIndices)
    weights = grid.GetAreaWeights()
    weights[:] = 0.0
    assert all(w > 0 for w in grid.areaWeightsAngstrom2)

    # JSON serializable without NaN.
    dumped = json.loads(json.dumps(grid.ToDict(), allow_nan=False))
    assert dumped['surface_model'] == 'rdkit_vdw_fibonacci_shell'
    assert dumped['units']['area_weights_a2'] == 'angstrom^2'

    # Exact total area for an isolated single atom sphere.
    neon = Chem.MolFromSmiles('[Ne]')
    conf = Chem.Conformer(1)
    conf.SetAtomPosition(0, Point3D(0.0, 0.0, 0.0))
    conf.Set3D(True)
    neon.AddConformer(conf, assignId=True)
    r_ne = Chem.GetPeriodicTable().GetRvdw(10)
    ne_grid = GenerateESPSurfaceGrid(neon, scales=(1.4, 1.8), pointsPerAtom=80)
    expected_area = 4.0 * np.pi * ((1.4 * r_ne) ** 2 + (1.8 * r_ne) ** 2)
    assert len(ne_grid.pointsAngstrom) == 160
    assert sum(ne_grid.areaWeightsAngstrom2) == pytest.approx(expected_area, rel=1e-12)


def test_surface_grid_respects_conf_id_and_rejects_bad_inputs():
    mol = _make_mol('O')
    conf1 = Chem.Conformer(mol.GetConformer(0))
    for i in range(mol.GetNumAtoms()):
        p = conf1.GetAtomPosition(i)
        conf1.SetAtomPosition(i, Point3D(p.x + 5.0, p.y, p.z))
    cid = mol.AddConformer(conf1, assignId=True)
    g0 = GenerateESPSurfaceGrid(mol, 0, scales=(1.4,), pointsPerAtom=20)
    g1 = GenerateESPSurfaceGrid(mol, cid, scales=(1.4,), pointsPerAtom=20)
    assert g1.settings['conf_id'] == cid
    np.testing.assert_allclose(
        g1.GetPoints() - g0.GetPoints(),
        np.broadcast_to(np.array([5.0, 0.0, 0.0]), g0.GetPoints().shape),
        atol=1e-12,
    )

    with pytest.raises(TypeError):
        GenerateESPSurfaceGrid(None)
    for kwargs in (
        {'confId': True},
        {'confId': -2},
        {'scales': ()},
        {'scales': (1.0,)},
        {'scales': (1.4, 1.4)},
        {'scales': (np.nan,)},
        {'scales': '1.4'},
        {'pointsPerAtom': True},
        {'pointsPerAtom': 5},
        {'pointsPerAtom': 100001},
    ):
        with pytest.raises(ValueError):
            GenerateESPSurfaceGrid(mol, **kwargs)


def test_analytical_politzer_murray_unweighted_and_weighted():
    mol = _make_two_atom_mol()
    # Two points near O (atom 0 at x=0) with negative V, two points near C (atom 1 at x=2) with positive V.
    points = np.array([
        [-2.0, 0.0, 0.0],
        [-2.5, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.5, 0.0, 0.0],
    ])
    # V_neg = [-0.04, -0.02] -> mean = -0.03, sigma2_neg = 0.0001
    # V_pos = [ 0.02,  0.06] -> mean =  0.04, sigma2_pos = 0.0004
    # V_mean = 0.005
    # Pi = mean(|[-0.045, -0.025, 0.015, 0.055]|) = 0.140 / 4 = 0.035
    # sigma2_tot = 0.0005
    # nu = (0.0004 * 0.0001) / (0.0005 ** 2) = 0.16
    # nu_sigma2 = 0.16 * 0.0005 = 0.00008
    values = np.array([-0.04, -0.02, 0.02, 0.06])
    res = _synthetic_result(mol, points, values)

    desc = CalcESPSurfaceDescriptors(mol, res)
    assert desc['n_points'] == 4
    assert desc['n_pos_points'] == 2
    assert desc['n_neg_points'] == 2
    assert desc['v_min_au'] == pytest.approx(-0.04, abs=1e-15)
    assert desc['v_max_au'] == pytest.approx(0.06, abs=1e-15)
    assert desc['v_mean_au'] == pytest.approx(0.005, abs=1e-15)
    assert desc['v_pos_mean_au'] == pytest.approx(0.04, abs=1e-15)
    assert desc['v_neg_mean_au'] == pytest.approx(-0.03, abs=1e-15)
    assert desc['v_pos_frac'] == pytest.approx(0.5, abs=1e-15)
    assert desc['v_neg_frac'] == pytest.approx(0.5, abs=1e-15)
    assert desc['pi_esp_au'] == pytest.approx(0.035, abs=1e-15)
    assert desc['sigma2_pos_au2'] == pytest.approx(0.0004, abs=1e-15)
    assert desc['sigma2_neg_au2'] == pytest.approx(0.0001, abs=1e-15)
    assert desc['sigma2_tot_au2'] == pytest.approx(0.0005, abs=1e-15)
    assert desc['nu'] == pytest.approx(0.16, abs=1e-15)
    assert desc['nu_sigma2_au2'] == pytest.approx(0.00008, abs=1e-15)

    # Matching kcal/mol conversions.
    c = HARTREE_TO_KCAL_MOL
    assert desc['v_max_kcalmol'] == pytest.approx(0.06 * c, rel=1e-12)
    assert desc['v_min_kcalmol'] == pytest.approx(-0.04 * c, rel=1e-12)
    assert desc['v_mean_kcalmol'] == pytest.approx(0.005 * c, rel=1e-12)
    assert desc['v_pos_mean_kcalmol'] == pytest.approx(0.04 * c, rel=1e-12)
    assert desc['v_neg_mean_kcalmol'] == pytest.approx(-0.03 * c, rel=1e-12)
    assert desc['pi_esp_kcalmol'] == pytest.approx(0.035 * c, rel=1e-12)
    assert desc['sigma2_pos_kcalmol2'] == pytest.approx(0.0004 * c * c, rel=1e-12)
    assert desc['sigma2_neg_kcalmol2'] == pytest.approx(0.0001 * c * c, rel=1e-12)
    assert desc['sigma2_tot_kcalmol2'] == pytest.approx(0.0005 * c * c, rel=1e-12)
    assert desc['nu_sigma2_kcalmol2'] == pytest.approx(0.00008 * c * c, rel=1e-12)

    # Per-atom checks in original RDKit order.
    assert len(desc['atoms']) == 2
    a0, a1 = desc['atoms']
    assert a0['atom_index'] == 0 and a0['n_points'] == 2
    assert a0['v_min_au'] == pytest.approx(-0.04)
    assert a0['v_max_au'] == pytest.approx(-0.02)
    assert a0['v_mean_au'] == pytest.approx(-0.03)
    assert a0['v_neg_mean_au'] == pytest.approx(-0.03)
    assert a0['v_pos_mean_au'] is None
    assert a0['v_pos_mean_kcalmol'] is None
    assert a0['sigma2_au2'] == pytest.approx(0.0001)
    assert a0['missing'] == {
        'v_pos_mean_au': 'no_positive_surface_points',
        'v_pos_mean_kcalmol': 'no_positive_surface_points',
    }

    assert a1['atom_index'] == 1 and a1['n_points'] == 2
    assert a1['v_mean_au'] == pytest.approx(0.04)
    assert a1['v_pos_mean_au'] == pytest.approx(0.04)
    assert a1['v_neg_mean_au'] is None
    assert a1['sigma2_au2'] == pytest.approx(0.0004)
    assert a1['missing']['v_neg_mean_au'] == 'no_negative_surface_points'

    # Weighted analytical check: weights = [1, 3, 2, 2] -> W=8, W_neg=4, W_pos=4.
    # V_neg weighted mean = (1*(-0.04) + 3*(-0.02))/4 = -0.025
    # sigma2_neg = (1*(-0.015)^2 + 3*(0.005)^2)/4 = (0.000225 + 0.000075)/4 = 0.000075
    # V_pos weighted mean = 0.04, sigma2_pos = 0.0004
    # V_mean = (4*(-0.025) + 4*(0.04))/8 = 0.0075
    weights = np.array([1.0, 3.0, 2.0, 2.0])
    wdesc = CalcESPSurfaceDescriptors(mol, res, weights=weights, includeKcalMol=False)
    assert 'v_max_kcalmol' not in wdesc
    assert 'v_max_kcalmol' not in wdesc['units']
    assert wdesc['v_mean_au'] == pytest.approx(0.0075, abs=1e-15)
    assert wdesc['v_neg_mean_au'] == pytest.approx(-0.025, abs=1e-15)
    assert wdesc['sigma2_neg_au2'] == pytest.approx(0.000075, abs=1e-15)
    assert wdesc['sigma2_tot_au2'] == pytest.approx(0.000475, abs=1e-15)
    assert wdesc['nu'] == pytest.approx((0.0004 * 0.000075) / (0.000475 ** 2), abs=1e-15)
    assert wdesc['atoms'][0]['v_mean_au'] == pytest.approx(-0.025, abs=1e-15)
    assert wdesc['atoms'][0]['sigma2_au2'] == pytest.approx(0.000075, abs=1e-15)
    assert json.loads(json.dumps(desc, allow_nan=False)) == desc


def test_balance_maximum_and_zero_potential_points():
    mol = _make_two_atom_mol()
    points = np.array([
        [-2.0, 0.0, 0.0],
        [-2.5, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [4.5, 0.0, 0.0],
        [1.0, 3.0, 0.0],
    ])
    # Equal positive and negative variances -> nu == 0.25 exactly; one V=0 point.
    values = np.array([-0.03, -0.01, 0.01, 0.03, 0.0])
    res = _synthetic_result(mol, points, values)
    desc = CalcESPSurfaceDescriptors(mol, res)
    assert desc['v_pos_frac'] == pytest.approx(0.4)
    assert desc['v_neg_frac'] == pytest.approx(0.4)
    assert desc['sigma2_pos_au2'] == pytest.approx(0.0001)
    assert desc['sigma2_neg_au2'] == pytest.approx(0.0001)
    assert desc['nu'] == pytest.approx(0.25, abs=1e-15)
    assert desc['nu_sigma2_au2'] == pytest.approx(0.25 * 0.0002, abs=1e-15)


@pytest.mark.parametrize(
    'values,expected_reason,expect_pos,expect_neg,expect_tot',
    [
        ([0.02, 0.06, 0.04, 0.08], 'single_sign_surface', True, False, True),
        ([-0.02, -0.06, -0.04, -0.08], 'single_sign_surface', False, True, True),
        ([0.05, 0.05, -0.05, -0.05], 'zero_total_variance', True, True, True),
        ([0.0, 0.0, 0.0, 0.0], 'no_nonzero_surface_points', False, False, False),
    ],
)
def test_edge_cases_single_sign_zero_variance_and_all_zero(
    values, expected_reason, expect_pos, expect_neg, expect_tot
):
    mol = _make_two_atom_mol()
    points = np.array([[-2.0, 0.0, 0.0], [-2.5, 0.0, 0.0], [4.0, 0.0, 0.0], [4.5, 0.0, 0.0]])
    res = _synthetic_result(mol, points, np.array(values))
    desc = CalcESPSurfaceDescriptors(mol, res)

    assert desc['nu'] is None
    assert desc['nu_sigma2_au2'] is None
    assert desc['nu_sigma2_kcalmol2'] is None
    assert desc['missing']['nu'] == expected_reason
    assert desc['missing']['nu_sigma2_au2'] == expected_reason
    assert desc['missing']['nu_sigma2_kcalmol2'] == expected_reason
    assert (desc['v_pos_mean_au'] is not None) is expect_pos
    assert (desc['sigma2_pos_au2'] is not None) is expect_pos
    assert (desc['v_neg_mean_au'] is not None) is expect_neg
    assert (desc['sigma2_neg_au2'] is not None) is expect_neg
    assert (desc['sigma2_tot_au2'] is not None) is expect_tot
    if not expect_pos:
        assert desc['missing']['v_pos_mean_au'] == 'no_positive_surface_points'
    if not expect_neg:
        assert desc['missing']['v_neg_mean_au'] == 'no_negative_surface_points'
    if not expect_tot:
        assert desc['missing']['sigma2_tot_au2'] == 'no_nonzero_surface_points'
    json.dumps(desc, allow_nan=False)


def test_buried_atom_and_zero_weight_atom_return_none_with_reasons():
    mol = _make_buried_mol()
    grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=30)
    assert 0 not in grid.atomIndices

    rng = np.random.default_rng(7)
    values = rng.normal(0.0, 0.05, size=len(grid.pointsAngstrom))
    res = _synthetic_result(mol, grid.pointsAngstrom, values)

    desc_grid = CalcESPSurfaceDescriptors(mol, res, surfaceGrid=grid)
    desc_explicit = CalcESPSurfaceDescriptors(
        mol, res, atomIndices=grid.atomIndices, weights=grid.areaWeightsAngstrom2
    )
    assert desc_grid['v_mean_au'] == pytest.approx(desc_explicit['v_mean_au'], abs=1e-15)
    assert desc_grid['sigma2_tot_au2'] == pytest.approx(desc_explicit['sigma2_tot_au2'], abs=1e-15)

    buried = desc_grid['atoms'][0]
    assert buried['atom_index'] == 0
    assert buried['n_points'] == 0
    for key in ('v_min_au', 'v_max_au', 'v_mean_au', 'v_pos_mean_au', 'v_neg_mean_au', 'sigma2_au2'):
        assert buried[key] is None
        assert buried['missing'][key] == 'no_surface_points_assigned'
    assert desc_grid['missing']['atoms']['0'] == buried['missing']

    # Also verify zero-weight assigned atom returns None with zero_surface_weight_assigned.
    mol2 = _make_two_atom_mol()
    pts2 = np.array([[-2.0, 0.0, 0.0], [4.0, 0.0, 0.0], [4.5, 0.0, 0.0]])
    res2 = _synthetic_result(mol2, pts2, np.array([-0.5, 0.02, 0.06]))
    desc2 = CalcESPSurfaceDescriptors(mol2, res2, atomIndices=[0, 1, 1], weights=[0.0, 1.0, 1.0])
    assert desc2['n_points'] == 2
    assert desc2['n_points_total'] == 3
    assert desc2['v_min_au'] == pytest.approx(0.02)
    assert desc2['atoms'][0]['n_points'] == 0
    assert desc2['atoms'][0]['v_mean_au'] is None
    assert desc2['atoms'][0]['missing']['v_mean_au'] == 'zero_surface_weight_assigned'


def test_rigid_motion_point_permutation_and_atom_renumbering_invariance():
    mol = _make_mol('C=O', seed=11)
    grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=40)
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    points = grid.GetPoints()
    # Analytical Coulomb-like potential from fixed charges (strictly rigid-motion invariant).
    charges = np.array([0.4, -0.5, 0.05, 0.05])
    dists = np.linalg.norm(points[:, None, :] - coords[None, :, :], axis=2)
    values = (1.0 / dists) @ charges

    res = _synthetic_result(mol, points, values)
    base_nearest = CalcESPSurfaceDescriptors(mol, res)
    base_grid = CalcESPSurfaceDescriptors(mol, res, surfaceGrid=grid)

    # 1. Rigid rotation + translation with co-transformed grid.
    theta = 0.73
    rot_z = np.array([[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])
    phi = -0.41
    rot_x = np.array([[1.0, 0.0, 0.0], [0.0, np.cos(phi), -np.sin(phi)], [0.0, np.sin(phi), np.cos(phi)]])
    rotation = rot_z @ rot_x
    translation = np.array([3.5, -2.1, 7.4])
    moved_mol = _transform_mol(mol, rotation, translation)
    moved_points = points @ rotation.T + translation
    moved_res = _synthetic_result(moved_mol, moved_points, values)

    moved_nearest = CalcESPSurfaceDescriptors(moved_mol, moved_res)
    moved_explicit = CalcESPSurfaceDescriptors(
        moved_mol, moved_res, atomIndices=grid.atomIndices, weights=grid.areaWeightsAngstrom2
    )
    for key in (
        'v_max_au', 'v_min_au', 'v_mean_au', 'v_pos_mean_au', 'v_neg_mean_au',
        'v_pos_frac', 'v_neg_frac', 'pi_esp_au', 'sigma2_pos_au2', 'sigma2_neg_au2',
        'sigma2_tot_au2', 'nu', 'nu_sigma2_au2',
    ):
        assert moved_nearest[key] == pytest.approx(base_nearest[key], abs=1e-12)
        assert moved_explicit[key] == pytest.approx(base_grid[key], abs=1e-12)
    for i in range(mol.GetNumAtoms()):
        assert moved_nearest['atoms'][i]['n_points'] == base_nearest['atoms'][i]['n_points']
        assert moved_nearest['atoms'][i]['v_mean_au'] == pytest.approx(base_nearest['atoms'][i]['v_mean_au'], abs=1e-12)
        assert moved_explicit['atoms'][i]['v_mean_au'] == pytest.approx(base_grid['atoms'][i]['v_mean_au'], abs=1e-12)

    # 2. Point permutation invariance.
    perm_pts = np.random.default_rng(19).permutation(len(points))
    perm_res = _synthetic_result(mol, points[perm_pts], values[perm_pts])
    perm_desc = CalcESPSurfaceDescriptors(
        mol,
        perm_res,
        atomIndices=grid.GetAtomIndices()[perm_pts],
        weights=grid.GetAreaWeights()[perm_pts],
    )
    for key in ('v_mean_au', 'pi_esp_au', 'sigma2_tot_au2', 'nu', 'nu_sigma2_au2'):
        assert perm_desc[key] == pytest.approx(base_grid[key], abs=1e-12)
    for i in range(mol.GetNumAtoms()):
        assert perm_desc['atoms'][i]['v_mean_au'] == pytest.approx(base_grid['atoms'][i]['v_mean_au'], abs=1e-12)

    # 3. RDKit atom-renumbering invariance (new_idx -> old_idx = atom_perm[new_idx]).
    atom_perm = [2, 0, 3, 1]
    inv_perm = {old: new for new, old in enumerate(atom_perm)}
    renumbered_mol = Chem.RenumberAtoms(mol, atom_perm)
    renumbered_res = _synthetic_result(renumbered_mol, points, values)
    renum_nearest = CalcESPSurfaceDescriptors(renumbered_mol, renumbered_res)
    renum_explicit = CalcESPSurfaceDescriptors(
        renumbered_mol,
        renumbered_res,
        atomIndices=[inv_perm[idx] for idx in grid.atomIndices],
        weights=grid.areaWeightsAngstrom2,
    )
    for key in ('v_max_au', 'v_min_au', 'v_mean_au', 'pi_esp_au', 'sigma2_tot_au2', 'nu'):
        assert renum_nearest[key] == pytest.approx(base_nearest[key], abs=1e-12)
        assert renum_explicit[key] == pytest.approx(base_grid[key], abs=1e-12)
    for new_idx, old_idx in enumerate(atom_perm):
        assert renum_nearest['atoms'][new_idx]['n_points'] == base_nearest['atoms'][old_idx]['n_points']
        assert renum_nearest['atoms'][new_idx]['v_mean_au'] == pytest.approx(
            base_nearest['atoms'][old_idx]['v_mean_au'], abs=1e-12
        )
        assert renum_explicit['atoms'][new_idx]['v_mean_au'] == pytest.approx(
            base_grid['atoms'][old_idx]['v_mean_au'], abs=1e-12
        )


def test_invalid_input_guards():
    mol = _make_two_atom_mol()
    grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=12)
    values = np.linspace(-0.02, 0.03, len(grid.pointsAngstrom))
    res = _synthetic_result(mol, grid.pointsAngstrom, values)

    with pytest.raises(TypeError):
        CalcESPSurfaceDescriptors(None, res)
    with pytest.raises(TypeError):
        CalcESPSurfaceDescriptors(mol, mol)
    with pytest.raises(ValueError, match='includeKcalMol'):
        CalcESPSurfaceDescriptors(mol, res, includeKcalMol='yes')
    with pytest.raises(ValueError, match='failed'):
        CalcESPSurfaceDescriptors(mol, replace(res, metadata={**res.metadata, 'calculation_status': 'error'}))
    with pytest.raises(ValueError, match='provenance'):
        CalcESPSurfaceDescriptors(mol, replace(res, metadata={'calculation_status': 'success'}))

    # Geometry mismatch.
    moved = _transform_mol(mol, np.eye(3), np.array([1.0, 0.0, 0.0]))
    with pytest.raises(ValueError, match='coordinates'):
        CalcESPSurfaceDescriptors(moved, res)

    # Atom order mismatch.
    reordered = Chem.RenumberAtoms(mol, [1, 0])
    with pytest.raises(ValueError, match='structure or atom order'):
        CalcESPSurfaceDescriptors(reordered, res)

    # Unsampled vs malformed ESP.
    with pytest.raises(DescriptorUnavailableError):
        CalcESPSurfaceDescriptors(mol, replace(res, espPointsAngstrom=(), espValuesAU=()))
    with pytest.raises(ValueError):
        CalcESPSurfaceDescriptors(mol, replace(res, espValuesAU=(np.nan,) * len(values)))
    with pytest.raises(ValueError):
        CalcESPSurfaceDescriptors(mol, replace(res, espValuesAU=(0.1,)))

    # Bad atomIndices and weights.
    n = len(values)
    for bad_idx in ([0, 1], [True] * n, [0] + [True] * (n - 1), [0.0] * n, [-1] * n, [2] * n, '0' * n):
        with pytest.raises(ValueError):
            CalcESPSurfaceDescriptors(mol, res, atomIndices=bad_idx)
    for bad_w in ([1.0], [1.0] + [True] * (n - 1), [-1.0] * n, [0.0] * n, [np.nan] * n):
        with pytest.raises(ValueError):
            CalcESPSurfaceDescriptors(mol, res, weights=bad_w)

    # Bad surfaceGrid combinations.
    with pytest.raises(TypeError):
        CalcESPSurfaceDescriptors(mol, res, surfaceGrid=grid.pointsAngstrom)
    with pytest.raises(ValueError, match='not both'):
        CalcESPSurfaceDescriptors(mol, res, surfaceGrid=grid, atomIndices=grid.atomIndices)
    with pytest.raises(ValueError, match='not both'):
        CalcESPSurfaceDescriptors(mol, res, surfaceGrid=grid, weights=grid.areaWeightsAngstrom2)
    other_grid = GenerateESPSurfaceGrid(mol, scales=(1.6,), pointsPerAtom=12)
    with pytest.raises(ValueError, match='do not match'):
        CalcESPSurfaceDescriptors(mol, res, surfaceGrid=other_grid)


@pytest.mark.integration
def test_real_xtb_polar_ion_buried_and_invariance():
    _require_esp_libxtb()
    calc = XTBCalculator(method='GFN2-xTB')

    # 1. Polar neutral molecule (water): both positive & negative surface regions.
    water = _make_mol('O', seed=42)
    grid = GenerateESPSurfaceGrid(water, scales=(1.4,), pointsPerAtom=50)
    res = calc.Calculate(water, espPointsAngstrom=grid.pointsAngstrom)
    desc = CalcESPSurfaceDescriptors(water, res, surfaceGrid=grid)
    assert desc['v_min_au'] < 0.0 < desc['v_max_au']
    assert 0.0 < desc['v_pos_frac'] < 1.0
    assert 0.0 < desc['v_neg_frac'] < 1.0
    assert desc['pi_esp_au'] > 0.0
    assert desc['sigma2_tot_au2'] > 0.0
    assert 0.0 < desc['nu'] <= 0.25
    assert desc['nu_sigma2_au2'] == pytest.approx(desc['nu'] * desc['sigma2_tot_au2'], rel=1e-12)
    # Oxygen (idx 0) has negative mean ESP; hydrogens (idx 1, 2) have positive mean ESP.
    assert desc['atoms'][0]['v_mean_au'] < 0.0
    assert desc['atoms'][1]['v_mean_au'] > 0.0
    assert desc['atoms'][2]['v_mean_au'] > 0.0
    json.dumps(desc, allow_nan=False)

    # 2. Rigid rotation + translation invariance with real xTB ESP on co-rotated grid.
    rotation = np.array([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    translation = np.array([2.5, -1.3, 4.2])
    moved_water = _transform_mol(water, rotation, translation)
    moved_points = grid.GetPoints() @ rotation.T + translation
    moved_res = calc.Calculate(moved_water, espPointsAngstrom=moved_points)
    moved_desc = CalcESPSurfaceDescriptors(
        moved_water, moved_res, atomIndices=grid.atomIndices, weights=grid.areaWeightsAngstrom2
    )
    for key in ('v_max_au', 'v_min_au', 'v_mean_au', 'pi_esp_au', 'sigma2_tot_au2', 'nu', 'nu_sigma2_au2'):
        assert moved_desc[key] == pytest.approx(desc[key], abs=1e-9)
    for i in range(water.GetNumAtoms()):
        assert moved_desc['atoms'][i]['v_mean_au'] == pytest.approx(desc['atoms'][i]['v_mean_au'], abs=1e-9)

    # 3. Single-sign ion ([NH4+]): all positive surface -> nu is None.
    ammonium = _make_mol('[NH4+]', seed=42)
    ion_grid = GenerateESPSurfaceGrid(ammonium, scales=(1.4,), pointsPerAtom=40)
    ion_res = calc.Calculate(ammonium, espPointsAngstrom=ion_grid.pointsAngstrom)
    ion_desc = CalcESPSurfaceDescriptors(ammonium, ion_res, surfaceGrid=ion_grid)
    assert ion_desc['v_min_au'] > 0.0
    assert ion_desc['v_pos_frac'] == pytest.approx(1.0)
    assert ion_desc['v_neg_frac'] == pytest.approx(0.0)
    assert ion_desc['v_neg_mean_au'] is None
    assert ion_desc['sigma2_neg_au2'] is None
    assert ion_desc['nu'] is None
    assert ion_desc['nu_sigma2_au2'] is None
    assert ion_desc['missing']['nu'] == 'single_sign_surface'

    # 4. Real molecule with a buried central atom (CCl4 central carbon).
    ccl4 = _make_mol('C(Cl)(Cl)(Cl)Cl', seed=42)
    ccl4_grid = GenerateESPSurfaceGrid(ccl4, scales=(1.4,), pointsPerAtom=40)
    ccl4_res = calc.Calculate(ccl4, espPointsAngstrom=ccl4_grid.pointsAngstrom)
    ccl4_desc = CalcESPSurfaceDescriptors(ccl4, ccl4_res, surfaceGrid=ccl4_grid)
    assert ccl4_desc['atoms'][0]['n_points'] == 0
    assert ccl4_desc['atoms'][0]['v_mean_au'] is None
    assert ccl4_desc['atoms'][0]['missing']['v_mean_au'] == 'no_surface_points_assigned'

    # 5. Real xTB atom-renumbering invariance on formaldehyde.
    hcho = _make_mol('C=O', seed=42)
    hcho_grid = GenerateESPSurfaceGrid(hcho, scales=(1.4,), pointsPerAtom=40)
    hcho_res = calc.Calculate(hcho, espPointsAngstrom=hcho_grid.pointsAngstrom)
    hcho_desc = CalcESPSurfaceDescriptors(hcho, hcho_res)
    perm = [2, 1, 3, 0]
    hcho_renum = Chem.RenumberAtoms(hcho, perm)
    hcho_renum_res = calc.Calculate(hcho_renum, espPointsAngstrom=hcho_grid.pointsAngstrom)
    hcho_renum_desc = CalcESPSurfaceDescriptors(hcho_renum, hcho_renum_res)
    for key in ('v_max_au', 'v_min_au', 'v_mean_au', 'pi_esp_au', 'sigma2_tot_au2', 'nu'):
        assert hcho_renum_desc[key] == pytest.approx(hcho_desc[key], abs=1e-9)
    for new_idx, old_idx in enumerate(perm):
        assert hcho_renum_desc['atoms'][new_idx]['v_mean_au'] == pytest.approx(
            hcho_desc['atoms'][old_idx]['v_mean_au'], abs=1e-9
        )


@pytest.mark.integration
@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['O', '[NH4+]'])
def test_native_and_compiled_esp_surface_parity(method, smiles):
    _require_esp_libxtb()
    rdXTBTools = pytest.importorskip('eled.rdXTBTools')

    mol = _make_mol(smiles, seed=42)
    grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=40)
    native = XTBCalculator(method=method).Calculate(mol, espPointsAngstrom=grid.pointsAngstrom)
    ok, compiled = rdXTBTools.RunMol(mol, method=method, espPointsAngstrom=grid.pointsAngstrom)
    assert ok, compiled.errorMessage

    d_native = CalcESPSurfaceDescriptors(mol, native, surfaceGrid=grid)
    d_compiled = CalcESPSurfaceDescriptors(mol, compiled, surfaceGrid=grid)

    assert d_native['missing'] == d_compiled['missing']
    for key in (
        'v_max_au', 'v_min_au', 'v_mean_au', 'v_pos_mean_au', 'v_neg_mean_au',
        'v_pos_frac', 'v_neg_frac', 'pi_esp_au', 'sigma2_pos_au2', 'sigma2_neg_au2',
        'sigma2_tot_au2', 'nu', 'nu_sigma2_au2',
        'v_max_kcalmol', 'sigma2_tot_kcalmol2', 'nu_sigma2_kcalmol2',
    ):
        if d_native[key] is None:
            assert d_compiled[key] is None
        else:
            assert d_compiled[key] == pytest.approx(d_native[key], abs=1e-10)
    for a_nat, a_cmp in zip(d_native['atoms'], d_compiled['atoms']):
        assert a_nat['atom_index'] == a_cmp['atom_index']
        assert a_nat['n_points'] == a_cmp['n_points']
        assert a_nat['missing'] == a_cmp['missing']
        for key in ('v_min_au', 'v_max_au', 'v_mean_au', 'v_pos_mean_au', 'v_neg_mean_au', 'sigma2_au2'):
            if a_nat[key] is None:
                assert a_cmp[key] is None
            else:
                assert a_cmp[key] == pytest.approx(a_nat[key], abs=1e-10)
