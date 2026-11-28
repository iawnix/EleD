"""Tests for bond electric fields, real-space density/ALIE/LEA, and 0.001 a.u. isosurface."""
from __future__ import annotations

from dataclasses import fields
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (
    CalcBondElectricField,
    CalcBondElectricFieldFromESP,
    CalcDensityIsosurface,
    CalcRealSpaceProperties,
    CalcRealSpaceSurfaceDescriptors,
    GenerateESPSurfaceGrid,
    XTBCalculator,
    XTBOptions,
    XTBResults,
)
from eled.electric_field import AU_EFIELD_TO_MV_CM, _bond_stencil
from eled.errors import DescriptorUnavailableError
from eled.esp import ANGSTROM_TO_BOHR
from eled.input import _input_metadata, _xyz_from_mol
from eled.runtime import _find_xtb_library, _native_realspace_build


def _make_mol(smiles: str, seed: int = 17):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _make_single_atom(z: int = 10):
    rw = Chem.RWMol()
    a = Chem.Atom(z)
    a.SetNoImplicit(True)
    rw.AddAtom(a)
    mol = rw.GetMol()
    mol.UpdatePropertyCache()
    conf = Chem.Conformer(1)
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def _rotate_translate_mol(mol, angle_rad: float = 0.85, shift=(1.5, -2.0, 0.7)):
    c, s = np.cos(angle_rad), np.sin(angle_rad)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    shift_arr = np.asarray(shift, dtype=float)
    out = Chem.Mol(mol)
    conf = out.GetConformer()
    for i in range(out.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        pos = rot @ np.array([p.x, p.y, p.z], dtype=float) + shift_arr
        conf.SetAtomPosition(i, pos.tolist())
    return out, rot, shift_arr


# ---------------------------------------------------------------------------
# Pure-Python unit tests: Bond Electric Field & Real-Space Surface Analysis
# ---------------------------------------------------------------------------


def test_bond_electric_field_synthetic_uniform_field_and_reversal():
    mol = _make_mol('O', seed=17)
    step = 1e-3
    bonds_meta, stencil = _bond_stencil(mol, mol.GetConformer().GetId(), step)
    # Uniform electric field E0 (a.u.): V(r_bohr) = -E0 . r_bohr -> -grad V = E0.
    e0 = np.array([0.03, -0.04, 0.05], dtype=float)
    v_au = -(stencil * ANGSTROM_TO_BOHR) @ e0
    meta = _input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))
    res = XTBResults(
        method='GFN2-xTB',
        atomicCharges=(-0.6, 0.3, 0.3),
        espPointsAngstrom=tuple(map(tuple, stencil.tolist())),
        espValuesAU=tuple(map(float, v_au.tolist())),
        metadata=meta,
    )
    out = CalcBondElectricFieldFromESP(mol, res, stepAngstrom=step)
    assert out['n_bonds'] == 2
    for idx, row in enumerate(out['bonds']):
        np.testing.assert_allclose(row['e_vec_au'], e0, atol=1e-11)
        expected_par = float(np.dot(e0, bonds_meta[idx]['unit_vector']))
        expected_norm = float(np.linalg.norm(e0))
        expected_perp = float(np.linalg.norm(e0 - expected_par * bonds_meta[idx]['unit_vector']))
        assert row['e_parallel_au'] == pytest.approx(expected_par, abs=1e-11)
        assert row['e_parallel_mv_cm'] == pytest.approx(expected_par * AU_EFIELD_TO_MV_CM, rel=1e-10)
        assert row['e_perp_au'] == pytest.approx(expected_perp, abs=1e-11)
        assert row['e_norm_au'] == pytest.approx(expected_norm, abs=1e-11)
    json.dumps(out, allow_nan=False)

    # Wrong stencil raises ValueError.
    with pytest.raises(ValueError, match='stencil'):
        CalcBondElectricFieldFromESP(mol, res, stepAngstrom=2e-3)


def test_bond_electric_field_zero_bond_atom_and_guards():
    ne = _make_single_atom(10)
    meta = _input_metadata(ne, -1, 0, 0, _xyz_from_mol(ne, -1))
    res = XTBResults(method='GFN2-xTB', atomicCharges=(0.0,), metadata=meta)
    out = CalcBondElectricFieldFromESP(ne, res)
    assert out['n_bonds'] == 0
    assert out['bonds'] == []
    assert out['missing']['bonds'] == 'molecule_has_no_bonds'

    with pytest.raises(ValueError, match='stepAngstrom'):
        CalcBondElectricField(ne, stepAngstrom=0.0)
    with pytest.raises(ValueError, match='stepAngstrom'):
        CalcBondElectricField(ne, stepAngstrom=0.5)
    with pytest.raises(ValueError, match='interface'):
        CalcBondElectricField(ne, interface='cli')


def test_realspace_options_and_pure_python_cutoff_missing_semantics():
    opt = XTBOptions()
    assert opt.realSpacePointsAngstrom == ()
    assert opt.realSpaceOccupationThreshold == 1e-7
    assert opt.realSpaceLEAWindowEV == 0.0
    d = opt.ToDict()
    assert d['realSpacePointsAngstrom'] == [] and set(d) == {f.name for f in fields(XTBOptions)}
    assert XTBOptions.FromDict(d) == opt

    opt2 = XTBOptions(
        realSpacePointsAngstrom=[[0.0, 0.0, 1.0]],
        realSpaceOccupationThreshold=1e-6,
        realSpaceLEAWindowEV=5.0,
    )
    assert XTBOptions.FromDict(json.loads(json.dumps(opt2.ToDict()))) == opt2

    with pytest.raises(ValueError):
        XTBOptions(realSpaceOccupationThreshold=-1e-3)
    with pytest.raises(ValueError):
        XTBOptions(realSpaceLEAWindowEV=-1.0)

    mol = _make_mol('O', seed=17)
    meta = _input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))
    unsampled = XTBResults(method='GFN2-xTB', atomicCharges=(-0.6, 0.3, 0.3), metadata=meta)
    with pytest.raises(DescriptorUnavailableError):
        CalcRealSpaceProperties(unsampled)

    # Point 0 valid, point 1 density below cutoff (must return None, not 0.0),
    # atom 2 has no assigned points (buried atom -> explicit missing).
    sampled = XTBResults(
        method='GFN2-xTB',
        atomicCharges=(-0.6, 0.3, 0.3),
        realSpacePointsAngstrom=((0.0, 0.0, 1.0), (100.0, 0.0, 0.0)),
        realSpaceDensityAU=(0.01, 0.0),
        realSpaceALIEEV=(12.5, 0.0),
        realSpaceVirtualDensityAU=(0.005, 0.0),
        realSpaceLEAEV=(-1.5, 0.0),
        metadata=meta,
    )
    props = CalcRealSpaceProperties(sampled, densityCutoffAU=1e-18)
    assert props['alie_ev'] == [12.5, None]
    assert props['lea_ev'] == [-1.5, None]
    assert props['missing']['alie_ev'] == 'density_below_cutoff_at_some_points'

    surf = CalcRealSpaceSurfaceDescriptors(mol, sampled, atomIndices=[0, 1], weights=[2.0, 1.0])
    assert surf['alie_mean_ev'] == pytest.approx(12.5)
    assert surf['atoms'][0]['alie_mean_ev'] == pytest.approx(12.5)
    assert surf['atoms'][1]['alie_mean_ev'] is None
    assert surf['atoms'][1]['missing']['alie_mean_ev'] == 'density_below_cutoff'
    assert surf['atoms'][2]['alie_mean_ev'] is None
    assert surf['atoms'][2]['missing']['alie_mean_ev'] == 'no_surface_points_assigned'
    json.dumps(surf, allow_nan=False)


# ---------------------------------------------------------------------------
# Integration tests (require verified libxtb with realspace + ESP extensions)
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def library():
    path = _find_xtb_library(None)
    if not _native_realspace_build(path):
        pytest.skip('verified realspace libxtb build required')
    return path


@pytest.mark.integration
def test_bond_electric_field_integration_parity_and_invariance(library):
    mol = _make_mol('O', seed=17)
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)
    ef_native = CalcBondElectricField(mol, options=opts, interface='native')
    ef_compiled = CalcBondElectricField(mol, options=opts, interface='compiled')

    assert ef_native['n_bonds'] == 2
    for bn, bc in zip(ef_native['bonds'], ef_compiled['bonds']):
        np.testing.assert_allclose(bn['e_vec_au'], bc['e_vec_au'], atol=1e-9)
        assert bn['e_parallel_au'] == pytest.approx(bc['e_parallel_au'], abs=1e-9)
        assert bn['e_norm_au'] > 0.05
    # Both O-H bonds in water are equivalent (at similar embedded bond lengths): same sign on e_parallel.
    assert ef_native['bonds'][0]['e_parallel_au'] * ef_native['bonds'][1]['e_parallel_au'] > 0.0

    # Rigid rotation + translation invariance (central FD O(h^2) + libxtb Cartesian ESP screening).
    mol_rot, rot, _ = _rotate_translate_mol(mol)
    ef_rot = CalcBondElectricField(mol_rot, options=opts, interface='native')
    for b0, b1 in zip(ef_native['bonds'], ef_rot['bonds']):
        assert b1['e_parallel_au'] == pytest.approx(b0['e_parallel_au'], rel=1e-4, abs=1e-4)
        assert b1['e_perp_au'] == pytest.approx(b0['e_perp_au'], rel=1e-4, abs=1e-4)
        assert b1['e_norm_au'] == pytest.approx(b0['e_norm_au'], rel=1e-4, abs=1e-4)
        np.testing.assert_allclose(b1['e_vec_au'], rot @ np.asarray(b0['e_vec_au']), rtol=1e-4, atol=1e-4)


@pytest.mark.integration
def test_realspace_native_compiled_parity_and_surface_alie_lea(library):
    from eled import rdXTBTools

    mol = _make_mol('O', seed=17)
    grid = GenerateESPSurfaceGrid(mol, scales=(1.4,), pointsPerAtom=30)
    pts = grid.pointsAngstrom

    calc = XTBCalculator(method='GFN2-xTB', xtbLibrary=library)
    res_native = calc.Calculate(mol, realSpacePointsAngstrom=pts)
    ok, res_compiled = rdXTBTools.RunMol(
        mol,
        method='GFN2-xTB',
        xtbLibrary=library,
        realSpacePointsAngstrom=pts,
    )
    assert ok
    assert res_native.realSpaceSampled is True
    assert res_compiled.realSpaceSampled is True
    assert res_compiled.realSpaceExportVerified is True

    np.testing.assert_allclose(res_native.GetRealSpaceDensity(), res_compiled.GetRealSpaceDensity(), atol=1e-12)
    np.testing.assert_allclose(res_native.GetRealSpaceALIE(), res_compiled.GetRealSpaceALIE(), atol=1e-12)
    np.testing.assert_allclose(
        res_native.GetRealSpaceVirtualDensity(), res_compiled.GetRealSpaceVirtualDensity(), atol=1e-12
    )
    np.testing.assert_allclose(res_native.GetRealSpaceLEA(), res_compiled.GetRealSpaceLEA(), atol=1e-12)

    surf_native = CalcRealSpaceSurfaceDescriptors(mol, res_native, surfaceGrid=grid)
    surf_compiled = CalcRealSpaceSurfaceDescriptors(mol, res_compiled, surfaceGrid=grid)
    assert surf_native['alie_mean_ev'] == pytest.approx(surf_compiled['alie_mean_ev'], abs=1e-10)
    assert surf_native['lea_mean_ev'] == pytest.approx(surf_compiled['lea_mean_ev'], abs=1e-10)
    assert 8.0 < surf_native['alie_mean_ev'] < 25.0
    assert surf_native['alie_min_ev'] < surf_native['alie_max_ev']
    json.dumps(surf_native, allow_nan=False)


@pytest.mark.integration
def test_density_isosurface_0001_au_volume_area_esp_and_alie(library):
    water = _make_mol('O', seed=17)
    ethanol = _make_mol('CCO', seed=17)
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)

    iso_h2o = CalcDensityIsosurface(
        water, isovalueAU=0.001, pointsPerAtom=80, includeESP=True, includeALIE=True, options=opts
    )
    iso_etoh = CalcDensityIsosurface(
        ethanol, isovalueAU=0.001, pointsPerAtom=60, includeESP=False, includeALIE=False, options=opts
    )

    assert iso_h2o['surface_model'] == 'electron_density_isosurface'
    assert iso_h2o['max_isovalue_error_au'] < 1e-5
    assert 12.0 < iso_h2o['volume_a3'] < 40.0
    assert 20.0 < iso_h2o['area_a2'] < 70.0
    assert 0.70 < iso_h2o['sphericity'] <= 1.05
    assert iso_h2o['esp_surface'] is not None
    assert iso_h2o['realspace_surface'] is not None
    assert iso_h2o['esp_surface']['analysis_settings']['surface_model'] == 'electron_density_isosurface'
    assert iso_h2o['realspace_surface']['analysis_settings']['surface_model'] == 'electron_density_isosurface'
    assert iso_h2o['esp_surface']['v_min_au'] < 0.0 < iso_h2o['esp_surface']['v_max_au']
    assert 8.0 < iso_h2o['realspace_surface']['alie_mean_ev'] < 25.0

    # Larger molecule has strictly larger volume and area.
    assert iso_etoh['volume_a3'] > iso_h2o['volume_a3'] * 2.0
    assert iso_etoh['area_a2'] > iso_h2o['area_a2'] * 1.5
    assert iso_etoh['esp_surface'] is None
    assert iso_etoh['realspace_surface'] is None
    assert iso_etoh['missing']['esp_surface'] == 'not_requested'
    assert iso_etoh['missing']['realspace_surface'] == 'not_requested'

    # Translation invariance (exact).
    water_shifted, _, _ = _rotate_translate_mol(water, angle_rad=0.0, shift=(3.5, -1.2, 2.4))
    iso_shift = CalcDensityIsosurface(
        water_shifted, isovalueAU=0.001, pointsPerAtom=80, includeESP=True, includeALIE=True, options=opts
    )
    assert iso_shift['volume_a3'] == pytest.approx(iso_h2o['volume_a3'], abs=1e-8)
    assert iso_shift['area_a2'] == pytest.approx(iso_h2o['area_a2'], abs=1e-8)
    assert iso_shift['esp_surface']['v_mean_au'] == pytest.approx(
        iso_h2o['esp_surface']['v_mean_au'], abs=1e-8
    )
    assert iso_shift['realspace_surface']['alie_mean_ev'] == pytest.approx(
        iso_h2o['realspace_surface']['alie_mean_ev'], abs=1e-8
    )
    json.dumps(iso_h2o, allow_nan=False)
