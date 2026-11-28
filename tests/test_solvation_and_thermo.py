"""Tests for Tier 3B implicit solvation (ALPB/GBSA) and vibrational thermochemistry."""
from __future__ import annotations

import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import OptimizeXTB, XTBCacheKey, XTBCalculator, XTBOptions
from eled.errors import XTBEnvironmentError
from eled.runtime import _find_xtb_library, _native_properties_build
from eled.solvation import CalcSolvationDescriptors
from eled.thermochemistry import CalcVibrationalThermochemistry


def _make_mol(smiles: str = 'O', seed: int = 17):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _rotate_translate_permute(mol, perm):
    theta = 0.65
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([-1.4, 2.1, 0.8])
    reordered = Chem.RenumberAtoms(mol, list(perm))
    conf = reordered.GetConformer()
    orig_conf = mol.GetConformer()
    for new_idx, old_idx in enumerate(perm):
        pt = orig_conf.GetAtomPosition(old_idx)
        pos = np.array([pt.x, pt.y, pt.z], dtype=float)
        conf.SetAtomPosition(new_idx, (rot @ pos + shift).tolist())
    return reordered


# ---------------------------------------------------------------------------
# Pure-Python option validation & guard tests
# ---------------------------------------------------------------------------

def test_solvation_options_defaults_roundtrip_and_guards():
    gas = XTBOptions()
    assert gas.solvent is None and gas.solvationModel == 'alpb' and gas.solvationGrid is None
    gas_dict = gas.ToDict()
    assert 'solvent' not in gas_dict and 'solvationModel' not in gas_dict and 'solvationGrid' not in gas_dict
    assert XTBOptions.FromDict(gas_dict) == gas

    solv = XTBOptions(solvent=' water ', solvationModel='GBSA', solvationGrid=110)
    assert solv.solvent == 'water' and solv.solvationModel == 'gbsa' and solv.solvationGrid == 110
    solv_dict = solv.ToDict()
    assert solv_dict['solvent'] == 'water' and solv_dict['solvationModel'] == 'gbsa' and solv_dict['solvationGrid'] == 110
    assert XTBOptions.FromDict(json.loads(json.dumps(solv_dict))) == solv

    for bad in (
        {'solvent': ''},
        {'solvent': '   '},
        {'solvent': 123},
        {'solvationModel': 'pcm'},
        {'solvationModel': ''},
        {'solvationGrid': 100},
        {'solvationGrid': True},
        {'solvationGrid': -1},
    ):
        with pytest.raises(ValueError):
            XTBOptions(**bad)


def test_solvation_and_thermo_input_guards():
    mol = _make_mol('O')
    with pytest.raises(ValueError, match='solvent'):
        CalcSolvationDescriptors(mol, None)
    with pytest.raises(ValueError, match='solvent'):
        CalcSolvationDescriptors(mol, '  ')
    with pytest.raises(ValueError, match='solvationModel'):
        CalcSolvationDescriptors(mol, 'water', solvationModel='smd')
    with pytest.raises(TypeError, match='XTBOptions'):
        CalcSolvationDescriptors(mol, 'water', options={'method': 'GFN2-xTB'})
    for loose in ({'method': 'GFN2-xTB'}, {'xtbLibrary': None}, {'confId': 0}):
        with pytest.raises(TypeError):
            CalcSolvationDescriptors(mol, 'water', **loose)
        with pytest.raises(TypeError):
            CalcVibrationalThermochemistry(mol, **loose)

    for bad in (
        {'temperatureK': 0.0},
        {'temperatureK': -10.0},
        {'pressureAtm': 0.0},
        {'stepBohr': 0.0},
        {'symmetryNumber': 0},
        {'symmetryNumber': True},
        {'rotorCutoffCm1': -5.0},
        {'stationaryGradientThresholdHaBohr': 0.0},
    ):
        with pytest.raises(ValueError):
            CalcVibrationalThermochemistry(mol, **bad)


# ---------------------------------------------------------------------------
# Integration tests (implicit solvation + RRHO/qRRHO thermochemistry)
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def library():
    path = _find_xtb_library(None)
    if not _native_properties_build(path):
        pytest.skip('verified extended properties libxtb build required')
    return path


@pytest.mark.integration
def test_solvation_cache_keys_and_native_compiled_parity(library):
    from eled import rdXTBTools
    mol = _make_mol('O', seed=17)

    # Cache keys: gas-phase default preserved; distinct solvents/models/grids get distinct keys.
    k_gas_default = XTBCacheKey(mol, xtbLibrary=library)
    k_gas_none = XTBCacheKey(mol, solvent=None, xtbLibrary=library)
    k_alpb_water = XTBCacheKey(mol, solvent='water', solvationModel='alpb', xtbLibrary=library)
    k_gbsa_water = XTBCacheKey(mol, solvent='water', solvationModel='gbsa', xtbLibrary=library)
    k_alpb_dmso = XTBCacheKey(mol, solvent='dmso', solvationModel='alpb', xtbLibrary=library)
    k_alpb_grid110 = XTBCacheKey(mol, solvent='water', solvationModel='alpb', solvationGrid=110, xtbLibrary=library)
    assert k_gas_default == k_gas_none
    assert len({k_gas_default, k_alpb_water, k_gbsa_water, k_alpb_dmso, k_alpb_grid110}) == 5

    # Native vs compiled parity for ALPB and GBSA on GFN2 and GFN1
    for method in ('GFN2-xTB', 'GFN1-xTB'):
        for model in ('alpb', 'gbsa'):
            nat = XTBCalculator(
                method=method, solvent='water', solvationModel=model, xtbLibrary=library
            ).Calculate(mol)
            ok, cmp_res = rdXTBTools.RunMol(
                mol, method=method, solvent='water', solvationModel=model, xtbLibrary=library
            )
            assert ok, cmp_res.errorMessage
            assert nat.solvent == cmp_res.solvent == 'water'
            assert nat.solvationModel == cmp_res.solvationModel == model
            assert nat.solvationGrid == cmp_res.solvationGrid == 230
            assert cmp_res.totalEnergyHartree == pytest.approx(nat.totalEnergyHartree, abs=1e-10)
            assert cmp_res.solvationEnergyHartree == pytest.approx(nat.solvationEnergyHartree, abs=1e-10)
            assert cmp_res.solvationBornEnergyHartree == pytest.approx(nat.solvationBornEnergyHartree, abs=1e-10)
            assert cmp_res.solvationSASAEnergyHartree == pytest.approx(nat.solvationSASAEnergyHartree, abs=1e-10)
            assert cmp_res.solvationHBEnergyHartree == pytest.approx(nat.solvationHBEnergyHartree, abs=1e-10)
            assert cmp_res.solvationShiftEnergyHartree == pytest.approx(nat.solvationShiftEnergyHartree, abs=1e-10)
            np.testing.assert_allclose(cmp_res.GetDipole(), nat.GetDipole(), atol=1e-9)
            np.testing.assert_allclose(cmp_res.GetAtomicCharges(), nat.GetAtomicCharges(), atol=1e-8)


@pytest.mark.integration
def test_calc_solvation_descriptors_physics_permutation_and_error_recovery(library):
    mol = _make_mol('O', seed=17)
    original_bin = mol.ToBinary()

    desc_water = CalcSolvationDescriptors(mol, 'water', solvationModel='alpb', options=XTBOptions(xtbLibrary=library))
    assert mol.ToBinary() == original_bin
    assert desc_water['solvent'] == 'water'
    assert desc_water['solvation_model'] == 'alpb'
    assert desc_water['delta_g_solv_ha'] == pytest.approx(
        desc_water['e_solv_ha'] - desc_water['e_gas_ha'], abs=1e-12
    )
    assert desc_water['delta_g_solv_ha'] < 0.0  # water in water is thermodynamically stabilized
    assert desc_water['delta_g_solv_kcalmol'] < -2.0
    assert desc_water['delta_dipole_debye'] > 0.0  # polar solvent enhances dipole moment
    assert desc_water['solute_polarization_energy_ha'] >= -1e-6
    assert desc_water['charge_shift_rms_e'] > 0.0
    assert len(desc_water['atoms']) == 3
    json.dumps(desc_water, allow_nan=False)

    # Atom-permutation invariance
    perm = [2, 0, 1]
    reordered = Chem.RenumberAtoms(mol, perm)
    desc_perm = CalcSolvationDescriptors(reordered, 'water', solvationModel='alpb', options=XTBOptions(xtbLibrary=library))
    assert desc_perm['delta_g_solv_ha'] == pytest.approx(desc_water['delta_g_solv_ha'], abs=1e-8)
    for new_idx, old_idx in enumerate(perm):
        assert desc_perm['atoms'][new_idx]['solv_delta_q_e'] == pytest.approx(
            desc_water['atoms'][old_idx]['solv_delta_q_e'], abs=1e-8
        )

    # GFN0 solvation: CM5 shifts must be None with explicit missing reasons
    desc_gfn0 = CalcSolvationDescriptors(mol, 'water', solvationModel='alpb', options=XTBOptions(method='GFN0-xTB', xtbLibrary=library))
    assert desc_gfn0['atoms'][0]['cm5_q_gas_e'] is None
    assert desc_gfn0['atoms'][0]['missing']['solv_delta_cm5_q_e'] == 'cm5_undefined_for_gfn0_eeq'
    assert desc_gfn0['missing']['cm5_charges'] == 'cm5_undefined_for_gfn0_eeq'
    assert desc_gfn0['g_born_ha'] is None
    assert desc_gfn0['missing']['g_born_ha'] == 'not_separately_returned_by_gfn0'
    json.dumps(desc_gfn0, allow_nan=False)

    # Unknown solvent fails cleanly and does not contaminate subsequent gas-phase calls
    gas_before = XTBCalculator(xtbLibrary=library).Calculate(mol)
    with pytest.raises(XTBEnvironmentError):
        XTBCalculator(solvent='unobtainium', xtbLibrary=library).Calculate(mol)
    gas_after = XTBCalculator(xtbLibrary=library).Calculate(mol)
    assert gas_after.totalEnergyHartree == pytest.approx(gas_before.totalEnergyHartree, abs=1e-12)


@pytest.mark.integration
def test_vibrational_thermochemistry_optimized_water_and_invariance(library):
    mol = _make_mol('O', seed=17)
    opt = OptimizeXTB(mol, method='GFN2-xTB', gradientTolerance=1e-5, xtbLibrary=library)
    assert opt.success
    opt_mol = opt.GetMolecule()
    original_bin = opt_mol.ToBinary()

    thermo = CalcVibrationalThermochemistry(
        opt_mol, temperatureK=298.15, pressureAtm=1.0, symmetryNumber=2, options=XTBOptions(xtbLibrary=library)
    )
    assert opt_mol.ToBinary() == original_bin
    assert thermo['is_stationary_point'] is True
    assert thermo['warnings'] == []
    assert thermo['is_atom'] is False and thermo['is_linear'] is False
    assert thermo['n_vibrational_modes'] == 3
    assert thermo['n_imaginary'] == 0
    freqs = thermo['frequencies_cm1']
    assert len(freqs) == 3
    # GFN2-xTB water harmonic frequencies: bend ~1500 cm^-1, stretches ~3600-3800 cm^-1
    assert 1200.0 < freqs[0] < 1900.0
    assert 3200.0 < freqs[1] < 4100.0
    assert 3200.0 < freqs[2] < 4100.0
    assert all(intensity >= 0.0 for intensity in thermo['ir_intensities_km_mol'])
    assert max(thermo['ir_intensities_km_mol']) > 1.0

    # Exact thermodynamic identities
    assert thermo['zpe_ha'] > 0.015
    assert thermo['thermal_enthalpy_correction_ha'] > thermo['thermal_energy_correction_ha'] > thermo['zpe_ha']
    assert thermo['entropy_cal_mol_k'] > 30.0
    assert thermo['total_gibbs_free_energy_ha'] < thermo['total_enthalpy_ha']
    assert thermo['gibbs_correction_ha'] == pytest.approx(
        thermo['thermal_enthalpy_correction_ha'] - thermo['ts_ha'], abs=1e-12
    )
    assert thermo['total_gibbs_free_energy_ha'] == pytest.approx(
        thermo['electronic_energy_ha'] + thermo['gibbs_correction_ha'], abs=1e-12
    )
    json.dumps(thermo, allow_nan=False)

    # Pure RRHO mode
    thermo_rrho = CalcVibrationalThermochemistry(
        opt_mol, rotorCutoffCm1=0.0, symmetryNumber=2, options=XTBOptions(xtbLibrary=library)
    )
    assert thermo_rrho['thermo_model'] == 'rrho'
    assert thermo_rrho['zpe_ha'] == pytest.approx(thermo['zpe_ha'], abs=1e-12)

    # Rigid motion + atom permutation invariance
    moved = _rotate_translate_permute(opt_mol, [2, 0, 1])
    thermo_moved = CalcVibrationalThermochemistry(
        moved, temperatureK=298.15, pressureAtm=1.0, symmetryNumber=2, options=XTBOptions(xtbLibrary=library)
    )
    np.testing.assert_allclose(thermo_moved['frequencies_cm1'], thermo['frequencies_cm1'], atol=0.5)
    assert thermo_moved['zpe_ha'] == pytest.approx(thermo['zpe_ha'], abs=1e-6)
    assert thermo_moved['entropy_j_mol_k'] == pytest.approx(thermo['entropy_j_mol_k'], abs=1e-2)
    assert thermo_moved['total_gibbs_free_energy_ha'] == pytest.approx(thermo['total_gibbs_free_energy_ha'], abs=1e-6)


@pytest.mark.integration
def test_vibrational_thermochemistry_linear_atom_and_nonstationary_flags(library):
    # Linear molecule: CO2
    co2 = Chem.MolFromSmiles('O=C=O')
    conf = Chem.Conformer(3)
    conf.Set3D(True)
    conf.SetAtomPosition(0, (-1.16, 0.0, 0.0))
    conf.SetAtomPosition(1, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(2, (1.16, 0.0, 0.0))
    co2.AddConformer(conf, assignId=True)
    opt_co2 = OptimizeXTB(co2, method='GFN2-xTB', gradientTolerance=1e-5, xtbLibrary=library)
    assert opt_co2.success
    th_co2 = CalcVibrationalThermochemistry(opt_co2.GetMolecule(), symmetryNumber=2, options=XTBOptions(xtbLibrary=library))
    assert th_co2['is_linear'] is True
    assert th_co2['n_vibrational_modes'] == 4  # 3N - 5
    assert len(th_co2['rotational_constants_cm1']) == 1
    assert th_co2['n_imaginary'] == 0
    assert th_co2['is_stationary_point'] is True

    # Monoatomic system: Ar
    ar = Chem.MolFromSmiles('[Ar]')
    ar_conf = Chem.Conformer(1)
    ar_conf.Set3D(True)
    ar_conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    ar.AddConformer(ar_conf, assignId=True)
    th_ar = CalcVibrationalThermochemistry(ar, options=XTBOptions(xtbLibrary=library))
    assert th_ar['is_atom'] is True
    assert th_ar['n_vibrational_modes'] == 0
    assert th_ar['zpe_ha'] == 0.0
    assert th_ar['entropy_components_j_mol_k']['rotational'] == 0.0
    assert th_ar['entropy_components_j_mol_k']['vibrational'] == 0.0
    # Argon standard molar entropy at 298.15 K, 1 atm is ~154.8 J/(mol*K)
    assert 153.0 < th_ar['entropy_j_mol_k'] < 156.5

    # Strongly distorted non-stationary geometry (linear stretched H2O) must be flagged!
    distorted = _make_mol('O', seed=17)
    dconf = distorted.GetConformer()
    dconf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    dconf.SetAtomPosition(1, (1.10, 0.0, 0.0))
    dconf.SetAtomPosition(2, (-1.10, 0.0, 0.0))
    th_bad = CalcVibrationalThermochemistry(distorted, options=XTBOptions(xtbLibrary=library))
    assert th_bad['is_stationary_point'] is False
    assert 'non_stationary_geometry_large_gradient' in th_bad['warnings']
    assert th_bad['max_atom_gradient_ha_bohr'] > 1e-3
    assert th_bad['n_imaginary'] >= 1
    assert 'imaginary_vibrational_modes_present' in th_bad['warnings']
    json.dumps(th_bad, allow_nan=False)


def test_top_level_exports_catalog_and_serialization_roundtrip(tmp_path):
    import eled
    from eled import (
        DescriptorCatalog,
        DescriptorStatus,
        XTBDiskCache,
        XTBResults,
    )

    expected_exports = (
        'ESPSurfaceGrid',
        'GenerateESPSurfaceGrid',
        'CalcESPSurfaceDescriptors',
        'CalcOrbitalWeightedFukui',
        'CalcAtomicOrbitalEnergies',
        'CalcBondOverlapDescriptors',
        'CalcLocalizedOrbitals',
        'CalcEnergyDecomposition',
        'CalcMultipoleDescriptors',
        'CalcNuclearESP',
        'CalcPolarizabilityDescriptors',
        'CalcSolvationDescriptors',
        'CalcVibrationalThermochemistry',
    )
    for name in expected_exports:
        assert name in eled.__all__
        assert hasattr(eled, name)
    for removed in (
        'CalcFiniteFieldPolarizability',
        'CalcMayerBondDescriptors',
        'CalcEHTDescriptors',
        'CalcOrbitalDescriptors',
        'DescriptorSensitivity',
    ):
        assert removed not in eled.__all__
        assert not hasattr(eled, removed)
    assert 'eas' not in eled.__all__

    catalog = DescriptorCatalog()
    for section in ('esp_surface', 'localized_orbitals', 'properties', 'solvation', 'thermochemistry'):
        assert section in catalog and len(catalog[section]) >= 2
    status = DescriptorStatus()
    assert set(status) == {'recommended', 'diagnostic', 'exploratory'}

    synthetic = XTBResults(
        method='GFN2-xTB',
        totalEnergyHartree=-5.07,
        atomicCharges=(-0.5, 0.25, 0.25),
        solvent='water',
        solvationModel='alpb',
        solvationGrid=230,
        energyDecompositionFlags=1 | 2 | 4 | 8 | 16 | 32 | 64 | 256 | 512 | 1024,
        electronicEnergyHartree=-5.5,
        repulsionEnergyHartree=0.43,
        electrostaticEnergyHartree=-0.1,
        anisotropicElectrostaticEnergyHartree=-0.02,
        anisotropicXCEnergyHartree=-0.01,
        dispersionEnergyHartree=-0.005,
        solvationEnergyHartree=-0.012,
        solvationBornEnergyHartree=-0.010,
        solvationSASAEnergyHartree=0.002,
        solvationHBEnergyHartree=-0.003,
        solvationShiftEnergyHartree=-0.001,
        fermiLevelAlphaHartree=-0.25,
        fermiLevelBetaHartree=-0.25,
        fermiLevelAlphaEV=-6.8028,
        fermiLevelBetaEV=-6.8028,
        molecularPolarizabilityAU=8.5,
        molecularC6AU=45.0,
        molecularC8AU=1100.0,
        atomicPolarizabilitiesAU=(5.5, 1.5, 1.5),
        atomicC6AU=(15.0, 2.5, 2.5),
        atomicDipolesAU=((0.1, 0.0, 0.0), (-0.05, 0.02, 0.0), (-0.05, -0.02, 0.0)),
        atomicQuadrupolesAU=((0.1, 0.0, -0.05, 0.0, 0.0, -0.05),) * 3,
        molecularQuadrupoleAU=(0.3, 0.0, -0.15, 0.0, 0.0, -0.15),
        nuclearESPAU=(-22.3, -1.1, -1.1),
        metadata={'calculation_status': 'success', 'atom_count': 3},
    )
    restored = XTBResults.FromDict(json.loads(json.dumps(synthetic.ToDict(), allow_nan=False)))
    assert restored == synthetic
    assert isinstance(restored.atomicDipolesAU[0], tuple)
    assert isinstance(restored.atomicQuadrupolesAU[0], tuple)
    assert isinstance(restored.nuclearESPAU, tuple)

    disk_cache = XTBDiskCache(tmp_path / 'props_cache.sqlite')
    disk_cache.Set('solv_water', synthetic)
    loaded = disk_cache.Get('solv_water')
    assert loaded == synthetic
    assert loaded.solvent == 'water' and loaded.solvationModel == 'alpb' and loaded.solvationGrid == 230


@pytest.mark.integration
def test_cache_keys_disk_cache_and_spawned_workers_solvation_properties(library, tmp_path):
    from eled import RunXTBBatch, XTBCacheKey, XTBDiskCache

    mol = _make_mol('O', seed=17)
    k_gas = XTBCacheKey(mol, xtbLibrary=library)
    k_gas_default = XTBCacheKey(mol, solvent=None, solvationModel='alpb', solvationGrid=None, xtbLibrary=library)
    k_alpb = XTBCacheKey(mol, solvent='water', solvationModel='alpb', xtbLibrary=library)
    k_gbsa = XTBCacheKey(mol, solvent='water', solvationModel='gbsa', xtbLibrary=library)
    k_grid110 = XTBCacheKey(mol, solvent='water', solvationModel='alpb', solvationGrid=110, xtbLibrary=library)
    assert k_gas == k_gas_default
    assert len({k_gas, k_alpb, k_gbsa, k_grid110}) == 4

    cache = XTBDiskCache(tmp_path / 'batch_solv.sqlite')
    serial = RunXTBBatch([mol], solvent='water', solvationModel='alpb', xtbLibrary=library)
    parallel = RunXTBBatch(
        [mol, _make_mol('CO', seed=17)],
        solvent='water',
        solvationModel='alpb',
        xtbLibrary=library,
        maxWorkers=2,
        timeout=30,
        cache=cache,
    )
    assert len(cache) == 2
    assert serial[0].solvent == 'water' and parallel[0].solvent == 'water'
    assert parallel[0].totalEnergyHartree == pytest.approx(serial[0].totalEnergyHartree, abs=1e-10)
    assert parallel[0].solvationEnergyHartree == pytest.approx(serial[0].solvationEnergyHartree, abs=1e-10)
    assert parallel[0].molecularPolarizabilityAU == pytest.approx(serial[0].molecularPolarizabilityAU, abs=1e-10)
    np.testing.assert_allclose(parallel[0].GetNuclearESP(), serial[0].GetNuclearESP(), atol=1e-10)
    np.testing.assert_allclose(parallel[0].GetAtomicDipoles(), serial[0].GetAtomicDipoles(), atol=1e-10)

    # Disk cache hit must return identical extended properties without invoking a worker.
    cached = RunXTBBatch(
        [mol],
        solvent='water',
        solvationModel='alpb',
        xtbLibrary=library,
        timeout=1e-6,
        cache=cache,
    )[0]
    assert cached == parallel[0]
