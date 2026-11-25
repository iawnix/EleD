"""Charge-fit inverse problems and sampled-ESP backend/cache contracts."""
from dataclasses import replace
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (CalcESP, DescriptorUnavailableError, FitESPCharges, FitRESPCharges,
                  GenerateESPGrid, XTBCalculator, XTBCacheKey, XTBDiskCache, XTBOptions, XTBResults)
from eled.esp import ANGSTROM_TO_BOHR


def synthetic():
    coords = np.array([[0., 0., 0.], [.9, .5, 0.], [-.9, .5, 0.]])
    points = np.random.default_rng(42).normal(size=(100, 3))
    points *= 3 / np.linalg.norm(points, axis=1)[:, None]
    q = np.array([-.8, .4, .4])
    design = 1/(np.linalg.norm(points[:, None]-coords, axis=2)*ANGSTROM_TO_BOHR)
    result = XTBResults('GFN2-xTB', None, -1., None, 6, 8, 0,
        espPointsAngstrom=tuple(map(tuple, points)), espValuesAU=tuple(design@q),
        metadata={'coordinates_angstrom': coords.tolist(), 'atomic_numbers': [8, 1, 1],
                  'charge': 0, 'calculation_status': 'success'})
    return result, q, design


def test_esp_recovers_known_charges_constraints_and_json():
    result, q, _ = synthetic()
    for kwargs in ({}, {'equivalenceGroups': [(1, 2)]}, {'fixedCharges': {0: -.8}},
                   {'fixedCharges': dict(enumerate(q))}, {'weights': np.linspace(.2, 2, 100)}):
        fit = FitESPCharges(result, **kwargs)
        np.testing.assert_allclose(fit.charges, q, atol=1e-12)
        np.testing.assert_allclose(fit.GetCharges(), q, atol=1e-12)
        owned = fit.GetCharges(); owned[:] = 999
        np.testing.assert_allclose(fit.GetCharges(), q, atol=1e-12)
        assert fit.constraintError < 1e-12 and fit.weightedRMSAU < 1e-12
        assert json.loads(json.dumps(fit.ToDict(), allow_nan=False))['charge_model'] == 'esp'
    charged = replace(result, espValuesAU=tuple(np.array(result.espValuesAU)+synthetic()[2][:, 0]),
                      metadata={**result.metadata, 'charge': 1})
    np.testing.assert_allclose(FitESPCharges(charged).charges, q+[1, 0, 0], atol=1e-12)


def test_resp_stationarity_shrinkage_and_constraint_conservation():
    result, original, design = synthetic()
    fit = FitRESPCharges(result, restraintStrength=.03, tolerance=1e-12)
    q = np.array(fit.charges)
    assert abs(q[0]) < abs(original[0]) and abs(sum(q)) < 1e-12
    # Independent derivative of the documented objective, projected onto sum(q)=0.
    gradient = design.T@(design@q-np.array(result.espValuesAU))
    gradient[0] += .03*q[0]/np.sqrt(q[0]**2+.1**2)
    assert np.max(np.abs(gradient-gradient.mean())) < 1e-10
    constrained = FitRESPCharges(result, equivalenceGroups=[(1, 2)], fixedCharges={0: -.7})
    np.testing.assert_allclose(constrained.charges, [-.7, .35, .35], atol=1e-12)
    with pytest.raises(RuntimeError, match='converge'):
        FitRESPCharges(result, restraintStrength=.1, maxIterations=1, tolerance=1e-14)


def test_rank_deficiency_and_resp_identifiability():
    result, _, _ = synthetic()
    result = replace(result, espPointsAngstrom=((0., 3., 0.),)*4, espValuesAU=(.1,)*4)
    with pytest.raises(ValueError, match='rank deficient'):
        FitESPCharges(result)
    fit = FitRESPCharges(result, restrainHydrogens=True)
    assert fit.dataRank < fit.freeParameters and fit.conditionNumber is None
    assert sum(fit.charges) == pytest.approx(0, abs=1e-12)
    with pytest.raises(ValueError, match='determine'):
        FitRESPCharges(result)  # equivalent H columns, neither restrained


@pytest.mark.parametrize('kwargs', [
    {'fixedCharges': {0: 1., 1: 1., 2: 1.}},
    {'equivalenceGroups': [(0, 1)], 'fixedCharges': {0: -.8, 1: .4}},
    {'equivalenceGroups': [(0, 0)]}, {'fixedCharges': {True: 0}},
    {'fixedCharges': {4: 0}}, {'fixedCharges': {0: np.nan}},
    {'weights': [-1.]*100}, {'weights': [0.]*100}, {'weights': [1.]}, {'rcond': 0.},
])
def test_bad_fit_constraints_and_weights_fail(kwargs):
    with pytest.raises(ValueError):
        FitESPCharges(synthetic()[0], **kwargs)


def test_fit_is_invariant_to_rigid_motion_atom_and_point_permutation():
    result, q, _ = synthetic()
    rotation = np.array([[.6, -.8, 0], [.8, .6, 0], [0, 0, 1]])
    permutation = [2, 0, 1]
    coords = np.array(result.metadata['coordinates_angstrom'])[permutation]@rotation.T+1.2
    points = np.array(result.espPointsAngstrom)[::-1]@rotation.T+1.2
    changed = replace(result, espPointsAngstrom=tuple(map(tuple, points)),
                      espValuesAU=result.espValuesAU[::-1], metadata={**result.metadata,
                      'coordinates_angstrom': coords.tolist(), 'atomic_numbers': [1, 8, 1]})
    np.testing.assert_allclose(FitESPCharges(changed).charges, q[permutation], atol=1e-12)
    np.testing.assert_allclose(FitRESPCharges(changed).charges,
                               np.array(FitRESPCharges(result).charges)[permutation], atol=1e-12)


@pytest.mark.parametrize('points', [[[0, 0]], [[0, 0, np.nan]], [[0, 0, 1j]], [[True, False, True]],
                                   ['1', '2', '3'], np.empty((2, 0))])
def test_invalid_sampling_options_fail_before_calculation(points):
    with pytest.raises(ValueError, match='espPointsAngstrom'):
        XTBOptions(espPointsAngstrom=points)


def test_sampling_options_are_owned_and_roundtrip():
    points = np.array([[2., 3., 4.]])
    options = XTBOptions(espPointsAngstrom=points)
    points[:] = 0
    assert options.espPointsAngstrom == ((2., 3., 4.),)
    assert json.loads(json.dumps(options.ToDict())) == options.ToDict()
    assert XTBOptions.FromDict(json.loads(json.dumps(options.ToDict()))) == options
    result = synthetic()[0]
    with pytest.raises(DescriptorUnavailableError):
        CalcESP(replace(result, espPointsAngstrom=(), espValuesAU=()))
    for bad in (replace(result, espValuesAU=(np.nan,)*100),
                replace(result, espValuesAU=(0.,)), replace(result, espPointsAngstrom=())):
        with pytest.raises(ValueError): CalcESP(bad)


def molecule(smiles='O'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    return mol


def test_grid_exposed_shells_and_reproducibility():
    mol = molecule()
    points = GenerateESPGrid(mol, scales=(1.4,), pointsPerAtom=50)
    assert points == GenerateESPGrid(mol, scales=(1.4,), pointsPerAtom=50)
    coordinates = mol.GetConformer().GetPositions()
    radii = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    distances = np.linalg.norm(np.array(points)[:, None]-coordinates, axis=2)
    assert np.all(distances >= 1.4*radii-1e-10)
    np.testing.assert_allclose(np.min(abs(distances-1.4*radii), axis=1), 0, atol=1e-10)
    for kwargs in ({'scales': (1.,)}, {'scales': (1.4, 1.4)}, {'pointsPerAtom': True}):
        with pytest.raises(ValueError): GenerateESPGrid(mol, **kwargs)


@pytest.mark.integration
@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['O', 'S', '[NH4+]', '[OH]'])
def test_native_compiled_sample_and_fit_parity(method, smiles):
    from eled import rdXTBTools
    mol = molecule(smiles)
    points = GenerateESPGrid(mol, pointsPerAtom=30)
    native = XTBCalculator(method=method).Calculate(mol, espPointsAngstrom=points)
    ok, compiled = rdXTBTools.RunMol(mol, method=method, espPointsAngstrom=points)
    assert ok, compiled.errorMessage
    assert compiled.espSampled
    np.testing.assert_allclose(CalcESP(native), CalcESP(compiled), atol=1e-12, rtol=0)
    np.testing.assert_array_equal(compiled.GetESPPoints(), points)
    view = compiled.GetESPValues(); view[:] = 0
    np.testing.assert_allclose(compiled.GetESPValues(), CalcESP(native), atol=1e-12, rtol=0)
    for function in [FitESPCharges, FitRESPCharges]:
        a, b = function(native), function(compiled)
        np.testing.assert_allclose(a.charges, b.charges, atol=1e-10, rtol=0)
        assert sum(a.charges) == pytest.approx(Chem.GetFormalCharge(mol), abs=1e-10)
    assert native.totalEnergyHartree == pytest.approx(XTBCalculator(method=method).Calculate(mol).totalEnergyHartree, abs=1e-12)


@pytest.mark.integration
def test_esp_cache_keys_and_disk_roundtrip(tmp_path):
    mol = molecule()
    points = GenerateESPGrid(mol, pointsPerAtom=10)
    keys = [XTBCacheKey(mol, espPointsAngstrom=p) for p in [(), points, points[::-1]]]
    assert len(set(keys)) == 3
    result = XTBCalculator().Calculate(mol, espPointsAngstrom=points)
    cache = XTBDiskCache(tmp_path/'esp.sqlite')
    cache.Set(keys[1], result)
    restored = XTBDiskCache(tmp_path/'esp.sqlite').Get(keys[1])
    assert restored == result
    assert FitESPCharges(restored).charges == FitESPCharges(result).charges


def test_esp_capability_rejects_unverified_extension(monkeypatch):
    from eled.backends.native import _NativeXTB
    monkeypatch.setattr('eled.runtime._native_esp_build', lambda path: False)
    dummy = object.__new__(_NativeXTB); dummy.library = '/no-library'
    with pytest.raises(DescriptorUnavailableError, match='verified ESP'):
        dummy.bind_esp_api()


@pytest.mark.parametrize('read', [CalcESP, FitESPCharges, FitRESPCharges])
def test_sample_readers_distinguish_missing_and_malformed_data(read):
    result = synthetic()[0]
    with pytest.raises(DescriptorUnavailableError):
        read(replace(result, espPointsAngstrom=None, espValuesAU=None))
    for points, values in [((), 1.), ((), (1.,)), (((2., 0., 0.),), None)]:
        with pytest.raises(ValueError) as error:
            read(replace(result, espPointsAngstrom=points, espValuesAU=values))
        assert not isinstance(error.value, DescriptorUnavailableError)


@pytest.mark.parametrize('fit', [FitESPCharges, FitRESPCharges])
def test_fit_requires_consistent_source_atoms(fit):
    result = synthetic()[0]
    for fields in ({'atom_count': 2}, {'atomic_numbers': None}):
        with pytest.raises(ValueError, match='atom'):
            fit(replace(result, metadata={**result.metadata, **fields}))
