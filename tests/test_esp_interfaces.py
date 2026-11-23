"""Charge-fit inverse problems and sampled-ESP backend/cache contracts."""
from dataclasses import replace
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (CalcESP, DescriptorUnavailableError, FitESPCharges, FitRESPCharges,
                  GenerateESPGrid, XTBCalculator, XTBCacheKey, XTBDiskCache, XTBOptions, XTBResult)
from eled.esp import ANGSTROM_TO_BOHR


def synthetic():
    coords = np.array([[0., 0., 0.], [.9, .5, 0.], [-.9, .5, 0.]])
    points = np.random.default_rng(42).normal(size=(100, 3))
    points *= 3 / np.linalg.norm(points, axis=1)[:, None]
    q = np.array([-.8, .4, .4])
    design = 1/(np.linalg.norm(points[:, None]-coords, axis=2)*ANGSTROM_TO_BOHR)
    result = XTBResult('GFN2-xTB', None, -1., None, 6, 8, 0,
        esp_points_angstrom=tuple(map(tuple, points)), esp_values_au=tuple(design@q),
        metadata={'coordinates_angstrom': coords.tolist(), 'atomic_numbers': [8, 1, 1],
                  'charge': 0, 'calculation_status': 'success'})
    return result, q, design


def test_esp_recovers_known_charges_constraints_and_json():
    result, q, _ = synthetic()
    for kwargs in ({}, {'equivalence_groups': [(1, 2)]}, {'fixed_charges': {0: -.8}},
                   {'fixed_charges': dict(enumerate(q))}, {'weights': np.linspace(.2, 2, 100)}):
        fit = FitESPCharges(result, **kwargs)
        np.testing.assert_allclose(fit.charges_e, q, atol=1e-12)
        assert fit.constraint_error_e < 1e-12 and fit.weighted_rms_au < 1e-12
        assert json.loads(json.dumps(fit.to_dict(), allow_nan=False))['charge_model'] == 'esp'
    charged = replace(result, esp_values_au=tuple(np.array(result.esp_values_au)+synthetic()[2][:, 0]),
                      metadata={**result.metadata, 'charge': 1})
    np.testing.assert_allclose(FitESPCharges(charged).charges_e, q+[1, 0, 0], atol=1e-12)


def test_resp_stationarity_shrinkage_and_constraint_conservation():
    result, original, design = synthetic()
    fit = FitRESPCharges(result, restraint_strength=.03, tolerance=1e-12)
    q = np.array(fit.charges_e)
    assert abs(q[0]) < abs(original[0]) and abs(sum(q)) < 1e-12
    # Independent derivative of the documented objective, projected onto sum(q)=0.
    gradient = design.T@(design@q-np.array(result.esp_values_au))
    gradient[0] += .03*q[0]/np.sqrt(q[0]**2+.1**2)
    assert np.max(np.abs(gradient-gradient.mean())) < 1e-10
    constrained = FitRESPCharges(result, equivalence_groups=[(1, 2)], fixed_charges={0: -.7})
    np.testing.assert_allclose(constrained.charges_e, [-.7, .35, .35], atol=1e-12)
    with pytest.raises(RuntimeError, match='converge'):
        FitRESPCharges(result, restraint_strength=.1, max_iterations=1, tolerance=1e-14)


def test_rank_deficiency_and_resp_identifiability():
    result, _, _ = synthetic()
    result = replace(result, esp_points_angstrom=((0., 3., 0.),)*4, esp_values_au=(.1,)*4)
    with pytest.raises(ValueError, match='rank deficient'):
        FitESPCharges(result)
    fit = FitRESPCharges(result, restrain_hydrogens=True)
    assert fit.data_rank < fit.free_parameters and fit.condition_number is None
    assert sum(fit.charges_e) == pytest.approx(0, abs=1e-12)
    with pytest.raises(ValueError, match='determine'):
        FitRESPCharges(result)  # equivalent H columns, neither restrained


@pytest.mark.parametrize('kwargs', [
    {'fixed_charges': {0: 1., 1: 1., 2: 1.}},
    {'equivalence_groups': [(0, 1)], 'fixed_charges': {0: -.8, 1: .4}},
    {'equivalence_groups': [(0, 0)]}, {'fixed_charges': {True: 0}},
    {'fixed_charges': {4: 0}}, {'fixed_charges': {0: np.nan}},
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
    points = np.array(result.esp_points_angstrom)[::-1]@rotation.T+1.2
    changed = replace(result, esp_points_angstrom=tuple(map(tuple, points)),
                      esp_values_au=result.esp_values_au[::-1], metadata={**result.metadata,
                      'coordinates_angstrom': coords.tolist(), 'atomic_numbers': [1, 8, 1]})
    np.testing.assert_allclose(FitESPCharges(changed).charges_e, q[permutation], atol=1e-12)
    np.testing.assert_allclose(FitRESPCharges(changed).charges_e,
                               np.array(FitRESPCharges(result).charges_e)[permutation], atol=1e-12)


@pytest.mark.parametrize('points', [[[0, 0]], [[0, 0, np.nan]], [[0, 0, 1j]], [[True, False, True]],
                                   ['1', '2', '3'], np.empty((2, 0))])
def test_invalid_sampling_options_fail_before_calculation(points):
    with pytest.raises(ValueError, match='esp_points'):
        XTBOptions(esp_points_angstrom=points)


def test_sampling_options_are_owned_and_roundtrip():
    points = np.array([[2., 3., 4.]])
    options = XTBOptions(esp_points_angstrom=points)
    points[:] = 0
    assert options.esp_points_angstrom == ((2., 3., 4.),)
    assert json.loads(json.dumps(options.to_dict())) == options.to_dict()
    assert XTBOptions(**json.loads(json.dumps(options.to_dict()))) == options
    result = synthetic()[0]
    with pytest.raises(DescriptorUnavailableError):
        CalcESP(replace(result, esp_points_angstrom=(), esp_values_au=()))
    for bad in (replace(result, esp_values_au=(np.nan,)*100),
                replace(result, esp_values_au=(0.,)), replace(result, esp_points_angstrom=())):
        with pytest.raises(ValueError): CalcESP(bad)


def molecule(smiles='O'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=42) == 0
    return mol


def test_grid_exposed_shells_and_reproducibility():
    mol = molecule()
    points = GenerateESPGrid(mol, scales=(1.4,), points_per_atom=50)
    assert points == GenerateESPGrid(mol, scales=(1.4,), points_per_atom=50)
    coordinates = mol.GetConformer().GetPositions()
    radii = np.array([Chem.GetPeriodicTable().GetRvdw(a.GetAtomicNum()) for a in mol.GetAtoms()])
    distances = np.linalg.norm(np.array(points)[:, None]-coordinates, axis=2)
    assert np.all(distances >= 1.4*radii-1e-10)
    np.testing.assert_allclose(np.min(abs(distances-1.4*radii), axis=1), 0, atol=1e-10)
    for kwargs in ({'scales': (1.,)}, {'scales': (1.4, 1.4)}, {'points_per_atom': True}):
        with pytest.raises(ValueError): GenerateESPGrid(mol, **kwargs)


@pytest.mark.integration
@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['O', 'S', '[NH4+]', '[OH]'])
def test_native_compiled_sample_and_fit_parity(method, smiles):
    from eled import rdXTBTools
    mol = molecule(smiles)
    points = GenerateESPGrid(mol, points_per_atom=30)
    native = XTBCalculator(method=method).calculate(mol, esp_points_angstrom=points)
    ok, compiled = rdXTBTools.RunMol(mol, method=method, espPointsAngstrom=points)
    assert ok, compiled.errorMessage
    assert compiled.espSampled
    np.testing.assert_allclose(CalcESP(native), CalcESP(compiled), atol=1e-12, rtol=0)
    np.testing.assert_array_equal(compiled.GetESPPoints(), points)
    view = compiled.GetESPValues(); view[:] = 0
    np.testing.assert_allclose(compiled.GetESPValues(), CalcESP(native), atol=1e-12, rtol=0)
    for function in [FitESPCharges, FitRESPCharges]:
        a, b = function(native), function(compiled)
        np.testing.assert_allclose(a.charges_e, b.charges_e, atol=1e-10, rtol=0)
        assert sum(a.charges_e) == pytest.approx(Chem.GetFormalCharge(mol), abs=1e-10)
    assert native.total_energy_hartree == pytest.approx(XTBCalculator(method=method).calculate(mol).total_energy_hartree, abs=1e-12)


@pytest.mark.integration
def test_esp_cache_keys_and_disk_roundtrip(tmp_path):
    mol = molecule()
    points = GenerateESPGrid(mol, points_per_atom=10)
    keys = [XTBCacheKey(mol, esp_points_angstrom=p) for p in [(), points, points[::-1]]]
    assert len(set(keys)) == 3
    result = XTBCalculator().calculate(mol, esp_points_angstrom=points)
    cache = XTBDiskCache(tmp_path/'esp.sqlite')
    cache.set(keys[1], result)
    restored = XTBDiskCache(tmp_path/'esp.sqlite').get(keys[1])
    assert restored == result
    assert FitESPCharges(restored).charges_e == FitESPCharges(result).charges_e


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
        read(replace(result, esp_points_angstrom=None, esp_values_au=None))
    for points, values in [((), 1.), ((), (1.,)), (((2., 0., 0.),), None)]:
        with pytest.raises(ValueError) as error:
            read(replace(result, esp_points_angstrom=points, esp_values_au=values))
        assert not isinstance(error.value, DescriptorUnavailableError)


@pytest.mark.parametrize('fit', [FitESPCharges, FitRESPCharges])
def test_fit_requires_consistent_source_atoms(fit):
    result = synthetic()[0]
    for fields in ({'atom_count': 2}, {'atomic_numbers': None}):
        with pytest.raises(ValueError, match='atom'):
            fit(replace(result, metadata={**result.metadata, **fields}))
