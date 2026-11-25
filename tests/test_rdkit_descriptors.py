"""RDKit naming, single-call wrappers, PropertyFunctor registry and result contracts."""
from dataclasses import replace
import gc
import math

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors

import eled
from eled import (
    AnnotateMol,
    CalcEHTDescriptors,
    CalcMolXTBDescriptors,
    CalcOrbitalDescriptors,
    CalcXTBDescriptors,
    RegisterRDKitProperties,
    XTBCache,
    XTBCalculator,
    XTBNativeNotFoundError,
    XTBResult,
    XTBResults,
)
from eled.input import _input_metadata, _xyz_from_mol
from eled.runtime import _extended_libxtb_build, _find_xtb_library


def _water(conf_id=0):
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    mol.GetConformer().SetId(conf_id)
    return mol


def _synthetic_water_result(mol):
    return XTBResults(
        'GFN2-xTB',
        'test',
        -5.07,
        3.0,
        3,
        8,
        0,
        orbitalEnergiesEV=(-12.0, -5.0, 2.0),
        orbitalOccupancies=(2.0, 6.0, 0.0),
        atomicCharges=(-0.6, 0.3, 0.3),
        cm5Charges=(-0.5, 0.25, 0.25),
        dipoleAtomicUnits=(0.0, 0.0, 0.8),
        wibergBondOrders=((0.0, 0.95, 0.95), (0.95, 0.0, 0.01), (0.95, 0.01, 0.0)),
        orbitalCoefficients=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        overlapMatrix=((1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)),
        aoAtomMap=(0, 1, 2),
        aoShellMap=(0, 1, 2),
        shellAtomMap=(0, 1, 2),
        shellAngularMomenta=(0, 0, 0),
        aoPDirections=((0.0, 0.0, 0.0),) * 3,
        gradientHartreePerBohr=((0.01, 0.0, 0.0), (-0.005, 0.0, 0.0), (-0.005, 0.0, 0.0)),
        espPointsAngstrom=((2.0, 0.0, 0.0),),
        espValuesAU=(-0.05,),
        metadata=_input_metadata(mol, mol.GetConformer().GetId(), 0, 0, _xyz_from_mol(mol, -1)),
    )


def test_top_level_namespace_exports_rdkit_helpers_and_excludes_eas_dataset_tools():
    import eled.eas as eas_module
    assert XTBResult is XTBResults
    assert callable(CalcMolXTBDescriptors)
    assert callable(RegisterRDKitProperties)
    assert callable(CalcEHTDescriptors)
    assert callable(CalcOrbitalDescriptors)
    for eas_name in (
        'NormalizeEASRecords',
        'AssignEASSplits',
        'AssignEASFolds',
        'PrepareEASDataset',
        'LoadFrozenProtocol',
        'FrozenProtocolSHA256',
        'ClusterKeys',
    ):
        assert eas_name not in eled.__all__
        assert not hasattr(eled, eas_name)
        assert hasattr(eas_module, eas_name)


def test_passing_mol_directly_to_result_descriptor_gives_actionable_type_error():
    mol = _water()
    with pytest.raises(TypeError, match='CalcMolXTBDescriptors'):
        CalcXTBDescriptors(mol)


def test_annotate_mol_private_props_follow_rdkit_underscore_convention():
    mol = _water()
    result = replace(_synthetic_water_result(mol), orbitalOccupancies=(2.0, 2.0, 0.0), numElectrons=4)
    private_mol = AnnotateMol(mol, result, privateProps=True)
    assert private_mol.GetAtomWithIdx(0).GetDoubleProp('_xtb_partial_charge') == pytest.approx(-0.6)
    assert private_mol.GetAtomWithIdx(0).GetDoubleProp('_xtb_wbo_sum') == pytest.approx(1.9)
    assert private_mol.GetBondWithIdx(0).GetDoubleProp('_xtb_wbo') == pytest.approx(0.95)
    assert private_mol.GetDoubleProp('_xtb_gap_ev') == pytest.approx(7.0)
    assert private_mol.HasProp('_xtb_metadata_json')
    assert '_xtb_gap_ev' not in private_mol.GetPropsAsDict()
    assert '_xtb_gap_ev' in private_mol.GetPropsAsDict(includePrivate=True)
    assert '_xtb_partial_charge' not in private_mol.GetAtomWithIdx(0).GetPropsAsDict(includePrivate=False)
    assert '_xtb_partial_charge' in private_mol.GetAtomWithIdx(0).GetPropsAsDict(includePrivate=True)


def test_python_xtb_results_getters_return_owned_numpy_copies_and_guard_failures():
    mol = _water(conf_id=4)
    result = _synthetic_water_result(mol)
    assert result.success and result.conformerId == 4 and result.numAtoms == 3 and result.charge == 0
    assert result.overlapRetained and result.espSampled
    for getter in (
        result.GetAtomicCharges,
        result.GetCM5Charges,
        result.GetOrbitalEnergies,
        result.GetOrbitalOccupancies,
        result.GetDipole,
        result.GetGradient,
        result.GetWibergBondOrders,
        result.GetOrbitalCoefficients,
        result.GetOverlapMatrix,
        result.GetAOAtomMap,
        result.GetAOShellMap,
        result.GetShellAtomMap,
        result.GetShellAngularMomenta,
        result.GetAOPDirections,
        result.GetESPPoints,
        result.GetESPValues,
    ):
        arr = getter()
        assert isinstance(arr, np.ndarray) and arr.flags.owndata
        snapshot = arr.copy()
        arr[...] = 0
        np.testing.assert_array_equal(getter(), snapshot)

    failed = replace(result, metadata={'calculation_status': 'error', 'error_message': 'SCF did not converge'})
    assert not failed.success
    with pytest.raises(RuntimeError, match='No numerical results'):
        failed.GetAtomicCharges()


def test_calc_mol_xtb_descriptors_missing_val_silent_and_cache(monkeypatch):
    mol = _water()
    sample = replace(_synthetic_water_result(mol), orbitalOccupancies=(2.0, 2.0, 0.0), numElectrons=4)
    calls = 0

    def fake_run(target, **kwargs):
        nonlocal calls
        _xyz_from_mol(target, kwargs.get('confId', -1))
        calls += 1
        return sample

    monkeypatch.setattr('eled.calculator._run_xtb', fake_run)
    monkeypatch.setattr(
        'eled.batch._prepare',
        lambda m, opts: (replace(opts, xtbLibrary='/fake/libxtb.so'), sample.metadata),
    )
    monkeypatch.setattr('eled.batch._key', lambda opts, meta: 'water-key')

    cache = XTBCache()
    first = CalcMolXTBDescriptors(mol, missingVal=float('nan'), cache=cache)
    second = CalcMolXTBDescriptors(mol, missingVal=float('nan'), cache=cache)
    assert calls == 1
    assert first['xtb_energy_ha'] == pytest.approx(-5.07)
    assert math.isnan(first['xtb_homo_m2_ev'])
    assert second['xtb_energy_ha'] == pytest.approx(-5.07)

    unembedded = Chem.MolFromSmiles('O')
    with pytest.raises(ValueError):
        CalcMolXTBDescriptors(unembedded, silent=False)
    suppressed = CalcMolXTBDescriptors(unembedded, missingVal=-999.0, silent=True)
    assert suppressed['xtb_energy_ha'] == -999.0
    assert suppressed['xtb_cm5_q_sd_e'] == -999.0


def test_register_rdkit_properties_deduplicates_and_survives_gc(monkeypatch):
    mol = _water()
    sample = replace(_synthetic_water_result(mol), orbitalOccupancies=(2.0, 2.0, 0.0), numElectrons=4)
    calls = 0

    def fake_run(target, **kwargs):
        nonlocal calls
        _xyz_from_mol(target, kwargs.get('confId', -1))
        calls += 1
        return sample

    monkeypatch.setattr('eled.calculator._run_xtb', fake_run)
    names = RegisterRDKitProperties(prefix='unit_')
    gc.collect()
    assert 'unit_xtb_energy_ha' in names
    assert 'unit_xtb_cm5_q_sd_e' in names
    assert 'unit_xtb_method' not in names

    props = rdMolDescriptors.Properties(list(names))
    values = dict(zip(names, props.ComputeProperties(mol)))
    assert calls == 1
    assert values['unit_xtb_energy_ha'] == pytest.approx(-5.07)
    assert values['unit_xtb_gap_ev'] == pytest.approx(7.0)
    assert math.isnan(values['unit_xtb_homo_m2_ev'])


@pytest.mark.integration
def test_real_rdkit_properties_and_single_call_wrapper_parity():
    try:
        library = _find_xtb_library()
    except XTBNativeNotFoundError:
        pytest.skip('extended libxtb unavailable')
    if not _extended_libxtb_build(library):
        pytest.skip('reviewed extended libxtb required')

    mol = _water()
    direct = CalcXTBDescriptors(XTBCalculator(xtbLibrary=library).Calculate(mol))
    wrapped = CalcMolXTBDescriptors(mol, xtbLibrary=library)
    assert wrapped == direct

    names = RegisterRDKitProperties(xtbLibrary=library, prefix='real_')
    computed = dict(zip(names, rdMolDescriptors.Properties(list(names)).ComputeProperties(mol)))
    for field in ('xtb_energy_ha', 'xtb_gap_ev', 'xtb_dip_debye', 'xtb_cm5_q_sd_e'):
        assert computed[f'real_{field}'] == pytest.approx(direct[field], abs=1e-12)
