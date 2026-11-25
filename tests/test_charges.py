"""Charge-model semantics and CM5 conservation contracts."""
from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from eled import CalcCM5Charges, CalcReactivityDescriptors, DescriptorUnavailableError, XTBResults
from eled.charges import (CHARGE_MODEL_CM5, CHARGE_MODEL_EEQ, CHARGE_MODEL_SCC,
                          _calc_cm5_charges, _charge_model_for_method)
from eled.input import _input_metadata, _xyz_from_mol


def test_engine_charge_models_are_method_specific_and_cm5_is_distinct():
    assert _charge_model_for_method('GFN0-xTB') == CHARGE_MODEL_EEQ
    assert _charge_model_for_method('GFN1-xTB') == CHARGE_MODEL_SCC
    assert _charge_model_for_method('GFN2-xTB') == CHARGE_MODEL_SCC
    q = _calc_cm5_charges((6, 1), ((0., 0., 0.), (1.09, 0., 0.)), (0., 0.))
    assert q[0] != 0. and q[1] != 0.
    assert sum(q) == pytest.approx(0.)


def test_cm5_accessor_requires_explicit_derived_values():
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    conf = Chem.Conformer(3)
    conf.Set3D(True)
    for i, xyz in enumerate(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.))):
        conf.SetAtomPosition(i, xyz)
    mol.AddConformer(conf)
    metadata = _input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))
    result = XTBResults('GFN2-xTB', None, -1., 1., 3, 8, 0,
                        orbitalEnergiesEV=(-2., 1., 2.), orbitalOccupancies=(2., 2., 0.),
                        atomicCharges=(-.6, .2, .4), chargeModel=CHARGE_MODEL_SCC,
                        metadata=metadata)
    with pytest.raises(DescriptorUnavailableError, match='unavailable'):
        CalcCM5Charges(replace(result, cm5Charges=()))
    cm5 = _calc_cm5_charges(metadata['atomic_numbers'], metadata['coordinates_angstrom'], result.atomicCharges)
    enriched = replace(result, cm5Charges=cm5)
    np.testing.assert_allclose(CalcCM5Charges(enriched), cm5)
    assert sum(CalcCM5Charges(enriched)) == pytest.approx(0.)


def test_reactivity_cm5_scheme_is_explicit_and_conserves_fukui_sums():
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    conf = Chem.Conformer(3)
    conf.Set3D(True)
    for i, xyz in enumerate(((0., 0., 0.), (1., 0., 0.), (0., 1., 0.))):
        conf.SetAtomPosition(i, xyz)
    mol.AddConformer(conf)
    rows = []
    for charge, energy, values, occupations, spin in [
        (0, -10., (-.6, .3, .3), (2., 2., 2., 2., 0.), 0),
        (1, -9.6, (-.3, .65, .65), (2., 2., 2., 1., 0.), 1),
        (-1, -10.1, (-.8, -.1, -.1), (2., 2., 2., 2., 1.), 1),
    ]:
        metadata = _input_metadata(mol, -1, charge, spin, _xyz_from_mol(mol, -1))
        metadata.update(engine={'sha256': 'engine'}, parameter_file={'sha256': 'params'},
                        calculation_settings={'electronic_temperature': 0.,
                                              'effective_electronic_temperature': 0.},
                        calculation_status='success')
        q = _calc_cm5_charges(metadata['atomic_numbers'], metadata['coordinates_angstrom'], values)
        rows.append(XTBResults('GFN2-xTB', None, energy, None, 5, int(sum(occupations)), spin,
                               orbitalEnergiesEV=(-10., -8., -6., -4., 1.),
                               orbitalOccupancies=occupations, atomicCharges=values,
                               chargeModel=CHARGE_MODEL_SCC, cm5Charges=q, metadata=metadata))
    report = CalcReactivityDescriptors(mol, *rows, populationScheme='cm5')
    assert report['analysis_settings']['population_scheme'] == 'cm5'
    assert report['analysis_settings']['charge_model'] == CHARGE_MODEL_CM5
    assert sum(atom['f_plus'] for atom in report['atoms']) == pytest.approx(1.)
    assert sum(atom['f_minus'] for atom in report['atoms']) == pytest.approx(1.)


@pytest.mark.parametrize('numbers', [(5, 14), (6, 17), (1, 9), (6, 8), (7, 8)])
def test_cm5_element_order_invariance(numbers):
    q_ab = _calc_cm5_charges(numbers, ((0.0, 0.0, 0.0), (2.0, 0.0, 0.0)), (0.2, -0.2))
    q_ba = _calc_cm5_charges(numbers[::-1], ((2.0, 0.0, 0.0), (0.0, 0.0, 0.0)), (-0.2, 0.2))
    np.testing.assert_allclose(q_ab, (q_ba[1], q_ba[0]))
    assert sum(q_ab) == pytest.approx(0.0)


def test_cm5_same_element_zero_distance_raises_value_error():
    with pytest.raises(ValueError, match='distinct atom coordinates'):
        _calc_cm5_charges((6, 6), ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0)), (0.0, 0.0))


@pytest.mark.parametrize('method', ['GFN1-xTB', 'GFN2-xTB'])
def test_cm5_fukui_requires_one_charge_per_atom(method):
    from test_reactivity import states
    mol, results = states()
    results = [replace(r, method=method, cm5Charges=(r.metadata['charge'],)) for r in results]
    with pytest.raises(ValueError, match='atom count'):
        CalcReactivityDescriptors(mol, *results, populationScheme='cm5')


@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
def test_report_with_missing_overlap_and_cm5(method):
    from eled import CalcDescriptorReport
    from test_reactivity import states
    mol, results = states()
    result = replace(results[0], method=method, chargeModel=_charge_model_for_method(method),
                     overlapMatrix=())
    report = CalcDescriptorReport(mol, result)
    assert report['values']['xtb_energy_ha'] == result.totalEnergyHartree
    assert report['missing']['populations'] == 'basis_not_retained'
    for suffix in ('sum', 'min', 'max', 'mean', 'sd', 'abs'):
        key = f'xtb_cm5_q_{suffix}_e'
        assert report['values'][key] is None
        assert report['missing'][key] == ('cm5_unsupported_method' if method == 'GFN0-xTB'
                                          else 'cm5_not_reported')
    if method == 'GFN0-xTB':
        invalid = replace(result, cm5Charges=result.atomicCharges)
        with pytest.raises(ValueError, match='unavailable'):
            CalcCM5Charges(invalid)
        with pytest.raises(ValueError, match='unavailable'):
            CalcDescriptorReport(mol, invalid)


@pytest.mark.parametrize('elements,radii,pair_parameter', [
    ((1, 8), (0.32, 0.64), 0.1671),
    ((5, 14), (0.84, 1.14), -0.1030 - (-0.0790)),
])
def test_cm5_reference_pair_values(elements, radii, pair_parameter):
    from math import exp
    distance = 1.7
    coords = ((0., 0., 0.), (distance, 0., 0.))
    correction = pair_parameter * exp(-2.474 * (distance - sum(radii)))
    expected = (0.2 + correction, -0.2 - correction)
    for order in (slice(None), slice(None, None, -1)):
        q = _calc_cm5_charges(elements[order], coords[order], (0.2, -0.2)[order])
        assert q == pytest.approx(expected[order], abs=1e-12)
        assert sum(q) == pytest.approx(0., abs=1e-12)


@pytest.mark.parametrize('values,message', [((0., 0.), 'expected 3, got 2'),
                                           ((float('nan'), 0., 0.), 'finite')])
def test_corrupt_cm5_rejected_across_interfaces(monkeypatch, values, message):
    from eled import (AssessDescriptorQuality, CalcDescriptorReport, CalcXTBDescriptors,
                      CalcXTBAtomDescriptors, CalcXTBFragmentDescriptors, RunDescriptorBatch)
    from test_reactivity import states
    mol, rows = states()
    rows = [replace(r, cm5Charges=r.atomicCharges) for r in rows]
    bad = replace(rows[0], cm5Charges=values)
    for read in (CalcCM5Charges, CalcXTBDescriptors, CalcXTBAtomDescriptors,
                 AssessDescriptorQuality, lambda r: CalcXTBFragmentDescriptors(r, {'O': [0]}),
                 lambda r: CalcDescriptorReport(mol, r),
                 lambda r: CalcReactivityDescriptors(mol, r, *rows[1:], populationScheme='cm5')):
        with pytest.raises(ValueError, match=message):
            read(bad)
    monkeypatch.setattr('eled.batch.RunXTBBatch', lambda *a, **kw: [bad, rows[0]])
    failed, valid = RunDescriptorBatch([mol, mol])
    assert message in failed['error']['message']
    assert valid['error'] is None
