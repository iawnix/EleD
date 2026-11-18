"""Charge-model semantics and CM5 conservation contracts."""
from dataclasses import replace

import numpy as np
import pytest
from rdkit import Chem

from eled import CalcCM5Charges, CalcReactivityDescriptors, XTBResult
from eled.charges import (CHARGE_MODEL_CM5, CHARGE_MODEL_EEQ, CHARGE_MODEL_SCC,
                          calc_cm5_charges, charge_model_for_method)
from eled.input import _input_metadata, _xyz_from_mol


def test_engine_charge_models_are_method_specific_and_cm5_is_distinct():
    assert charge_model_for_method('GFN0-xTB') == CHARGE_MODEL_EEQ
    assert charge_model_for_method('GFN1-xTB') == CHARGE_MODEL_SCC
    assert charge_model_for_method('GFN2-xTB') == CHARGE_MODEL_SCC
    q = calc_cm5_charges((6, 1), ((0., 0., 0.), (1.09, 0., 0.)), (0., 0.))
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
    result = XTBResult('GFN2-xTB', None, -1., 1., 3, 8, 0,
                       orbital_energies_ev=(-2., 1., 2.), orbital_occupations=(2., 2., 0.),
                       partial_charges=(-.6, .2, .4), charge_model=CHARGE_MODEL_SCC,
                       metadata=metadata)
    with pytest.raises(ValueError, match='unavailable'):
        CalcCM5Charges(replace(result, cm5_charges=()))
    cm5 = calc_cm5_charges(metadata['atomic_numbers'], metadata['coordinates_angstrom'], result.partial_charges)
    enriched = replace(result, cm5_charges=cm5)
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
                        CalculationSettings={'electronic_temperature': 0.,
                                             'effective_electronic_temperature': 0.},
                        calculation_status='success')
        q = calc_cm5_charges(metadata['atomic_numbers'], metadata['coordinates_angstrom'], values)
        rows.append(XTBResult('GFN2-xTB', None, energy, None, 5, int(sum(occupations)), spin,
                              orbital_energies_ev=(-10., -8., -6., -4., 1.),
                              orbital_occupations=occupations, partial_charges=values,
                              charge_model=CHARGE_MODEL_SCC, cm5_charges=q, metadata=metadata))
    report = CalcReactivityDescriptors(mol, *rows, population_scheme='cm5')
    assert report['analysis_settings']['population_scheme'] == 'cm5'
    assert report['analysis_settings']['charge_model'] == CHARGE_MODEL_CM5
    assert sum(atom['f_plus'] for atom in report['atoms']) == pytest.approx(1.)
    assert sum(atom['f_minus'] for atom in report['atoms']) == pytest.approx(1.)

