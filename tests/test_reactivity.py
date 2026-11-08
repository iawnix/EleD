"""Analytical signs, conventions, provenance and charge-state failure behavior."""
from copy import deepcopy
from dataclasses import replace
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem

from eled import (XTBOptions, XTBResult, XTBRunError, calc_reactivity_descriptors,
                       run_reactivity, ChargeStateCalculationError)
from eled.input import _input_metadata,_xyz_from_mol
from eled.options import calculation_settings
from eled.backends.native import HARTREE_TO_EV


def states(ip=10.,ea=2.):
    mol=Chem.AddHs(Chem.MolFromSmiles('O'));conf=Chem.Conformer(3);conf.Set3D(True)
    for i,xyz in enumerate(((0,0,0),(1,0,0),(0,1,0))):conf.SetAtomPosition(i,xyz)
    mol.AddConformer(conf)
    results=[]
    for charge,spin,ne,energy,q,occ in [
        (0,0,8,-10.,(-.6,.3,.3),(2,2,2,2,0)),
        (1,1,7,-10.+ip/HARTREE_TO_EV,(-.3,.65,.65),(2,2,2,1,0)),
        (-1,1,9,-10.-ea/HARTREE_TO_EV,(-.8,-.1,-.1),(2,2,2,2,1))]:
        meta=_input_metadata(mol,-1,charge,spin,_xyz_from_mol(mol,-1))
        meta.update(engine={'sha256':'same-engine'},parameter_file={'sha256':'same-parameters'},
                    calculation_settings=calculation_settings(electronic_temperature=0.),calculation_status='success')
        results.append(XTBResult('GFN2-xTB',None,energy,None,5,ne,spin,
            orbital_energies_ev=(-10.,-8.,-6.,-4.,1.),orbital_occupations=occ,partial_charges=q,metadata=meta))
    return mol,results


def test_vertical_formulas_fukui_signs_and_fragment_conservation():
    mol,r=states(); report=calc_reactivity_descriptors(mol,*r,fragments={'O':[0],'H':[1,2]})
    assert report['values']==pytest.approx({'vertical_ip_ev':10.,'vertical_ea_ev':2.,
        'chemical_potential_ev':-6.,'electronegativity_ev':6.,'hardness_ev':4.,
        'softness_per_ev':.125,'electrophilicity_ev':4.5})
    np.testing.assert_allclose([a['f_plus'] for a in report['atoms']],[.2,.4,.4])
    np.testing.assert_allclose([a['f_minus'] for a in report['atoms']],[.3,.35,.35])
    np.testing.assert_allclose([a['dual'] for a in report['atoms']],[-.1,.05,.05])
    assert sum(v['f_plus'] for v in report['fragments'].values())==pytest.approx(1)
    assert sum(v['dual'] for v in report['fragments'].values())==pytest.approx(0)
    assert report['quality']['spin_states_verified_as_ground_states'] is False


def test_negative_ea_is_retained_and_nonpositive_hardness_not_divided():
    mol,r=states(ea=-2.);report=calc_reactivity_descriptors(mol,*r)
    assert report['values']['vertical_ea_ev']==pytest.approx(-2.)
    assert 'negative_vertical_ea' in {f['code'] for f in report['quality']['flags']}
    mol,r=states(ip=2.,ea=4.);report=calc_reactivity_descriptors(mol,*r)
    assert report['values']['hardness_ev']==pytest.approx(-1)
    assert report['values']['softness_per_ev'] is None
    assert report['values']['electrophilicity_ev'] is None


@pytest.mark.parametrize('damage',['geometry','charge','method','engine','parameters','settings','electron','conservation'])
def test_mixed_or_invalid_charge_states_are_rejected(damage):
    mol,r=states(); meta=deepcopy(r[1].metadata)
    if damage=='geometry': meta['geometry_sha256']='different'
    elif damage=='charge':meta['charge']=2
    elif damage=='engine':meta['engine']['sha256']='different'
    elif damage=='parameters':meta['parameter_file']['sha256']='different'
    elif damage=='settings':meta['calculation_settings']['accuracy']=.1
    r[1]=replace(r[1],metadata=meta)
    if damage=='method':r[1]=replace(r[1],method='GFN1-xTB')
    if damage=='electron':r[1]=replace(r[1],number_of_electrons=8)
    if damage=='conservation':r[1]=replace(r[1],partial_charges=(0.,0.,0.))
    with pytest.raises(ValueError):calc_reactivity_descriptors(mol,*r)


def test_failed_required_state_never_yields_partial_reactivity():
    mol,r=states()
    with patch('eled.reactivity._resolve_options',side_effect=lambda x:x), \
         patch('eled.reactivity.XTBCalculator.calculate',side_effect=[r[0],XTBRunError('SCF failure')]) as calculate:
        report=run_reactivity(mol,return_errors=True,initialization="independent")
    assert calculate.call_count==2
    assert report['states']['attached']['status']=='not_run'
    assert report['states']['reference']['energy_hartree']==-10.
    assert report['quality']['status']=='error'
    assert report['analysis_settings']['unpaired_electrons']=={'reference':0,'ionized':1,'attached':1}
    assert report['analysis_settings']['requested_options']['electronic_temperature']==300.
    assert all(v is None for v in report['values'].values())
    with patch('eled.reactivity._resolve_options',side_effect=lambda x:x), \
         patch('eled.reactivity.XTBCalculator.calculate',side_effect=[r[0],XTBRunError('SCF failure')]):
        with pytest.raises(ChargeStateCalculationError) as error:run_reactivity(mol,initialization="independent")
    assert error.value.state=='ionized'


def test_different_valid_conformers_are_not_vertical_charge_states():
    mol,r=states()
    conf=Chem.Conformer(mol.GetConformer());conf.SetAtomPosition(1,(1.1,0.,0.))
    index=mol.AddConformer(conf,assignId=True)
    meta=deepcopy(r[1].metadata)
    meta.update(_input_metadata(mol,index,1,1,_xyz_from_mol(mol,index)))
    r[1]=replace(r[1],metadata=meta)
    with pytest.raises(ValueError,match='fixed geometry'):
        calc_reactivity_descriptors(mol,*r)


@pytest.mark.parametrize('change',[{'number_of_unpaired_electrons':3},
                                 {'metadata_status':'error'}])
def test_invalid_spin_or_failed_state_is_rejected(change):
    mol,r=states()
    if 'metadata_status' in change:
        r[1]=replace(r[1],metadata={**r[1].metadata,'calculation_status':change['metadata_status']})
    else:
        r[1]=replace(r[1],**change)
    with pytest.raises(ValueError): calc_reactivity_descriptors(mol,*r)


def test_open_shell_and_gfn0_policies_are_checked_before_engine_lookup():
    mol,r=states()
    with patch('eled.reactivity._resolve_options',side_effect=AssertionError('must not load')):
        with pytest.raises(ValueError,match='open-shell'):
            run_reactivity(mol,options=XTBOptions(charge=1,uhf=1))
        with pytest.raises(ValueError,match='GFN1/2'):
            run_reactivity(mol,options=XTBOptions(method='gfn0'))
        with pytest.raises(ValueError,match='electron count'):
            run_reactivity(mol,ionized_uhf=0)
        with pytest.raises(ValueError,match='initialization'):
            run_reactivity(mol,initialization='auto')


def test_default_protocol_uses_300k_reference_copies_and_records_version():
    mol,r=states()
    r=[replace(v,metadata={**v.metadata,'calculation_settings':calculation_settings()}) for v in r]
    with patch('eled.reactivity._resolve_options',side_effect=lambda x:x), \
         patch('eled.reactivity._charge_state_results',return_value=iter(r)) as calculate:
        # A generator is required because early returns explicitly close it.
        calculate.return_value=(v for v in r)
        report=run_reactivity(mol)
    _,configs,interface,initialization=calculate.call_args.args
    assert all(v.electronic_temperature==300. for v in configs)
    assert initialization=='reference' and interface=='native'
    assert report['analysis_settings']['protocol_version']=='2.0'
    assert report['analysis_settings']['initialization']=='reference'
