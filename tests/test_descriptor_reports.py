"""Quality semantics, fixed table schemas and transparent missing/failure rows."""
from dataclasses import replace
import csv
import json
from unittest.mock import patch

import pytest
from rdkit import Chem

from eled import (
    XTBOptions, XTBResults, XTBRunError, DescriptorColumns, CalcDescriptorReport,
    AssessDescriptorQuality, ExportDescriptorReports, RunDescriptorBatch,
)
from eled.input import _input_metadata, _xyz_from_mol
from eled.options import _validate_calculation_settings


def sample():
    mol = Chem.AddHs(Chem.MolFromSmiles('O')); conf = Chem.Conformer(3); conf.Set3D(True)
    for i, xyz in enumerate(((0,0,0),(1,0,0),(0,1,0))): conf.SetAtomPosition(i,xyz)
    mol.AddConformer(conf)
    meta = _input_metadata(mol,-1,0,0,_xyz_from_mol(mol,-1))
    meta.update(calculation_settings=_validate_calculation_settings(),calculation_status='success')
    r = XTBResults('GFN2-xTB',None,-1.,3.,3,4,0,orbitalEnergiesEV=(-10.,-2.,1.),
        orbitalOccupancies=(2.,2.,0.),atomicCharges=(-.6,.3,.3),dipoleAtomicUnits=(0.,0.,1.),
        wibergBondOrders=((0.,1.,1.),(1.,0.,0.),(1.,0.,0.)),metadata=meta)
    return mol,r


def test_report_records_missing_basis_and_exports_fixed_columns(tmp_path):
    mol,r = sample(); report = CalcDescriptorReport(mol,r,recordId='water')
    assert report['quality']['status']=='partial'
    assert report['missing']['populations']=='basis_not_retained'
    assert report['values']['homo_pr'] is None
    assert set(report['values'])==set(DescriptorColumns())
    assert report['units']['xtb_energy_ha']=='Hartree'
    assert report['units']['xtb_homo_ev']=='eV'
    assert report['quality']['scf_iterations'] is None
    for kind in ['csv','jsonl']:
        p=ExportDescriptorReports([report],tmp_path/f'data.{kind}')
        if kind=='csv':
            row=list(csv.DictReader(p.open()))[0]
            assert row['record_id']=='water' and row['homo_pr']==''
            assert json.loads(row['missing_json'])['frontier']=='basis_not_retained'
            assert json.loads(row['provenance_json'])['geometry_sha256']==r.metadata['geometry_sha256']
        else: assert json.loads(p.read_text())==report
    report['provenance']['charge']=8
    assert r.metadata['charge']==0


def test_report_records_mci_request_when_basis_is_not_retained():
    mol, r = sample()
    report = CalcDescriptorReport(mol, r, mciRings={'candidate': [0, 1, 2]})
    assert report['mci'] is None
    assert report['missing']['mci'] == 'basis_not_retained'
    assert report['analysis_settings']['mciRings'] == {'candidate': [0, 1, 2]}


def test_quality_flags_fractional_and_smeared_frontiers():
    _, r = sample()
    quality=AssessDescriptorQuality(replace(r,orbitalOccupancies=(2.,2.,5e-7)))
    assert 'fractional_frontier_occupation' not in {f['code'] for f in quality['flags']}
    assert quality['frontier_occupations']=={'homo':2.,'lumo':5e-7}
    quality=AssessDescriptorQuality(replace(r,orbitalOccupancies=(2.,1.9,.1)))
    assert 'fractional_frontier_occupation' in {f['code'] for f in quality['flags']}
    assert quality['frontier_occupations']=={'homo':1.9,'lumo':.1}
    quality=AssessDescriptorQuality(replace(r,orbitalOccupancies=(2.,1.5,.5)))
    assert quality['fractional_orb_idx']==[1,2]
    assert quality['status']=='warning'
    q=AssessDescriptorQuality(replace(r,atomicCharges=(1.,1.,1.)))
    assert q['status']=='error' and q['charge_sum_error_e']==3
    q=AssessDescriptorQuality(replace(r,metadata={}))
    assert 'calculation_settings_unknown' in {f['code'] for f in q['flags']}
    with pytest.raises(ValueError,match='exceed two'):
        AssessDescriptorQuality(replace(r,orbitalOccupancies=(3.,1.,0.)))


def test_export_rejects_schema_mismatch_before_overwriting(tmp_path):
    mol,r=sample(); a=CalcDescriptorReport(mol,r); b=CalcDescriptorReport(mol,r,frontierDepth=3)
    path=tmp_path/'data.csv';path.write_text('keep')
    with pytest.raises(ValueError,match='incompatible'): ExportDescriptorReports([a,b],path)
    assert path.read_text()=='keep'
    a['values']['xtb_gap_ev']=float('nan')
    with pytest.raises(ValueError): ExportDescriptorReports([a],path)
    assert path.read_text()=='keep'


def test_batch_failures_preserve_order_ids_and_null_features():
    mol,r=sample()
    with patch('eled.batch.RunXTBBatch',return_value=[r,XTBRunError('no convergence'),r]):
        reports=RunDescriptorBatch([mol]*3,recordIds=['a','b','c'])
    assert [v['record_id'] for v in reports]==['a','b','c']
    assert reports[1]['quality']['status']=='error'
    assert reports[1]['error']['type']=='XTBRunError'
    assert all(v is None for v in reports[1]['values'].values())
    assert reports[0]['units']==reports[1]['units']
    assert reports[0]['analysis_settings']==reports[1]['analysis_settings']
    assert reports[1]['provenance']['requested_options']['accuracy']==1.


@pytest.mark.parametrize('kwargs',[{'typo':1},{'degeneracyEV':-1},
                                  {'occupationTolerance':float('nan')}, {'includeOrbitals':1}])
def test_batch_rejects_global_analysis_errors_before_calculation(kwargs):
    mol,_=sample()
    with patch('eled.batch.RunXTBBatch',side_effect=AssertionError('must not calculate')):
        with pytest.raises((TypeError,ValueError)):
            RunDescriptorBatch([mol],**kwargs)


@pytest.mark.parametrize('kwargs',[
    {'maxIterations':True},{'maxIterations':2.5},{'maxIterations':0},
    {'maxIterations':2**32},{'accuracy':True},{'accuracy':1e-5},
    {'accuracy':1001.},{'accuracy':float('nan')},{'electronicTemperature':-1},
    {'electronicTemperature':float('inf')},{'electronicTemperature':True},
])
def test_calculation_options_reject_clamping_and_invalid_numbers(kwargs):
    with pytest.raises(ValueError): XTBOptions(**kwargs)


def test_zero_temperature_records_actual_engine_floor():
    settings=_validate_calculation_settings(electronicTemperature=0)
    assert settings['electronic_temperature']==0
    assert settings['effective_electronic_temperature']==1e-6


@pytest.mark.parametrize('kwargs',[{'maxIterations':True}, {'maxIterations':2**32},
                                  {'accuracy':1e-5}, {'electronicTemperature':float('nan')}])
def test_compiled_settings_reject_invalid_values_before_library_lookup(kwargs):
    module=pytest.importorskip('eled.rdXTBTools')
    mol,_=sample()
    with pytest.raises(ValueError):
        module.RunMol(mol,xtbLibrary='/missing/libxtb.so',**kwargs)
