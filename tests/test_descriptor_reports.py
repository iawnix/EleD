"""Quality semantics, fixed table schemas and transparent missing/failure rows."""
from dataclasses import replace
import csv
import json
from pathlib import Path
from unittest.mock import patch

import pytest
from rdkit import Chem

from rdkit_eth import (
    XTBOptions, XTBResult, XTBRunError, descriptor_columns, calc_descriptor_report,
    assess_descriptor_quality, export_descriptor_reports, run_descriptor_batch,
)
from rdkit_eth.input import _input_metadata, _xyz_from_mol
from rdkit_eth.options import calculation_settings


def sample():
    mol = Chem.AddHs(Chem.MolFromSmiles('O')); conf = Chem.Conformer(3); conf.Set3D(True)
    for i, xyz in enumerate(((0,0,0),(1,0,0),(0,1,0))): conf.SetAtomPosition(i,xyz)
    mol.AddConformer(conf)
    meta = _input_metadata(mol,-1,0,0,_xyz_from_mol(mol,-1))
    meta.update(calculation_settings=calculation_settings(),calculation_status='success')
    r = XTBResult('GFN2-xTB',None,-1.,3.,3,4,0,orbital_energies_ev=(-10.,-2.,1.),
        orbital_occupations=(2.,2.,0.),partial_charges=(-.6,.3,.3),dipole_e_bohr=(0.,0.,1.),
        wiberg_bond_orders=((0.,1.,1.),(1.,0.,0.),(1.,0.,0.)),metadata=meta)
    return mol,r


def test_report_records_missing_basis_and_exports_fixed_columns(tmp_path):
    mol,r = sample(); report = calc_descriptor_report(mol,r,record_id='water')
    assert report['quality']['status']=='partial'
    assert report['missing']['populations']=='basis_not_retained'
    assert report['values']['frontier_homo_participation_ratio'] is None
    assert set(report['values'])==set(descriptor_columns())
    assert report['units']['xtb_total_energy_hartree']=='Hartree'
    assert report['units']['xtb_homo_ev']=='eV'
    assert report['quality']['scf_iterations'] is None
    for kind in ['csv','jsonl']:
        p=export_descriptor_reports([report],tmp_path/f'data.{kind}')
        if kind=='csv':
            row=list(csv.DictReader(p.open()))[0]
            assert row['record_id']=='water' and row['frontier_homo_participation_ratio']==''
            assert json.loads(row['missing_json'])['frontier']=='basis_not_retained'
            assert json.loads(row['provenance_json'])['geometry_sha256']==r.metadata['geometry_sha256']
        else: assert json.loads(p.read_text())==report
    report['provenance']['charge']=8
    assert r.metadata['charge']==0


def test_quality_flags_fractional_and_threshold_sensitive_frontiers():
    _, r = sample()
    mixed = replace(r,orbital_occupations=(2.,1.9999995,5e-7))
    quality=assess_descriptor_quality(mixed)
    assert 'frontier_threshold_sensitive' in {f['code'] for f in quality['flags']}
    quality=assess_descriptor_quality(replace(r,orbital_occupations=(2.,1.5,.5)))
    assert quality['fractional_orbital_indices']==[1,2]
    assert quality['status']=='warning'
    q=assess_descriptor_quality(replace(r,partial_charges=(1.,1.,1.)))
    assert q['status']=='error' and q['charge_sum_error_e']==3
    q=assess_descriptor_quality(replace(r,metadata={}))
    assert 'calculation_settings_unknown' in {f['code'] for f in q['flags']}
    with pytest.raises(ValueError,match='exceed two'):
        assess_descriptor_quality(replace(r,orbital_occupations=(3.,1.,0.)))


def test_export_rejects_schema_mismatch_before_overwriting(tmp_path):
    mol,r=sample(); a=calc_descriptor_report(mol,r); b=calc_descriptor_report(mol,r,frontier_depth=3)
    path=tmp_path/'data.csv';path.write_text('keep')
    with pytest.raises(ValueError,match='incompatible'): export_descriptor_reports([a,b],path)
    assert path.read_text()=='keep'
    a['values']['xtb_gap_ev']=float('nan')
    with pytest.raises(ValueError): export_descriptor_reports([a],path)
    assert path.read_text()=='keep'


def test_batch_failures_preserve_order_ids_and_null_features():
    mol,r=sample()
    with patch('rdkit_eth.batch.run_xtb_batch',return_value=[r,XTBRunError('no convergence'),r]):
        reports=run_descriptor_batch([mol]*3,record_ids=['a','b','c'])
    assert [v['record_id'] for v in reports]==['a','b','c']
    assert reports[1]['quality']['status']=='error'
    assert reports[1]['error']['type']=='XTBRunError'
    assert all(v is None for v in reports[1]['values'].values())
    assert reports[0]['units']==reports[1]['units']
    assert reports[0]['analysis_settings']==reports[1]['analysis_settings']
    assert reports[1]['provenance']['requested_options']['accuracy']==1.


@pytest.mark.parametrize('kwargs',[{'typo':1},{'degeneracy_ev':-1},
                                  {'occupation_threshold':float('nan')}, {'include_orbitals':1}])
def test_batch_rejects_global_analysis_errors_before_calculation(kwargs):
    mol,_=sample()
    with patch('rdkit_eth.batch.run_xtb_batch',side_effect=AssertionError('must not calculate')):
        with pytest.raises((TypeError,ValueError)):
            run_descriptor_batch([mol],**kwargs)


@pytest.mark.parametrize('kwargs',[
    {'max_iterations':True},{'max_iterations':2.5},{'max_iterations':0},
    {'max_iterations':2**32},{'accuracy':True},{'accuracy':1e-5},
    {'accuracy':1001.},{'accuracy':float('nan')},{'electronic_temperature':-1},
    {'electronic_temperature':float('inf')},{'electronic_temperature':True},
])
def test_calculation_options_reject_clamping_and_invalid_numbers(kwargs):
    with pytest.raises(ValueError): XTBOptions(**kwargs)


def test_zero_temperature_records_actual_engine_floor():
    settings=calculation_settings(electronic_temperature=0)
    assert settings['electronic_temperature']==0
    assert settings['effective_electronic_temperature']==1e-6


@pytest.mark.parametrize('kwargs',[{'maxIterations':True}, {'maxIterations':2**32},
                                  {'accuracy':1e-5}, {'electronicTemperature':float('nan')}])
def test_compiled_settings_reject_invalid_values_before_library_lookup(kwargs):
    module=pytest.importorskip('rdkit_eth.rdXTBTools')
    mol,_=sample()
    with pytest.raises(ValueError):
        module.RunMol(mol,xtbLibrary='/missing/libxtb.so',**kwargs)
