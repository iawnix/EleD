"""Interruptions, protocol isolation, transaction integrity and fatal errors."""
from copy import deepcopy
from dataclasses import replace
import json
import sqlite3
from unittest.mock import patch

import pytest
from rdkit import Chem

from rdkit_eth import (XTBOptions, XTBResult, XTBDiskCache, ReactivityCheckpoint,
    ReactivityJob, run_reactivity_batch, XTBCacheError, XTBEnvironmentError,
    run_reactivity, export_descriptor_reports)
from rdkit_eth.reactivity import _empty_report, REACTIVITY_UNITS


def molecule(smi='O'):
    mol=Chem.AddHs(Chem.MolFromSmiles(smi))
    conf=Chem.Conformer(mol.GetNumAtoms());conf.Set3D(True)
    for i in range(mol.GetNumAtoms()):conf.SetAtomPosition(i,(i,i%2,.1*i))
    mol.AddConformer(conf)
    return mol


def report(failed=False):
    r=_empty_report(None,'charges',{k:{'status':'success'} for k in ('reference','ionized','attached')})
    r['quality']={'status':'error' if failed else 'ok','flags':[]}
    if failed:
        r['error']={'type':'ChargeStateCalculationError','state':'ionized','message':'SCF failed'}
        r['states']['ionized']['status']='error';r['states']['attached']['status']='not_run'
    else:r['values']={k:1. for k in REACTIVITY_UNITS}
    return r


@pytest.fixture
def environment(monkeypatch):
    import rdkit_eth.reactivity_batch as batch
    monkeypatch.setattr(batch,'_resolve_options',lambda c:replace(c,xtb_library='/fake/libxtb.so'))
    monkeypatch.setattr(batch,'_engine_identity',lambda *a:{'sha256':'engine-v1','orbital_basis':True})
    monkeypatch.setattr(batch,'_parameter_file',lambda *a:{'sha256':'parameters-v1'})
    monkeypatch.setattr(batch,'_parameter_context',lambda :{})
    monkeypatch.setattr(batch,'_implementation_identity',lambda *a:{'sha256':'adapter-v1'})


def test_disk_cache_roundtrips_arrays_across_instances_and_rejects_corruption(tmp_path):
    path=tmp_path/'cache.sqlite'
    result=XTBResult('GFN2-xTB',None,-1.,2.,2,2,0,orbital_energies_ev=(-1.,1.),
        orbital_occupations=(2.,0.),overlap_matrix=((1.,0.),(0.,1.)),
        metadata={'calculation_status':'success','nested':{'value':1}})
    XTBDiskCache(path).set('single',result)
    loaded=XTBDiskCache(path).get('single')
    assert loaded==result and isinstance(loaded.overlap_matrix[0],tuple)
    loaded.metadata['nested']['value']=9
    assert XTBDiskCache(path).get('single').metadata['nested']['value']==1
    with sqlite3.connect(path) as db:
        db.execute("UPDATE records SET payload='{}' WHERE key='single'")
    with pytest.raises(XTBCacheError,match='checksum'):XTBDiskCache(path).get('single')


def test_cache_refuses_failed_nonfinite_and_unrelated_database(tmp_path):
    cache=XTBDiskCache(tmp_path/'cache.sqlite')
    with pytest.raises(ValueError):cache.set('failed',report(True))
    bad=report();bad['values']['vertical_ip_ev']=float('nan')
    with pytest.raises(ValueError):cache.set('bad',bad)
    assert len(cache)==0
    path=tmp_path/'unrelated.sqlite'
    with sqlite3.connect(path) as db:db.execute('CREATE TABLE keep_me (v TEXT)')
    with pytest.raises(XTBCacheError):XTBDiskCache(path)
    with sqlite3.connect(path) as db:assert db.execute('SELECT * FROM keep_me').fetchall()==[]


def test_interruption_commits_only_completed_rows_and_resumes(environment,tmp_path):
    jobs=[ReactivityJob('water',molecule()),ReactivityJob('ammonia',molecule('N'))]
    checkpoint=tmp_path/'checkpoint.sqlite'
    with patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=[report(),KeyboardInterrupt]):
        with pytest.raises(KeyboardInterrupt):run_reactivity_batch(jobs,checkpoint=checkpoint)
    assert len(ReactivityCheckpoint(checkpoint).reports())==1
    with patch('rdkit_eth.reactivity_batch.run_reactivity',return_value=report()) as calculate:
        rows=run_reactivity_batch(jobs,checkpoint=checkpoint)
    assert calculate.call_count==1
    assert [r['record_id'] for r in rows]==['water','ammonia']
    assert rows[0]['execution']['source']=='checkpoint'
    with patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=AssertionError('recalculated')):
        run_reactivity_batch(jobs,checkpoint=checkpoint)
    with pytest.raises(XTBCacheError,match='differs'):
        run_reactivity_batch(jobs[::-1],checkpoint=checkpoint)
    with pytest.raises(XTBCacheError,match='differs'):
        run_reactivity_batch(jobs,checkpoint=checkpoint,initialization='independent')


def test_failure_history_is_preserved_and_only_explicitly_retried(environment,tmp_path):
    jobs=[ReactivityJob('water',molecule())];checkpoint=ReactivityCheckpoint(tmp_path/'checkpoint.sqlite')
    cache=XTBDiskCache(tmp_path/'cache.sqlite')
    with patch('rdkit_eth.reactivity_batch.run_reactivity',return_value=report(True)):
        rows=run_reactivity_batch(jobs,checkpoint=checkpoint,cache=cache)
    assert rows[0]['error'] and len(cache)==0
    with patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=AssertionError('retried')):
        assert run_reactivity_batch(jobs,checkpoint=checkpoint)[0]['error']
    with patch('rdkit_eth.reactivity_batch.run_reactivity',return_value=report()) as calculate:
        rows=run_reactivity_batch(jobs,checkpoint=checkpoint,cache=cache,retry_failed=True)
    assert calculate.call_count==1 and len(cache)==1 and rows[0]['execution']['attempt']==2
    history=checkpoint.history('water')
    assert len(history)==2 and history[0]['report']['error'] and history[1]['report']['error'] is None


def test_cached_reports_are_independent_of_record_ids_and_export_execution(environment,tmp_path):
    cache=XTBDiskCache(tmp_path/'cache.sqlite')
    jobs=[ReactivityJob('a',molecule()),ReactivityJob('b',molecule())]
    with patch('rdkit_eth.reactivity_batch.run_reactivity',return_value=report()) as calculate:
        rows=run_reactivity_batch(jobs,cache=cache)
    assert calculate.call_count==1 and rows[1]['execution']['source']=='cache'
    assert [r['record_id'] for r in rows]==['a','b']
    assert list(rows[0]['units'])==list(rows[1]['units'])
    rows[0]['values']['vertical_ip_ev']=100
    assert rows[1]['values']['vertical_ip_ev']==1
    path=export_descriptor_reports(rows,tmp_path/'rows.csv')
    import csv
    exported=list(csv.DictReader(path.open()))
    assert json.loads(exported[1]['execution_json'])['source']=='cache'


@pytest.mark.parametrize('change',['temperature','accuracy','iterations','method','initialization',
                                  'geometry','spin','fragments','scheme','interface','engine','parameters','adapter'])
def test_protocol_changes_do_not_hit_old_cache(environment,tmp_path,monkeypatch,change):
    import rdkit_eth.reactivity_batch as batch
    cache=XTBDiskCache(tmp_path/'cache.sqlite');job=ReactivityJob('a',molecule())
    kwargs={}
    with patch.object(batch,'run_reactivity',side_effect=lambda *a,**k:report()) as calculate:
        run_reactivity_batch([job],cache=cache)
        if change=='temperature':kwargs['options']=XTBOptions(electronic_temperature=0.)
        if change=='accuracy':kwargs['options']=XTBOptions(accuracy=.1)
        if change=='iterations':kwargs['options']=XTBOptions(max_iterations=300)
        if change=='method':kwargs['options']=XTBOptions(method='gfn1')
        if change=='initialization':kwargs['initialization']='independent'
        if change=='scheme':kwargs['population_scheme']='lowdin'
        if change=='interface':kwargs['interface']='compiled'
        if change=='geometry':job.mol.GetConformer().SetAtomPosition(1,(3.,2.,1.))
        if change=='spin':job=replace(job,ionized_uhf=3)
        if change=='fragments':job=replace(job,fragments={'oxygen':[0]})
        if change=='engine':monkeypatch.setattr(batch,'_engine_identity',lambda *a:{'sha256':'engine-v2','orbital_basis':True})
        if change=='parameters':monkeypatch.setattr(batch,'_parameter_file',lambda *a:{'sha256':'parameters-v2'})
        if change=='adapter':monkeypatch.setattr(batch,'_implementation_identity',lambda *a:{'sha256':'adapter-v2'})
        run_reactivity_batch([job],cache=cache,**kwargs)
    assert calculate.call_count==2 and len(cache)==2


def test_invalid_geometry_is_a_row_but_spin_configuration_is_fatal(environment):
    missing=Chem.AddHs(Chem.MolFromSmiles('O'))
    jobs=[ReactivityJob('missing',missing),ReactivityJob('parse-failed',None),ReactivityJob('empty',Chem.Mol())]
    with patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=AssertionError('must not calculate')):
        rows=run_reactivity_batch(jobs)
        assert all(r['error']['stage']=='input' and all(v is None for v in r['values'].values()) for r in rows)
        with pytest.raises(ValueError,match='open-shell'):
            run_reactivity_batch([ReactivityJob('good',molecule()),ReactivityJob('radical',molecule('[CH3]'))])
        with pytest.raises(ValueError,match='unique'):
            run_reactivity_batch([ReactivityJob('same',molecule())]*2)


@pytest.mark.parametrize('error',[XTBEnvironmentError('allocation failed'),ImportError('bad extension'),
                                RuntimeError('unexpected failure')])
def test_runtime_failures_abort_without_recording_a_molecular_failure(environment,tmp_path,error):
    cp=ReactivityCheckpoint(tmp_path/'checkpoint.sqlite')
    with patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=error):
        with pytest.raises(type(error)):run_reactivity_batch([ReactivityJob('a',molecule())],checkpoint=cp)
    assert cp.reports()==[]


def test_checkpoint_concurrent_writer_and_tampered_row_are_rejected(environment,tmp_path):
    cp=ReactivityCheckpoint(tmp_path/'checkpoint.sqlite');jobs=[ReactivityJob('a',molecule())]
    with cp._locked():
        with pytest.raises(XTBCacheError,match='writer'):run_reactivity_batch(jobs,checkpoint=cp.path)
    with patch('rdkit_eth.reactivity_batch.run_reactivity',return_value=report()):run_reactivity_batch(jobs,checkpoint=cp)
    with sqlite3.connect(cp.path) as db:db.execute("UPDATE records SET payload='{}' WHERE namespace='rows'")
    with pytest.raises(XTBCacheError,match='checksum'):run_reactivity_batch(jobs,checkpoint=cp)


def test_single_reactivity_does_not_hide_runtime_errors():
    with patch('rdkit_eth.reactivity._resolve_options',side_effect=lambda c:c), \
         patch('rdkit_eth.reactivity.XTBCalculator.calculate',side_effect=ImportError('broken runtime')):
        with pytest.raises(ImportError):run_reactivity(molecule(),return_errors=True,initialization='independent')


def test_concurrent_disk_writes_are_complete(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    path=tmp_path/'concurrent.sqlite'
    def write(index):
        cache=XTBDiskCache(path)
        cache.set(str(index),report())
        assert cache.get(str(index))['error'] is None
    with ThreadPoolExecutor(max_workers=4) as executor:
        list(executor.map(write,range(12)))
    assert len(XTBDiskCache(path))==12


def test_checkpoint_write_error_aborts_and_keeps_prior_commit(environment,tmp_path):
    cp=ReactivityCheckpoint(tmp_path/'checkpoint.sqlite')
    jobs=[ReactivityJob('water',molecule()),ReactivityJob('ammonia',molecule('N'))]
    original=cp._save
    def save(identifier,*args):
        if identifier=='ammonia':raise sqlite3.OperationalError('disk full')
        original(identifier,*args)
    with patch.object(cp,'_save',save),patch('rdkit_eth.reactivity_batch.run_reactivity',side_effect=lambda *a,**k:report()):
        with pytest.raises(sqlite3.OperationalError,match='disk full'):
            run_reactivity_batch(jobs,checkpoint=cp)
    assert [r['record_id'] for r in cp.reports()]==['water']
