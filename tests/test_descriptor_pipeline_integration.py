"""Real-engine parameter parity, finite-temperature sensitivity and charge states."""
from helpers import calculate_xtb
from dataclasses import replace
import json
import os
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (XTBOptions, XTBNativeNotFoundError, RunXTBBatch, XTBCache,
    CalcDescriptorReport, RunReactivity, CalcReactivityDescriptors, AssessDescriptorQuality)
from eled.runtime import _find_xtb_library, _native_basis_build

pytestmark = pytest.mark.integration


@pytest.fixture(scope='module')
def library():
    try: lib = _find_xtb_library()
    except XTBNativeNotFoundError: pytest.skip('libxtb unavailable')
    if not _native_basis_build(lib): pytest.skip('native basis extension required')
    return lib


def molecule(smiles='O'):
    m=Chem.AddHs(Chem.MolFromSmiles(smiles)); assert AllChem.EmbedMolecule(m,randomSeed=1)==0
    return m


@pytest.mark.parametrize('method',['gfn0','gfn1','gfn2'])
@pytest.mark.parametrize('temperature',[0.,1200.])
def test_nondefault_settings_match_compiled_and_are_recorded(library,method,temperature):
    mod=pytest.importorskip('eled.rdXTBTools');m=molecule('CCO')
    a=calculate_xtb(m,method=method,xtb_library=library,accuracy=.1,max_iterations=333,
              electronic_temperature=temperature,keep_overlap=True)
    ok,b=mod.RunMol(m,method=method,xtbLibrary=library,accuracy=.1,maxIterations=333,
                   electronicTemperature=temperature,keepOverlap=True)
    assert ok
    ra,rb=CalcDescriptorReport(m,a),CalcDescriptorReport(m,b)
    assert ra['provenance']['CalculationSettings']==rb['provenance']['CalculationSettings']
    assert b.maxIterations==333 and b.accuracy==.1
    assert b.effectiveElectronicTemperature==max(temperature,1e-6)
    np.testing.assert_allclose(a.orbital_energies_ev,b.GetOrbitalEnergies(),atol=1e-9,rtol=0)
    for field in ra['values']:
        va,vb=ra['values'][field],rb['values'][field]
        if field=='xtb_backend': continue
        if va is None or isinstance(va,str): assert va==vb
        else: assert va==pytest.approx(vb,abs=1e-9)
    assert json.loads(json.dumps(ra,allow_nan=False))['quality']['calculation_status']=='success'


def test_cache_distinguishes_each_numerical_setting(library):
    m=molecule(); cache=XTBCache(); opts=XTBOptions(xtb_library=library)
    base=RunXTBBatch([m],options=opts,cache=cache)[0]
    for change in [{'accuracy':.2},{'max_iterations':300},{'electronic_temperature':1000.}]:
        result=RunXTBBatch([m],options=replace(opts,**change),cache=cache)[0]
        assert result.metadata['CalculationSettings'] != base.metadata['CalculationSettings']
    assert len(cache)==4


@pytest.mark.parametrize('method',['gfn1','gfn2'])
@pytest.mark.parametrize('smiles',['O','CCO','[CH3]','[NH4+]'])
@pytest.mark.parametrize('initialization',['reference','independent'])
def test_vertical_states_interfaces_and_population_definitions(library,method,smiles,initialization):
    pytest.importorskip('eled.rdXTBTools');m=molecule(smiles);original=m.ToBinary()
    opts=XTBOptions(method=method,xtb_library=library,accuracy=.1,electronic_temperature=300.)
    spins={'ionized_uhf':0,'attached_uhf':0} if smiles=='[CH3]' else {}
    a=RunReactivity(m,options=opts,population_scheme='mulliken',initialization=initialization,**spins)
    b=RunReactivity(m,options=opts,population_scheme='charges',interface='compiled',initialization=initialization,**spins)
    assert a['values']==pytest.approx(b['values'],abs=1e-8)
    for key in ['f_plus','f_minus','f_radical','dual']:
        np.testing.assert_allclose([v[key] for v in a['atoms']],[v[key] for v in b['atoms']],atol=1e-7,rtol=0)
    assert sum(r['f_plus'] for r in a['atoms'])==pytest.approx(1,abs=1e-7)
    assert sum(r['f_minus'] for r in a['atoms'])==pytest.approx(1,abs=1e-7)
    assert m.ToBinary()==original
    assert len({s['provenance']['geometry_sha256'] for s in a['states'].values()})==1
    for report in (a,b):
        assert report['analysis_settings']['initialization']==initialization
        for state,details in report['states'].items():
            assert details['provenance']['initialization']==('independent' if state=='reference' else initialization)


def test_reference_restart_is_copied_for_each_neighbor_and_freed(library):
    from eled.backends.native import _NativeXTB, _native_charge_states
    original=_NativeXTB._bind
    copies,deleted=[],[]
    def bind(instance,name,*args):
        fn=original(instance,name,*args)
        if name=='xtb_copyResults':
            def copy(source):
                result=fn(source);copies.append((source,result));return result
            return copy
        if name=='xtb_delResults':
            def delete(pointer):
                deleted.append(pointer._obj.value);return fn(pointer)
            return delete
        return fn
    m=molecule();base=XTBOptions(xtb_library=library)
    configs=[replace(base,charge=q,uhf=spin) for q,spin in [(0,0),(1,1),(-1,1)]]
    with patch.object(_NativeXTB,'_bind',bind):
        results=list(_native_charge_states(m,configs))
    assert len(copies)==3 and copies[1][0]==copies[2][0]==copies[0][1]
    assert copies[0][1] in deleted and len(deleted)==4
    # Reversing neighboring calculations must not turn this into a sequential restart.
    reversed_results=list(_native_charge_states(m,[configs[0],configs[2],configs[1]]))
    for a,b in zip(results,[reversed_results[i] for i in [0,2,1]]):
        np.testing.assert_allclose(a.partial_charges,b.partial_charges,atol=1e-10,rtol=0)
    copies.clear();deleted.clear()
    with patch.object(_NativeXTB,'_bind',bind):
        iterator=_native_charge_states(m,configs)
        next(iterator);iterator.close()
    assert copies[0][1] in deleted and len(deleted)==2


@pytest.mark.parametrize('fail_on,state',[(1,'reference'),(2,'ionized'),(3,'attached')])
def test_restart_copy_allocation_failure_is_reported_and_recovers(library,fail_on,state):
    from eled.backends.native import _NativeXTB
    original=_NativeXTB._bind
    count=0
    def bind(instance,name,*args):
        fn=original(instance,name,*args)
        if name=='xtb_copyResults':
            def copy(source):
                nonlocal count
                count+=1
                return None if count==fail_on else fn(source)
            return copy
        return fn
    from eled import XTBEnvironmentError
    with patch.object(_NativeXTB,'_bind',bind):
        with pytest.raises(XTBEnvironmentError,match='copy'):
            RunReactivity(molecule(),options=XTBOptions(xtb_library=library),return_errors=True)
    assert RunReactivity(molecule(),options=XTBOptions(xtb_library=library))['error'] is None


@pytest.mark.parametrize('interface',['native','compiled'])
def test_actual_scf_failure_keeps_protocol_and_does_not_poison_next_call(library,interface):
    if interface=='compiled': pytest.importorskip('eled.rdXTBTools')
    report=RunReactivity(molecule('CCO'), interface=interface, return_errors=True,
        options=XTBOptions(xtb_library=library,max_iterations=1,electronic_temperature=0.))
    assert report['error']['state']=='reference'
    assert 'converge' in report['error']['message']
    assert all(v is None for v in report['values'].values())
    assert report['states']['ionized']['status']=='not_run'
    assert report['analysis_settings']['requested_options']['max_iterations']==1
    successful=RunReactivity(molecule(),interface=interface,
        options=XTBOptions(xtb_library=library,electronic_temperature=0.))
    assert successful['error'] is None
    assert successful['values']['ip_v_ev'] > 0


def test_lowdin_fukui_and_finite_temperature_label(library):
    m=molecule('CCO')
    r=RunReactivity(m,options=XTBOptions(xtb_library=library,electronic_temperature=1000.),
                     population_scheme='lowdin',fragments={'all':list(range(m.GetNumAtoms()))})
    assert r['fragments']['all']['f_plus']==pytest.approx(1.,abs=1e-7)
    assert 'finite_electronic_temperature' in {f['code'] for f in r['quality']['flags']}


def test_report_does_not_recalculate_and_labels_thermal_occupations(library):
    m=molecule('c1ccccc1')
    result=calculate_xtb(m,xtb_library=library,electronic_temperature=50000.,keep_overlap=True)
    with patch('eled.xtb._run_xtb',side_effect=AssertionError('must not calculate')):
        report=CalcDescriptorReport(m,result)
    assert report['quality']['fractional_orb_idx']
    assert 'fractional_occupations' in {f['code'] for f in report['quality']['flags']}
