"""Sensitivity coverage, protocol isolation, error policy and physical invariance."""
from dataclasses import replace
import json
import os

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (XTBOptions, XTBResult, XTBRunError, XTBEnvironmentError,
                       XTBCacheError, XTBDiskCache, RunConformerAnalysis,
                       ConformerDescriptorColumns)
from eled import conformers
from eled.input import _input_metadata, _xyz_from_mol
from eled.options import CalculationSettings


@pytest.fixture
def synthetic(monkeypatch):
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    for identifier, shift in [(3, 0.), (8, .1), (12, .2)]:
        conf = Chem.Conformer(3); conf.Set3D(True); conf.SetId(identifier)
        for i, xyz in enumerate(((0, 0, 0), (1+shift, 0, 0), (0, 1, 0))): conf.SetAtomPosition(i, xyz)
        mol.AddConformer(conf, assignId=False)
    monkeypatch.setattr('eled.batch._resolve_options', lambda options: options)
    def calculate(self, molecule, **kwargs):
        identifier = kwargs.get('conf_id', self.options.conf_id)
        h = {3: -3., 8: -2., 12: -1.}[identifier]
        metadata = _input_metadata(molecule, identifier, 0, 0, _xyz_from_mol(molecule, identifier))
        metadata.update(CalculationSettings=CalculationSettings(), calculation_status='success',
                        engine={'sha256': 'fixture-engine'}, initialization='independent')
        return XTBResult('GFN2-xTB', None, -1. + h * .001, 3., 3, 4, 0,
            orbital_energies_ev=(-10., h, 1.), orbital_occupations=(2., 2., 0.),
            partial_charges=(-.6, .3, .3), dipole_e_bohr=(0., 0., 1.),
            wiberg_bond_orders=((0., 1., 1.), (1., 0., 0.), (1., 0., 0.)), metadata=metadata)
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', calculate)
    return mol, calculate


def test_statistics_order_coverage_and_input_unchanged(synthetic):
    mol, _ = synthetic; before = mol.ToBinary()
    report = RunConformerAnalysis(mol, record_id='water')
    assert report['conformer_ids'] == [3, 8, 12]
    assert report['status'] == 'partial'  # synthetic results lack basis arrays
    stats = report['statistics']['xtb_homo_ev']
    assert stats == dict(n_total=3, n_valid=3, n_missing=0, n_failed=0,
                         mean=-2., std=np.std([-3., -2., -1.]), min=-3., max=-1., range=2.)
    missing = report['statistics']['homo_pr']
    assert missing['n_missing'] == 3 and missing['mean'] is None and missing['n_failed'] == 0
    reverse = RunConformerAnalysis(mol, conf_ids=[12, 8, 3])
    assert reverse['statistics'] == report['statistics']
    assert mol.ToBinary() == before
    assert json.loads(json.dumps(report, allow_nan=False)) == report
    single = RunConformerAnalysis(mol, conf_ids=[8], columns=['xtb_homo_ev'])
    assert single['status'] == 'complete' and single['statistics']['xtb_homo_ev']['std'] == 0


def test_boltzmann_statistics_are_separate_and_energy_weighted(synthetic):
    mol, _ = synthetic
    report = RunConformerAnalysis(mol, columns=['xtb_homo_ev'],
                                  boltzmann_weighted=True,
                                  boltzmann_temperature=300.)
    weights = report['boltzmann']['weights']
    assert report['weighting'] == 'uniform_supplied_geometries'
    assert report['boltzmann']['enabled'] is True
    assert sum(weights) == pytest.approx(1.)
    assert weights[0] > weights[1] > weights[2]
    stats = report['boltzmann']['statistics']['xtb_homo_ev']
    assert stats['weighting'] == 'electronic_energy_not_free_energy'
    assert stats['mean'] < -2.


def test_failure_counts_and_fatal_errors(synthetic, monkeypatch):
    mol, calculate = synthetic
    def fail_middle(self, mol, **kwargs):
        if kwargs.get('conf_id', self.options.conf_id) == 8: raise XTBRunError('SCF failed')
        return calculate(self, mol, **kwargs)
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', fail_middle)
    report = RunConformerAnalysis(mol, columns=['xtb_homo_ev'])
    assert report['n_failed'] == 1 and report['status'] == 'partial'
    assert report['statistics']['xtb_homo_ev']['n_valid'] == 2
    assert report['statistics']['xtb_homo_ev']['n_missing'] == 0
    assert report['reports'][1]['error']['message'] == 'SCF failed'
    with pytest.raises(XTBRunError): RunConformerAnalysis(mol, return_errors=False)
    for error in (XTBEnvironmentError('runtime'), RuntimeError('unexpected'), ValueError('analysis')):
        def fail(*a, **k): raise error
        monkeypatch.setattr(conformers.XTBCalculator, 'calculate', fail)
        with pytest.raises(type(error)): RunConformerAnalysis(mol)
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', lambda *a, **k: (_ for _ in ()).throw(XTBRunError('SCF')))
    report = RunConformerAnalysis(mol)
    assert report['status'] == 'failed' and report['n_failed'] == 3
    assert all(s['mean'] is None and s['n_failed'] == 3 for s in report['statistics'].values())


@pytest.mark.parametrize('kwargs', [dict(conf_ids=[]), dict(conf_ids=[3, 3]), dict(conf_ids=[True]),
    dict(conf_ids=[-1]), dict(conf_ids=[999]), dict(columns=['xtb_homo_idx']),
    dict(columns=['xtb_dip_x_ebohr']), dict(columns=[]), dict(degeneracy_ev=-1),
    dict(options=XTBOptions(conf_id=3)), dict(return_errors=1)])
def test_invalid_inputs_do_not_calculate(synthetic, monkeypatch, kwargs):
    mol, _ = synthetic
    def unexpected(*a, **k): raise AssertionError('must validate first')
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', unexpected)
    with pytest.raises(ValueError): RunConformerAnalysis(mol, **kwargs)


def test_bad_later_geometry_rejected_before_first_calculation(synthetic, monkeypatch):
    mol, _ = synthetic; mol.GetConformer(12).SetAtomPosition(0, (float('nan'), 0, 0))
    def unexpected(*a, **k): raise AssertionError('must validate all first')
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', unexpected)
    with pytest.raises(ValueError, match='finite'): RunConformerAnalysis(mol)


@pytest.mark.parametrize('field,value', [('engine', {'sha256': 'other'}), ('charge', 2),
                                        ('conformer_id', 3), ('CalculationSettings', {})])
def test_mixed_or_wrong_provenance_is_rejected(synthetic, monkeypatch, field, value):
    mol, calculate = synthetic
    def mixed(self, mol, **kwargs):
        r = calculate(self, mol, **kwargs)
        identifier = kwargs.get('conf_id', self.options.conf_id)
        return replace(r, metadata={**r.metadata, field: value}) if identifier == 8 else r
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', mixed)
    with pytest.raises(ValueError, match='provenance|protocol'): RunConformerAnalysis(mol)


def test_cache_errors_are_fatal(synthetic, monkeypatch):
    mol, _ = synthetic
    monkeypatch.setattr(conformers, '_key', lambda *a: 'key')
    class Broken:
        def get(self, key): raise XTBCacheError('corrupt')
    with pytest.raises(XTBCacheError): RunConformerAnalysis(mol, cache=Broken())


def test_missing_gfn0_parameters_are_environment_errors(monkeypatch):
    from pathlib import Path
    from eled.runtime import _parameter_file
    monkeypatch.setattr(Path, 'is_file', lambda self: False)
    with pytest.raises(XTBEnvironmentError, match='param_gfn0'):
        _parameter_file('GFN0-xTB', '/missing/lib/libxtb.so')


def test_charged_dipoles_and_cartesian_fields_excluded():
    assert 'xtb_dip_debye' in ConformerDescriptorColumns()
    assert 'xtb_dip_debye' not in ConformerDescriptorColumns(charge=1)
    assert not any('_centroid_x_' in k or k.endswith('_index') for k in ConformerDescriptorColumns())


@pytest.mark.integration
@pytest.mark.parametrize('smiles', ['O', '[NH4+]'])
def test_real_engine_rigid_transform_and_disk_cache(smiles, tmp_path, monkeypatch):
    if os.environ.get('ELED_EXPECT_BASIS') != '1': pytest.skip('extended libxtb required')
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles)); assert AllChem.EmbedMolecule(mol, randomSeed=49) == 0
    conf = Chem.Conformer(mol.GetConformer()); conf.SetId(7)
    xyz = conf.GetPositions(); transformed = xyz @ np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]]) + [3., -2., 5.]
    for i, pos in enumerate(transformed): conf.SetAtomPosition(i, pos)
    mol.AddConformer(conf, assignId=False)
    cache = XTBDiskCache(tmp_path/'conformers.sqlite')
    report = RunConformerAnalysis(mol, cache=cache, frontier_depth=0)
    assert report['status'] == 'complete' and len(cache) == 2
    for key, stats in report['statistics'].items(): assert stats['range'] < 1e-6, key
    def unexpected(*a, **k): raise AssertionError('cache should hit')
    monkeypatch.setattr(conformers.XTBCalculator, 'calculate', unexpected)
    assert RunConformerAnalysis(mol, cache=cache, frontier_depth=0) == report


def test_local_statistics_keep_entity_frontier_coverage_and_units(synthetic,monkeypatch):
    mol,_=synthetic
    original=conformers.CalcDescriptorReport
    def report(*args,**kwargs):
        out=original(*args,**kwargs)
        i=[3,8,12].index(int(out['record_id']))
        out['frontier_bonds']={'bonds':[
            {'bond_index':0,'frontier_group':'homo','orb_idx':[1],'mul_ovlp_e':float(i+1)},
            {'bond_index':0,'frontier_group':'lumo','orb_idx':[2],'mul_ovlp_e':float(-i-1)}],
            'units':{'mul_ovlp_e':'e'}}
        out['pi']={'atoms':[{'atom_index':0,'pi_n_e':float(i+1),'pi_q_e':None,'normal':None}],
                   'bonds':[{'bond_index':0,'pi_bo':(.2,.4,None)[i]}],
                   'units':{'pi_n_e':'e','pi_q_e':'e','pi_bo':'dimensionless'}}
        out['frontier']={'homo':{'orb_idx':[1] if i==0 else [0,1]},'lumo':{'orb_idx':[2]}}
        return out
    monkeypatch.setattr(conformers,'CalcDescriptorReport',report)
    out=RunConformerAnalysis(mol,columns=['xtb_homo_ev'],boltzmann_weighted=True)
    bonds=out['bond_statistics']['0']; atoms=out['atom_statistics']['0']
    assert bonds['frontier_bonds.homo.mul_ovlp_e']['mean']==2.
    assert bonds['frontier_bonds.lumo.mul_ovlp_e']['mean']==-2.
    assert bonds['frontier_bonds.homo.mul_ovlp_e']['n_valid']==3
    assert bonds['pi.pi_bo']['n_valid']==2 and bonds['pi.pi_bo']['n_missing']==1
    assert bonds['pi.pi_bo']['mean']==pytest.approx(.3)
    assert 'pi.pi_n_e' not in bonds and 'pi.pi_bo' not in atoms
    assert atoms['pi.pi_q_e']['n_valid']==0 and atoms['pi.pi_q_e']['mean'] is None
    assert atoms['pi.pi_n_e']['mean']==2.
    assert out['frontier_group_composition_changed']
    assert out['atom_units']['pi.pi_n_e']=='e'
    assert out['bond_units']['frontier_bonds.lumo.mul_ovlp_e']=='e'
    weights=np.asarray(out['boltzmann']['weights'])
    assert out['boltzmann']['atom_statistics']['0']['pi.pi_n_e']['mean']==pytest.approx(weights@[1,2,3])
    weighted=out['boltzmann']['bond_statistics']['0']['pi.pi_bo']
    assert weighted['weight_sum']==pytest.approx(weights[:2].sum())
    assert weighted['mean']==pytest.approx(weights[:2]@[.2,.4]/weights[:2].sum())
    assert weighted['n_missing']==1
    assert json.loads(json.dumps(out,allow_nan=False))==out


def test_weighted_missing_failed_and_underflow_are_distinct():
    # A finite feature remains observed even when its energy weight underflows.
    out=conformers._summary([None,2.,None],[False,False,True],
                           weights=np.array([1.,0.,0.]),energy_valid=np.array([True,True,False]))
    assert out['n_valid']==1 and out['n_missing']==1 and out['n_failed']==1
    assert out['weight_sum']==0 and out['mean'] is None
    assert out['missing_reason']=='no_positive_weight_for_available_values'


def test_sensitivity_requires_matched_per_feature_conformer_coverage(synthetic):
    from copy import deepcopy
    from eled import DescriptorSensitivity
    mol,_=synthetic
    a=RunConformerAnalysis(mol,columns=['xtb_homo_ev'])
    b=deepcopy(a); b['calculation_provenance']['method']='GFN1-xTB'
    baseline=DescriptorSensitivity(a,gfn1=b,gfn2=a)
    assert baseline['values']['xtb_homo_ev']['gfn12_delta']==0.
    assert baseline['atoms']['0']['atoms.q_e']['gfn12_delta']==0.
    # Same number of valid points, different geometries: comparing means is invalid.
    a['reports'][0]['atoms'][0]['q_e']=None
    b['reports'][1]['atoms'][0]['q_e']=None
    output=DescriptorSensitivity(a,gfn1=b,gfn2=a)['atoms']['0']['atoms.q_e']
    assert output['gfn12_delta'] is None
    assert output['missing_reasons']['gfn1_vs_gfn2']=='unmatched_conformer_coverage'
