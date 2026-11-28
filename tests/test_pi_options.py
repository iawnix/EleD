"""Pi analysis settings propagate through public workflows and JSON reports."""
from dataclasses import replace
import json

import pytest
from rdkit import Chem

from eled import (CalcDescriptorReport, CalcDensityDescriptors, CalcPiDescriptors,
                  RunConformerAnalysis, RunDescriptorBatch, XTBRunError)
from eled.input import _input_metadata, _xyz_from_mol
from test_descriptor_regressions import pi_system


def asymmetric_pi_system():
    mol, result = pi_system()
    occ = list(result.orbitalOccupancies)
    occ[2], occ[5] = 2., 0.
    return mol, replace(result, orbitalOccupancies=tuple(occ))


@pytest.mark.parametrize('scheme', ['mulliken', 'lowdin'])
def test_pi_options_match_standalone_and_replay_after_json(scheme):
    mol, result = asymmetric_pi_system()
    refs = {0: .5}
    options = dict(piScheme=scheme, piPlaneTolAngstrom=.02, piRef=refs)
    direct = CalcPiDescriptors(mol, result, scheme=scheme, planeTolAngstrom=.02, piRef=refs)
    report = CalcDescriptorReport(mol, result, **options)
    assert report['pi'] == direct == CalcDensityDescriptors(mol, result, **options)['pi']
    assert report['pi']['pi_ref'] == {'0': .5}
    assert report['pi']['atoms'][0]['pi_q_e'] == pytest.approx(.5-direct['atoms'][0]['pi_n_e'])
    opposite = CalcPiDescriptors(mol, result, scheme='lowdin' if scheme == 'mulliken' else 'mulliken')
    assert direct['atoms'][0]['pi_n_e'] != pytest.approx(opposite['atoms'][0]['pi_n_e'])
    loaded = json.loads(json.dumps(report, allow_nan=False))
    assert loaded == report
    assert CalcDescriptorReport(mol, result, **loaded['analysis_settings']) == report
    refs[0] = 2.
    assert report['analysis_settings']['piRef'] == {'0': .5}


def test_report_plane_tolerance_changes_projection_availability():
    mol, result = asymmetric_pi_system()
    mol.GetConformer().SetAtomPosition(0, (1.4, 0., .05))
    result = replace(result, metadata={**result.metadata,
        **_input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))})
    tight = CalcDescriptorReport(mol, result, piPlaneTolAngstrom=1e-6)
    loose = CalcDescriptorReport(mol, result, piPlaneTolAngstrom=.1)
    assert tight['pi']['atoms'][0]['missing']['pi_n_e'] == 'local_plane_nonplanar'
    assert loose['pi']['atoms'][0]['pi_n_e'] is not None


def test_batch_preserves_pi_settings_in_success_and_failure_rows(monkeypatch):
    mol, result = asymmetric_pi_system()
    monkeypatch.setattr('eled.batch.RunXTBBatch', lambda *a, **k: [result, XTBRunError('SCF'), XTBRunError('SCF')])
    reports = RunDescriptorBatch([mol]*3, piScheme='mulliken', piPlaneTolAngstrom=.02, piRef={0: .5})
    assert reports[0]['pi'] == CalcPiDescriptors(mol, result, scheme='mulliken', planeTolAngstrom=.02, piRef={0: .5})
    assert reports[0]['analysis_settings'] == reports[1]['analysis_settings'] == reports[2]['analysis_settings']
    reports[1]['analysis_settings']['piRef']['0'] = 2.
    assert reports[2]['analysis_settings']['piRef'] == {'0': .5}


def test_conformer_pi_statistics_and_settings_agree_with_direct_analysis(monkeypatch):
    mol, result = asymmetric_pi_system()
    second = Chem.Conformer(mol.GetConformer()); second.SetId(7)
    mol.AddConformer(second, assignId=False)
    monkeypatch.setattr('eled.batch._resolve_options', lambda options: options)
    def calculate(self, snapshot):
        conf_id = self.GetOptions().confId
        return replace(result, metadata={**result.metadata,
            **_input_metadata(snapshot, conf_id, 0, 0, _xyz_from_mol(snapshot, conf_id))})
    monkeypatch.setattr('eled.conformers.XTBCalculator.Calculate', calculate)
    out = RunConformerAnalysis(mol, piScheme='mulliken', piPlaneTolAngstrom=.02,
                               piRef={0: .5}, boltzmannWeighted=True)
    direct = CalcPiDescriptors(mol, result, scheme='mulliken', planeTolAngstrom=.02, piRef={0: .5})
    assert all(row['pi'] == direct for row in out['reports'])
    stats = out['atom_statistics']['0']['pi.pi_q_e']
    assert stats['mean'] == pytest.approx(direct['atoms'][0]['pi_q_e'])
    assert stats['n_valid'] == 2 and stats['std'] == 0.
    assert out['boltzmann']['atom_statistics']['0']['pi.pi_q_e']['mean'] == pytest.approx(stats['mean'])
    assert out['analysis_settings'] == out['reports'][0]['analysis_settings']
    assert json.loads(json.dumps(out, allow_nan=False)) == out


@pytest.mark.parametrize('options', [
    {'piScheme': 'invalid'}, {'piPlaneTolAngstrom': -1}, {'piPlaneTolAngstrom': True},
    {'piPlaneTolAngstrom': float('nan')}, {'piRef': []}, {'piRef': {True: 1}},
    {'piRef': {-1: 1}}, {'piRef': {12: 1}}, {'piRef': {'01': 1}},
    {'piRef': {0: 1, '0': 1}}, {'piRef': {0: 3}}, {'piRef': {0: float('nan')}},
])
def test_invalid_pi_options_rejected_before_engine_calls(monkeypatch, options):
    mol, result = pi_system()
    def unexpected(*args, **kwargs): raise AssertionError('must validate before calculation')
    monkeypatch.setattr('eled.batch.RunXTBBatch', unexpected)
    monkeypatch.setattr('eled.conformers.XTBCalculator.Calculate', unexpected)
    for run in (lambda: CalcDescriptorReport(mol, result, **options),
                lambda: RunDescriptorBatch([mol], **options),
                lambda: RunConformerAnalysis(mol, **options)):
        with pytest.raises((ValueError, TypeError)): run()


@pytest.mark.integration
@pytest.mark.parametrize('method', ['gfn0', 'gfn1', 'gfn2'])
def test_nondefault_pi_options_match_compiled_native_and_batch(method):
    import os
    if os.environ.get('ELED_EXPECT_BASIS') != '1':
        pytest.skip('set ELED_EXPECT_BASIS=1 for native pi integration')
    from rdkit.Chem import AllChem
    from eled import XTBCalculator, XTBOptions
    from eled.runtime import _find_xtb_library
    mod = pytest.importorskip('eled.rdXTBTools')
    mol = Chem.AddHs(Chem.MolFromSmiles('c1ccncc1'))
    assert AllChem.EmbedMolecule(mol, randomSeed=76) == 0
    lib = _find_xtb_library()
    options = dict(piScheme='mulliken', piPlaneTolAngstrom=.05, piRef={3: 1.})
    native = XTBCalculator(method=method, xtbLibrary=lib, keepOverlap=True).Calculate(mol)
    ok, compiled = mod.RunMol(mol, method=method, xtbLibrary=lib, keepOverlap=True)
    assert ok
    a, b = CalcDescriptorReport(mol, native, **options), CalcDescriptorReport(mol, compiled, **options)
    batch = RunDescriptorBatch([mol], options=XTBOptions(method=method, xtbLibrary=lib, keepOverlap=True), **options)[0]
    assert a['pi']['atoms'][3]['pi_q_e'] is not None
    for candidate in (b, batch):
        assert candidate['analysis_settings'] == a['analysis_settings']
        for entity, fields in [('atoms', ('pi_n_e', 'pi_q_e')), ('bonds', ('pi_bo',))]:
            for left, right in zip(a['pi'][entity], candidate['pi'][entity]):
                for field in fields:
                    if left[field] is None: assert right[field] is None
                    else: assert right[field] == pytest.approx(left[field], abs=1e-8)
