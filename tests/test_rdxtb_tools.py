"""Exercise the compiled ROMol binding against independent xTB entry points."""
from helpers import calculate_xtb

from concurrent.futures import ThreadPoolExecutor, TimeoutError
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem, rdEHTTools

from eled.runtime import _find_xtb_library, _extended_libxtb_build, _native_gfn0_build

if importlib.util.find_spec('eled.rdXTBTools') is None:
    if os.environ.get('ELED_REQUIRE_RDXTB') == '1':
        raise RuntimeError('The compiled rdXTBTools module is required for this run')
    pytest.skip('optional rdXTBTools extension is not built', allow_module_level=True)

# A present but broken extension must fail collection, never silently skip.
from eled import rdXTBTools

pytestmark = pytest.mark.integration


def molecule(smiles='CCO'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    mol.GetConformer().SetId(7)
    return mol


@pytest.fixture(scope='module')
def library():
    path = _find_xtb_library(None)
    if not _extended_libxtb_build(path):
        pytest.skip('reviewed extended libxtb is required')
    return path


@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['CCO', 'O', 'c1ccccc1', '[NH4+]', '[CH3]'])
def test_numerical_agreement_units_and_mapping(method, smiles, library):
    if method == 'GFN0-xTB' and not _native_gfn0_build(library):
        pytest.skip('native GFN0 requires a verified patched library')
    mol = molecule(smiles)
    original = mol.ToBinary()
    ok, result = rdXTBTools.RunMol(mol, confId=7, method=method, xtbLibrary=library,
                                 keepOrbitalCoefficients=True)
    assert ok and result.success and result.errorMessage == ''
    native = calculate_xtb(mol, conf_id=7, method=method, backend='native', xtb_library=library)
    assert result.totalEnergyHartree == pytest.approx(native.total_energy_hartree, abs=1e-10)
    np.testing.assert_allclose(result.GetOrbitalEnergies(), native.orbital_energies_ev, atol=2e-5, rtol=0)
    np.testing.assert_allclose(result.GetOrbitalOccupancies(), native.orbital_occupations, atol=2e-5, rtol=0)
    np.testing.assert_allclose(result.GetAtomicCharges(), native.partial_charges, atol=2e-5, rtol=0)
    np.testing.assert_allclose(result.GetDipole(), native.dipole_e_bohr, atol=1e-9, rtol=0)
    np.testing.assert_allclose(result.GetWibergBondOrders(), native.wiberg_bond_orders, atol=1e-9, rtol=0)
    # Native basis tests independently check C.T S C and density populations.
    np.testing.assert_allclose(np.abs(result.GetOrbitalCoefficients()),
                               np.abs(native.orbital_coefficients), atol=1e-8, rtol=0)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        assert result.GetWibergBondOrders()[i, j] == pytest.approx(native.wiberg_bond_orders[i][j], abs=2e-5)
    assert result.numAtoms == mol.GetNumAtoms()
    assert result.numOrbitals == len(native.orbital_energies_ev)
    assert result.numElectrons == native.number_of_electrons
    assert result.GetOrbitalOccupancies().sum() == pytest.approx(result.numElectrons, abs=1e-7)
    assert result.numTotalElectrons == sum(a.GetAtomicNum() for a in mol.GetAtoms()) - Chem.GetFormalCharge(mol)
    assert result.charge == Chem.GetFormalCharge(mol)
    assert result.numUnpairedElectrons == sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    assert result.GetAtomicCharges().sum() == pytest.approx(result.charge, abs=1e-7)
    assert result.conformerId == 7 and result.method == method
    assert result.apiVersion in (10000, 20000)
    assert result.librarySHA256 == hashlib.sha256(Path(library).read_bytes()).hexdigest()
    assert result.parameterSHA256 == hashlib.sha256(Path(result.parameterFile).read_bytes()).hexdigest()
    assert mol.ToBinary() == original


def test_eht_coexistence_and_result_ownership(library):
    mol = molecule('O')
    eht_ok, eht = rdEHTTools.RunMol(mol, confId=7)
    ok, xtb = rdXTBTools.RunMol(mol, confId=7, xtbLibrary=library)
    assert eht_ok and ok
    assert eht.numElectrons == xtb.numElectrons == 8
    assert xtb.numTotalElectrons == 10
    assert rdXTBTools.totalEnergyUnit == 'Hartree'
    assert rdXTBTools.orbitalEnergyUnit == 'eV'
    assert rdXTBTools.dipoleUnit == 'e*Bohr'
    assert rdXTBTools.coefficientAxes == ('AO', 'MO')
    assert xtb.GetGradient().shape == (mol.GetNumAtoms(), 3)
    assert hasattr(eht, 'totalEnergy') and not hasattr(xtb, 'totalEnergy')
    assert not hasattr(xtb, 'fermiEnergy')
    assert not hasattr(xtb, 'GetReducedOverlapPopulationMatrix')
    for getter in [xtb.GetAtomicCharges, xtb.GetOrbitalEnergies, xtb.GetOrbitalOccupancies,
                   xtb.GetWibergBondOrders, xtb.GetDipole]:
        value = getter()
        assert value.flags.owndata
        original = value.copy()
        value[:] = 999
        np.testing.assert_array_equal(getter(), original)
    with pytest.raises(AttributeError):
        xtb.totalEnergyHartree = 0
    with pytest.raises(ValueError, match='not retained'):
        xtb.GetOrbitalCoefficients()
    retained = xtb.GetAtomicCharges()
    del xtb
    assert np.isfinite(retained).all()


@pytest.mark.parametrize('change, options, message', [
    ('empty', {}, 'contain atoms'),
    ('no-conformer', {}, '3D conformer'),
    ('implicit-h', {}, 'Chem.AddHs'),
    ('h-count', {}, 'Chem.AddHs'),
    ('2d', {}, 'marked as 3D'),
    ('nonfinite', {}, 'finite'),
    ('dummy', {}, 'Dummy'),
    (None, {'confId': 0}, 'Conformer ID'),
    (None, {'confId': -2}, 'confId'),
    (None, {'method': 'GFN-FF'}, 'method'),
    (None, {'numUnpairedElectrons': 1}, 'inconsistent'),
    (None, {'charge': 0.5}, 'integer'),
    (None, {'charge': True}, 'integer'),
    (None, {'maxIterations': 0}, 'iteration'),
    (None, {'accuracy': float('nan')}, 'accuracy'),
    (None, {'electronicTemperature': -1}, 'temperature'),
])
def test_invalid_inputs_fail_before_loading_xtb(change, options, message):
    mol = molecule('O')
    if change == 'empty':
        mol = Chem.Mol()
    elif change == 'no-conformer':
        mol.RemoveAllConformers()
    elif change in ('implicit-h', 'h-count'):
        mol = Chem.MolFromSmiles('O' if change == 'implicit-h' else '[NH4+]')
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()))
    elif change == '2d':
        mol.GetConformer().Set3D(False)
    elif change == 'nonfinite':
        mol.GetConformer().SetAtomPosition(0, (float('nan'), 0, 0))
    elif change == 'dummy':
        mol.GetAtomWithIdx(0).SetAtomicNum(0)
    with pytest.raises(ValueError, match=message):
        rdXTBTools.RunMol(mol, xtbLibrary='/must/not/be/loaded', **options)


def test_explicit_charge_spin_and_library_errors(library):
    mol = molecule('O')
    ok, result = rdXTBTools.RunMol(mol, charge=np.int64(1), numUnpairedElectrons=1, xtbLibrary=library)
    assert ok and result.charge == 1 and result.numUnpairedElectrons == 1
    assert result.numTotalElectrons == 9 and result.numElectrons == 7
    assert result.GetAtomicCharges().sum() == pytest.approx(1, abs=1e-7)
    with pytest.raises(RuntimeError, match='Cannot load xTB'):
        rdXTBTools.RunMol(mol, xtbLibrary='/does/not/exist')
    with pytest.raises(ValueError, match='Parameter file not found'):
        rdXTBTools.RunMol(mol, xtbLibrary=library, parameterDirectory='/does/not/exist')


def test_scf_failure_is_inspectable_and_recovers(library):
    mol = molecule('O')
    ok, before = rdXTBTools.RunMol(mol, xtbLibrary=library)
    assert ok
    ok, failed = rdXTBTools.RunMol(mol, xtbLibrary=library, maxIterations=1)
    assert not ok and not failed.success and failed.errorMessage
    assert failed.method == 'GFN2-xTB' and failed.conformerId == 7
    for name in ['totalEnergyHartree', 'GetAtomicCharges', 'GetOrbitalEnergies',
                 'GetOrbitalOccupancies', 'GetWibergBondOrders', 'GetDipole', 'GetOrbitalCoefficients']:
        with pytest.raises(RuntimeError, match='No numerical results'):
            value = getattr(failed, name)
            if callable(value):
                value()
    ok, after = rdXTBTools.RunMol(mol, xtbLibrary=library)
    assert ok and after.totalEnergyHartree == pytest.approx(before.totalEnergyHartree, abs=1e-10)


def test_native_calculation_does_not_delegate_to_python_backends(library, monkeypatch):
    import eled.backends.native as native
    import eled.xtb as adapter
    def forbidden(*args, **kwargs):
        pytest.fail('compiled calculations must call libxtb directly')
    monkeypatch.setattr(native._NativeXTB, '__init__', forbidden)
    monkeypatch.setattr(native, 'run_xtb_native', forbidden)
    monkeypatch.setattr(adapter, '_run_xtb', forbidden)
    assert rdXTBTools.RunMol(molecule(), xtbLibrary=library)[0]


def test_mixed_threads_share_lock_and_preserve_method_state(library):
    from eled.backends.native import _native_lock
    mol = molecule()
    ok, baseline = rdXTBTools.RunMol(mol, xtbLibrary=library)
    assert ok
    with ThreadPoolExecutor(max_workers=4) as pool:
        with _native_lock:
            pending = pool.submit(rdXTBTools.RunMol, mol, xtbLibrary=library)
            with pytest.raises(TimeoutError):
                pending.result(timeout=.1)
        assert pending.result(timeout=30)[0]
        methods = ['gfn2', 'gfn1', 'gfn0' if _native_gfn0_build(library) else 'gfn2'] * 3
        def calculate(index):
            method = methods[index]
            if index % 2:
                return calculate_xtb(mol, method=method, backend='native', xtb_library=library).total_energy_hartree
            ok, result = rdXTBTools.RunMol(mol, method=method, xtbLibrary=library)
            assert ok
            return result.totalEnergyHartree
        values = list(pool.map(calculate, range(len(methods))))
    for value, method in zip(values, methods):
        if method == 'gfn2':
            assert value == pytest.approx(baseline.totalEnergyHartree, abs=1e-10)
    if not _native_gfn0_build(library):
        with pytest.raises(ValueError, match='verified patched'):
            rdXTBTools.RunMol(mol, method='gfn0', xtbLibrary=library)
    ok, final = rdXTBTools.RunMol(mol, xtbLibrary=library)
    assert ok
    np.testing.assert_allclose(final.GetOrbitalEnergies(), baseline.GetOrbitalEnergies(), atol=1e-8, rtol=0)


def test_gfn0_verification_rejects_missing_or_mismatched_manifest(library, tmp_path):
    evidence = _native_gfn0_build(library)
    if not evidence:
        pytest.skip('requires patched library provenance to exercise manifest changes')
    copied = tmp_path / 'lib/libxtb.so'
    copied.parent.mkdir()
    shutil.copy2(library, copied)
    manifest = tmp_path / 'share/xtb/eled/build.json'
    manifest.parent.mkdir(parents=True)
    record = json.loads(Path(evidence['manifest']).read_text())
    mol = molecule('O')
    for damage in ['missing', 'json', 'capability', 'patch', 'hash']:
        damaged = json.loads(json.dumps(record))
        if damage == 'capability':
            damaged['capabilities']['native_gfn0'] = False
        elif damage == 'patch':
            if 'gfn0_fix' in damaged['source']:
                damaged['source']['gfn0_fix']['file_sha256'] = 'different-source'
            else:
                damaged['source']['patch_sha256'] = 'different-patch'
        elif damage == 'hash':
            damaged['files']['lib/libxtb.so'] = 'different-binary'
        if damage == 'json':
            manifest.write_text('{')
        elif damage != 'missing':
            manifest.write_text(json.dumps(damaged))
        with pytest.raises(ValueError, match='verified patched'):
            rdXTBTools.RunMol(mol, method='gfn0', xtbLibrary=str(copied))


def test_common_descriptors_accept_compiled_results(library):
    from eled import (calc_xtb_descriptors, calc_xtb_atom_descriptors,
                           calc_xtb_fragment_descriptors, calc_orbital_populations,
                           DescriptorUnavailableError)
    mol = molecule()
    ok, compiled = rdXTBTools.RunMol(mol, xtbLibrary=library)
    assert ok
    native = calculate_xtb(mol, xtb_library=library)
    a, b = calc_xtb_descriptors(compiled), calc_xtb_descriptors(native)
    for key in a:
        if key in ('xtb_backend', 'xtb_version'): continue
        if a[key] is None or isinstance(a[key], str): assert a[key] == b[key]
        else: assert a[key] == pytest.approx(b[key], abs=1e-8)
    assert len(calc_xtb_atom_descriptors(compiled)) == mol.GetNumAtoms()
    fragments = {'heavy': [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1],
                 'hydrogen': [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]}
    assert calc_xtb_fragment_descriptors(compiled, fragments)['pairs'][0]['wiberg_sum'] == pytest.approx(
        calc_xtb_fragment_descriptors(native, fragments)['pairs'][0]['wiberg_sum'], abs=1e-8)
    with pytest.raises(DescriptorUnavailableError): calc_orbital_populations(compiled)
    ok, failed = rdXTBTools.RunMol(mol, xtbLibrary=library, maxIterations=1)
    assert not ok
    with pytest.raises(ValueError, match='failed'): calc_xtb_descriptors(failed)
