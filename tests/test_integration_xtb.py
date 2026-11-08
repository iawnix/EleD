"""Real shared-library calculations, with no xTB executable dependency."""
from helpers import calculate_xtb
import ctypes
import os
import subprocess
from unittest.mock import patch

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import AnnotateMol, RunXTBBatch, XTBCache
from eled.runtime import _find_xtb_library, _extended_libxtb_build, _native_gfn0_build
from eled.errors import XTBNativeNotFoundError

pytestmark = pytest.mark.integration


@pytest.fixture(scope='module')
def library():
    try:
        path = _find_xtb_library(None)
    except XTBNativeNotFoundError:
        if os.environ.get('XTB_LIBRARY'):
            raise
        pytest.skip('extended libxtb is not installed')
    if not _extended_libxtb_build(path):
        pytest.skip('reviewed extended libxtb is required')
    return path


@pytest.fixture
def no_executable(monkeypatch):
    monkeypatch.setenv('PATH', '')
    monkeypatch.setenv('XTB_PATH', '/removed/executable')
    def forbidden(*args, **kwargs):
        pytest.fail('calculations must not launch an external process')
    monkeypatch.setattr(subprocess, 'Popen', forbidden)


def molecule(smiles):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    return mol


@pytest.mark.parametrize('method', ['GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['CO', 'c1ccccc1', '[NH4+]', '[CH3]'])
def test_native_results_and_mapping_without_executable(smiles, method, library, no_executable):
    mol = molecule(smiles)
    original = mol.ToBinary()
    native = calculate_xtb(mol, method=method, xtb_library=library)
    assert native.backend == 'native'
    assert np.isfinite(native.total_energy_hartree) and native.total_energy_hartree < 0
    assert sum(native.partial_charges) == pytest.approx(Chem.GetFormalCharge(mol), abs=1e-6)
    assert sum(native.orbital_occupations) == pytest.approx(native.number_of_electrons, abs=1e-7)
    assert np.isfinite(native.orbital_energies_ev).all()
    nao = native.number_of_molecular_orbitals
    assert np.shape(native.orbital_coefficients) == (nao, nao)
    np.testing.assert_allclose(native.wiberg_bond_orders, np.array(native.wiberg_bond_orders).T, atol=1e-10)
    target = AnnotateMol(mol, native)
    for bond in mol.GetBonds():
        i, j = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        assert target.GetBondWithIdx(bond.GetIdx()).GetDoubleProp('xtb_wiberg_order') == native.wiberg_bond_orders[i][j]
    assert native.metadata['engine']['sha256']
    api_version = ctypes.CDLL(library).xtb_getAPIVersion
    api_version.argtypes = []
    api_version.restype = ctypes.c_int
    assert native.metadata['native_api_version'] == api_version()
    assert native.metadata['native_api_version'] // 10000 in (1, 2)
    assert mol.ToBinary() == original


def test_fixed_water_matches_archived_reference(library, no_executable):
    # Frozen geometry and energy from the 0.1.0 official-xTB CLI installation
    # check; the reference does not require executing that program at test time.
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    conf = Chem.Conformer(3)
    for i, xyz in enumerate(((0., 0., 0.), (.9572, 0., 0.), (-.239987, .927297, 0.))):
        conf.SetAtomPosition(i, xyz)
    conf.Set3D(True)
    mol.AddConformer(conf)
    result = calculate_xtb(mol, xtb_library=library)
    assert result.total_energy_hartree == pytest.approx(-5.07037168, abs=2e-7)


def test_gfn0_selection_does_not_contaminate_native_gfn2(library, no_executable):
    from eled.backends.native import run_xtb_native
    mol = molecule('CCO')
    verified = _native_gfn0_build(library)
    expected = os.environ.get('ELED_EXPECT_NATIVE_GFN0')
    if expected is not None:
        assert bool(verified) == (expected == '1')
    before = calculate_xtb(mol, xtb_library=library)
    if verified:
        gfn0 = calculate_xtb(mol, method='GFN0-xTB', xtb_library=library)
        assert gfn0.backend == 'native' and gfn0.method == 'GFN0-xTB'
        assert gfn0.orbital_energies_ev
        assert gfn0.metadata['engine']['native_gfn0']['library_sha256'] == gfn0.metadata['engine']['sha256']
        calculate_xtb(molecule('CO'), method='gfn1', xtb_library=library)
        run_xtb_native(molecule('O'), method='gfn0', xtb_library=library)
    else:
        for entry in (run_xtb_native):
            with pytest.raises(ValueError, match='verified patched'):
                entry(mol, method='GFN0-xTB', xtb_library=library)
    after = calculate_xtb(mol, xtb_library=library)
    assert after.total_energy_hartree == pytest.approx(before.total_energy_hartree, abs=1e-10)
    np.testing.assert_allclose(after.orbital_energies_ev, before.orbital_energies_ev, atol=1e-8, rtol=0)
    np.testing.assert_allclose(after.orbital_coefficients, before.orbital_coefficients, atol=1e-8, rtol=0)


@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN2-xTB'])
def test_native_batch_and_cache_without_executable(method, library, no_executable):
    if method == 'GFN0-xTB' and not _native_gfn0_build(library):
        pytest.skip('native GFN0 needs a verified patched installation')
    from eled.backends.native import run_xtb_native
    mol = molecule('CCO')
    cache = XTBCache()
    with patch('eled.xtb.run_xtb_native', wraps=run_xtb_native) as calculate:
        results = RunXTBBatch([mol, mol], method=method, cache=cache, xtb_library=library)
    assert calculate.call_count == 1
    assert len(cache) == 1
    assert results[0] == results[1]
    assert results[0].backend == 'native'


def test_missing_library_fails_without_fallback(no_executable):
    with pytest.raises(XTBNativeNotFoundError):
        calculate_xtb(molecule('O'), xtb_library='/does/not/exist')
