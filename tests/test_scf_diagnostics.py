"""SCC diagnostic semantics, invalidation, copies and binding agreement."""
from helpers import calculate_xtb
import ctypes
import json
from pathlib import Path

import pytest
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import calc_descriptor_report
from rdkit_eth.runtime import _find_xtb_library, _extended_libxtb_build, _SCF_DIAGNOSTICS_PATCH
from rdkit_eth.scf import scf_record
from rdkit_eth.errors import XTBEnvironmentError, XTBNativeNotFoundError


@pytest.fixture
def library():
    try: path = _find_xtb_library()
    except XTBNativeNotFoundError: pytest.skip('libxtb unavailable')
    if not _extended_libxtb_build(path): pytest.skip('reviewed extended libxtb unavailable')
    return path


def mol():
    molecule = Chem.AddHs(Chem.MolFromSmiles('O'))
    assert AllChem.EmbedMolecule(molecule, randomSeed=17) == 0
    return molecule


def test_reviewed_patch_and_missing_semantics():
    import hashlib
    root = Path(__file__).resolve().parents[1]
    lock = json.loads((root/'tools/xtb-source.json').read_text())['diagnostics_extension']
    assert hashlib.sha256((root/'tools'/lock['patch']).read_bytes()).hexdigest() == lock['patch_sha256'] == _SCF_DIAGNOSTICS_PATCH
    assert _SCF_DIAGNOSTICS_PATCH in (root/'cpp/rdXTBTools/XTBTools.cpp').read_text()
    assert scf_record('GFN0-xTB')['reason'] == 'not_applicable_non_scf_method'
    assert not scf_record('GFN2-xTB')['available']
    with pytest.raises(XTBEnvironmentError): scf_record('GFN2-xTB', 3, [float('nan'), 0., 1., 1.])


@pytest.mark.integration
@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
def test_success_diagnostics_and_compiled_parity(library, method):
    molecule = mol()
    result = calculate_xtb(molecule, method=method, xtb_library=library)
    diagnostic = result.metadata['scf']
    if method == 'GFN0-xTB': assert diagnostic['available'] is False
    else:
        assert 1 <= diagnostic['iterations'] <= 250 and diagnostic['criteria_satisfied']
        assert diagnostic['energy_change_hartree'] < diagnostic['energy_threshold_hartree']
        assert calc_descriptor_report(molecule, result, include_orbitals=False)['quality']['scf_iterations'] == diagnostic['iterations']
    try: from rdkit_eth import rdXTBTools
    except ImportError: return  # native coverage still runs without the optional extension
    if not hasattr(rdXTBTools.XTBResults, 'scfIterations'): return
    from rdkit_eth.electronic_descriptors import _as_result
    ok, compiled = rdXTBTools.RunMol(molecule, method=method, xtbLibrary=library)
    assert ok
    compiled_diagnostic = _as_result(compiled).metadata['scf']
    assert compiled_diagnostic['available'] == diagnostic['available']
    if diagnostic['available']:
        assert compiled_diagnostic['iterations'] == diagnostic['iterations']
        assert compiled_diagnostic['criteria_satisfied'] == diagnostic['criteria_satisfied']
        for key in ('energy_change_hartree', 'mixer_rms',
                    'energy_threshold_hartree', 'mixer_threshold'):
            # ctypes and C++ bindings invoke the same engine through separate
            # call sequences; final SCF diagnostics can differ by roundoff.
            assert compiled_diagnostic[key] == pytest.approx(diagnostic[key], rel=1e-4, abs=1e-12)


@pytest.mark.integration
def test_c_api_copy_failed_reuse_and_gfn0_invalidate(library):
    from rdkit_eth.backends.native import _NativeXTB
    api = _NativeXTB(library)
    h = ctypes.c_void_p
    getter = api._bind('xtb_rdkit_getSCFDiagnostics', ctypes.c_int,
                      [h, h, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)])
    env = api.new_environment(); api.set_verbosity(env, 0)
    n = ctypes.c_int(3); charge = ctypes.c_double(0.); spin = ctypes.c_int(0)
    numbers = (ctypes.c_int*3)(8, 1, 1)
    xyz = (ctypes.c_double*9)(0, 0, 0, 0, 0, 1.8, 1.7, 0, -.6)
    molecule = api.new_molecule(env, ctypes.byref(n), numbers, xyz, ctypes.byref(charge), ctypes.byref(spin), None, None)
    calc = api.new_calculator(); result = api.new_results(); copied = None
    def read(handle):
        count = ctypes.c_int(123456); data = (ctypes.c_double*4)(*([987.]*4))
        available = getter(env, handle, ctypes.byref(count), data)
        return available, count.value, list(data)
    try:
        assert read(result) == (0, 123456, [987.]*4)
        api.load_gfn2xtb(env, molecule, calc, None)
        api.singlepoint(env, molecule, calc, result); api.check(env, 'test singlepoint')
        original = read(result); assert original[0] == 1
        copied = api.copy_results(result); assert read(copied) == original
        api.set_max_iterations(env, calc, 1)
        api.singlepoint(env, molecule, calc, result)
        assert api.check_environment(env)
        assert read(result) == (0, 123456, [987.]*4)
        assert read(copied) == original
        # Consume and clear the environment error, then recover with non-SCF GFN0.
        with pytest.raises(Exception): api.check(env, 'expected iteration-limit failure')
        parameters = str(Path(library).resolve().parent.parent/'share/xtb/param_gfn0-xtb.txt').encode()
        api.load_gfn0xtb(env, molecule, calc, parameters)
        api.singlepoint(env, molecule, calc, result); api.check(env, 'GFN0 recovery')
        assert read(result) == (0, 123456, [987.]*4)
    finally:
        for destructor, handle in [(api.del_results, copied), (api.del_results, result),
                (api.del_calculator, calc), (api.del_molecule, molecule), (api.del_environment, env)]:
            api._delete(destructor, handle)
