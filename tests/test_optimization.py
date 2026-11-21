import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from eled import OptimizeXTB, PrepareConformers
from eled.errors import XTBNativeNotFoundError, XTBRunError
from eled.runtime import _extended_libxtb_build, _find_xtb_library


def molecule(smiles='CCO'):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=11) == 0
    return mol


@pytest.fixture
def library():
    try:
        path = _find_xtb_library()
    except XTBNativeNotFoundError:
        pytest.skip('extended libxtb unavailable')
    if not _extended_libxtb_build(path):
        pytest.skip('reviewed extended libxtb required')
    return path


def test_optimization_validates_fixed_atoms():
    with pytest.raises(ValueError, match='every atom'):
        OptimizeXTB(molecule('O'), max_iterations=1, fixed_atoms=[0, 1, 2])
    with pytest.raises(ValueError, match='scf_max_iterations'):
        OptimizeXTB(molecule('O'), scf_max_iterations=0)


@pytest.mark.parametrize('mode', ['recover', 'all_fail', 'mixed', 'no_descent'])
def test_line_search_handles_failed_trial_calculations(monkeypatch, mode):
    from eled import optimization, XTBResult

    mol = molecule('O')
    original = mol.ToBinary()
    positions = []
    error = XTBRunError('trial SCF failed')

    def run(working, **options):
        positions.append(working.GetConformer().GetPositions())
        call = len(positions)
        if call > 1 and (mode == 'all_fail' or
                         (mode == 'recover' and call == 2) or
                         (mode == 'mixed' and call > 2)):
            raise error
        recovered = mode == 'recover' and call > 2
        return XTBResult(
            'GFN2-xTB', None, -1.1 if recovered else -1.0, None, 3, 8, 0,
            metadata={'gradient_hartree_per_bohr':
                      (0.0,) * 9 if recovered else (0.2, 0.0, 0.0) * 3},
        )

    monkeypatch.setattr(optimization, '_run_xtb_native', run)
    if mode == 'all_fail':
        with pytest.raises(XTBRunError) as caught:
            OptimizeXTB(mol)
        assert caught.value is error
        assert len(positions) == 15
    else:
        result = OptimizeXTB(mol)
        if mode == 'recover':
            import numpy as np
            assert result.success
            assert len(positions) == 3
            np.testing.assert_allclose(positions[2] - positions[0],
                                       (positions[1] - positions[0]) / 2)
        else:
            assert not result.success
            assert result.report['reason'] == 'line_search_failed'
            assert result.report['n_iterations'] == 0
            assert result.result.total_energy_hartree == -1.0
            assert len(positions) == 15
    assert mol.ToBinary() == original


@pytest.mark.integration
def test_xtb_optimization_is_owned_and_reduces_energy(library):
    mol = molecule()
    before = mol.ToBinary()
    result = OptimizeXTB(mol, xtb_library=library, max_iterations=100)
    assert result.success
    assert mol.ToBinary() == before
    report = result.report
    assert report['final_energy_hartree'] < report['initial_energy_hartree']
    assert report['final_max_gradient_hartree_per_bohr'] < report['gradient_tolerance_hartree_per_bohr']
    assert report['initial_max_gradient_hartree_per_bohr'] > 0.0
    assert report['scf_max_iterations'] == 250
    assert report['n_iterations'] <= report['max_iterations']


@pytest.mark.integration
def test_prepare_conformers_can_use_xtb_optimization(library):
    result = PrepareConformers('CCO', optimize='xTB', xtb_library=library,
                                optimization_iterations=100,
                                attempts_per_strategy=1,
                                retry_random_coords=False)
    assert result.success
    candidate = result.report['attempts'][0]['candidates'][0]
    assert candidate['optimization_status'] == 0
    assert candidate['xtb_optimization']['converged']
