import pytest

from rdkit import Chem
from rdkit.Chem import AllChem

from rdkit_eth import optimize_xtb, prepare_conformers
from rdkit_eth.errors import XTBNativeNotFoundError
from rdkit_eth.runtime import _extended_libxtb_build, _find_xtb_library


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
        optimize_xtb(molecule('O'), max_iterations=1, fixed_atoms=[0, 1, 2])
    with pytest.raises(ValueError, match='scf_max_iterations'):
        optimize_xtb(molecule('O'), scf_max_iterations=0)


@pytest.mark.integration
def test_xtb_optimization_is_owned_and_reduces_energy(library):
    mol = molecule()
    before = mol.ToBinary()
    result = optimize_xtb(mol, xtb_library=library, max_iterations=100)
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
    result = prepare_conformers('CCO', optimize='xTB', xtb_library=library,
                                optimization_iterations=100,
                                attempts_per_strategy=1,
                                retry_random_coords=False)
    assert result.success
    candidate = result.report['attempts'][0]['candidates'][0]
    assert candidate['optimization_status'] == 0
    assert candidate['xtb_optimization']['converged']
