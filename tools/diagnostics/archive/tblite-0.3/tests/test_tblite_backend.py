"""Optional tblite comparison backend.

These tests are skipped when the optional tblite wheel is absent. They never
make tblite part of the default eled installation.
"""

import importlib.util
import importlib.metadata

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import XTBRunError, XTBNativeNotFoundError, run_xtb, XTBCacheKey
from eled.runtime import _find_xtb_library

if importlib.util.find_spec("tblite") is None:
    pytest.skip("optional tblite backend is not installed", allow_module_level=True)

pytestmark = pytest.mark.integration


def molecule(smiles="CCO"):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    return mol


@pytest.fixture(scope="module")
def xtb_library():
    try:
        return _find_xtb_library(None)
    except XTBNativeNotFoundError:
        pytest.skip("libxtb is required for cross-engine comparisons")


@pytest.mark.parametrize("method", ["GFN1-xTB", "GFN2-xTB"])
@pytest.mark.parametrize("smiles", ["O", "CCO", "c1ccccc1", "[NH4+]"])
def test_tblite_matches_libxtb_for_closed_shell_molecules(method, smiles, xtb_library):
    mol = molecule(smiles)
    tblite = run_xtb(mol, method=method, backend="tblite")
    native = run_xtb(mol, method=method, backend="native", xtb_library=xtb_library)
    assert tblite.backend == "tblite"
    assert tblite.version == importlib.metadata.version("tblite")
    assert tblite.total_energy_eh == pytest.approx(native.total_energy_eh, abs=2e-7)
    # tblite and libxtb are separate implementations. Small SCF differences
    # accumulate in orbital energies and population analyses, so these are
    # agreement tolerances rather than bit-for-bit identity requirements.
    np.testing.assert_allclose(tblite.orbital_energies_ev, native.orbital_energies_ev,
                               atol=1.5e-4, rtol=0)
    np.testing.assert_allclose(tblite.orbital_occupations, native.orbital_occupations,
                               atol=1e-10, rtol=0)
    np.testing.assert_allclose(tblite.partial_charges, native.partial_charges,
                               atol=5e-5, rtol=0)
    np.testing.assert_allclose(tblite.dipole_au, native.dipole_au,
                               atol=2e-4, rtol=0)
    np.testing.assert_allclose(tblite.wiberg_bond_orders, native.wiberg_bond_orders,
                               atol=5e-5, rtol=0)
    assert tblite.metadata["engine"]["backend"] == "tblite"
    assert tblite.metadata["coefficient_axes"] == ["AO", "MO"]
    assert len(tblite.raw_json["orbital_map"]) == tblite.number_of_molecular_orbitals
    assert max(tblite.raw_json["shell_map"]) < mol.GetNumAtoms()
    assert len(tblite.raw_json["angular_momenta"]) == len(tblite.raw_json["shell_map"])


def test_tblite_is_explicit_and_default_stays_native(xtb_library):
    mol = molecule("O")
    explicit = run_xtb(mol, backend="tblite")
    default = run_xtb(mol, xtb_library=xtb_library)
    assert explicit.backend == "tblite"
    assert default.backend == "native"


def test_tblite_restrictions_are_clear():
    mol = molecule("O")
    with pytest.raises(ValueError, match="does not provide GFN0"):
        run_xtb(mol, method="gfn0", backend="tblite")
    with pytest.raises(XTBRunError, match="spin-polarized"):
        run_xtb(molecule("[CH3]"), backend="tblite")
    assert XTBCacheKey(mol, backend="tblite")


@pytest.mark.parametrize('method', ['GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['O', 'CCO', '[NH4+]', 'CS(C)=O', 'CCl'])
def test_population_descriptors_recover_engine_charges(method, smiles):
    from eled import CalcOrbitalPopulations, CalcFrontierDescriptors
    mol = molecule(smiles)
    result = run_xtb(mol, backend='tblite', method=method, keep_overlap=True)
    assert result.metadata['overlap_source'] == 'same tblite calculation'
    assert len(result.ao_atom_indices) == result.number_of_molecular_orbitals
    lowdin = CalcOrbitalPopulations(result)
    mulliken = CalcOrbitalPopulations(result, scheme='mulliken')
    # Independent engine output checks AO-to-atom mapping and density weighting,
    # including polarization d shells, hydrogens and a charged molecule.
    valence = [{1: 1, 6: 4, 7: 5, 8: 6, 16: 6, 17: 7}[a.GetAtomicNum()] for a in mol.GetAtoms()]
    np.testing.assert_allclose(np.array(valence) - mulliken['electron_populations_e'],
                               result.partial_charges, atol=1e-8, rtol=0)
    for population in (lowdin, mulliken):
        assert sum(population['electron_populations_e']) == pytest.approx(result.number_of_electrons, abs=1e-8)
        np.testing.assert_allclose(np.sum(population['atom_populations'], axis=0), 1., atol=1e-10)
    frontier = CalcFrontierDescriptors(mol, result)
    for label in ('homo', 'lumo'):
        assert 1 - 1e-10 <= frontier[label]['participation_ratio'] <= mol.GetNumAtoms() + 1e-10
    assert 0 <= frontier['homo_lumo_atom_overlap'] <= 1 + 1e-10
    assert frontier['homo_lumo_atom_centroid_separation_angstrom'] >= 0


def test_overlap_is_opt_in_and_cache_distinguishes_it():
    from eled import XTBCache, RunXTBBatch, CalcOrbitalPopulations, DescriptorUnavailableError
    mol = molecule('O')
    cache = XTBCache()
    ordinary = RunXTBBatch([mol], backend='tblite', cache=cache)[0]
    enriched = RunXTBBatch([mol], backend='tblite', keep_overlap=True, cache=cache)[0]
    assert len(cache) == 2
    assert not ordinary.overlap_matrix
    assert enriched.overlap_matrix
    assert ordinary.total_energy_eh == pytest.approx(enriched.total_energy_eh, abs=1e-12)
    with pytest.raises(DescriptorUnavailableError):
        CalcOrbitalPopulations(ordinary)


def test_real_benzene_degenerate_frontier_is_rotation_invariant():
    from dataclasses import replace
    from eled import CalcFrontierDescriptors
    mol = Chem.AddHs(Chem.MolFromSmiles('c1ccccc1'))
    conf = Chem.Conformer(12)
    for i in range(6):
        angle = i * np.pi / 3
        for index, radius in ((i, 1.4), (i+6, 2.48)):
            conf.SetAtomPosition(index, (radius*np.cos(angle), radius*np.sin(angle), 0.))
    conf.Set3D(True)
    mol.AddConformer(conf)
    result = run_xtb(mol, backend='tblite', keep_overlap=True)
    before = CalcFrontierDescriptors(mol, result)
    c = np.array(result.orbital_coefficients)
    for label in ('homo', 'lumo'):
        indices = before[label]['orbital_indices']
        assert len(indices) == 2
        c[:, indices] = c[:, indices] @ np.array([[.6, -.8], [.8, .6]])
    after = CalcFrontierDescriptors(mol, replace(result, orbital_coefficients=tuple(map(tuple, c))))
    for label in ('homo', 'lumo'):
        np.testing.assert_allclose(before[label]['atom_populations'], after[label]['atom_populations'], atol=1e-12)
