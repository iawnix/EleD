"""Direct compiled-module comparisons, including orientation-sensitive d shells."""

import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem

from rdkit_eth import XTBNativeNotFoundError
from rdkit_eth.backends.tblite import run_xtb_tblite
from rdkit_eth.orbital_alignment import compare_orbitals, libxtb_to_tblite_rows
from rdkit_eth.runtime import _find_xtb_library

if (importlib.util.find_spec("tblite") is None
        or importlib.util.find_spec("rdkit_eth.rdXTBTools") is None):
    pytest.skip("compiled rdXTBTools and optional tblite are required", allow_module_level=True)

from rdkit_eth import rdXTBTools

pytestmark = pytest.mark.integration
CASES = {c["name"]: c for c in json.loads(
    (Path(__file__).parent / "data/orbital_alignment_cases.json").read_text())["cases"]}


@pytest.fixture(scope="module")
def library():
    try:
        return _find_xtb_library(None)
    except XTBNativeNotFoundError:
        pytest.skip("libxtb is required")


@pytest.mark.parametrize("method", ["GFN1-xTB", "GFN2-xTB"])
@pytest.mark.parametrize("name", ["water", "phosphine", "dimethyl_sulfoxide", "chloromethane", "benzene", "acetate"])
def test_compiled_orbitals_in_tblite_metric(library, method, name):
    case = CASES[name]
    mol = Chem.AddHs(Chem.MolFromSmiles(case["smiles"]))
    # Apply a nontrivial rigid rotation to exercise angular conventions beyond
    # the fixture orientation. No new conformer generation or optimization.
    rotation, _ = np.linalg.qr(np.array([[1., 2., 3.], [4., -2., 1.], [.5, .3, -.8]]))
    if np.linalg.det(rotation) < 0:
        rotation[:, 0] *= -1
    positions = np.asarray(case["positions_angstrom"]) @ rotation
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, xyz in enumerate(positions):
        conf.SetAtomPosition(i, xyz.tolist())
    mol.AddConformer(conf)
    before = mol.ToBinary()
    ok, native = rdXTBTools.RunMol(mol, method=method, xtbLibrary=library,
                                  keepOrbitalCoefficients=True, accuracy=.01)
    assert ok, native.errorMessage
    other = run_xtb_tblite(mol, method=method, save_integrals=True, accuracy=.01)
    maps = other.raw_json
    rows, phases, labels = libxtb_to_tblite_rows(
        case["atomic_numbers"], native.parameterFile, maps["orbital_map"],
        maps["shell_map"], maps["angular_momenta"], native.numOrbitals)
    c = native.GetOrbitalCoefficients()[rows, :] * phases[:, None]
    result = compare_orbitals(other.orbital_coefficients, c, maps["overlap_matrix"],
                              other.orbital_energies_ev, native.GetOrbitalEnergies(),
                              other.orbital_occupations, native.GetOrbitalOccupancies())
    assert abs(native.totalEnergyHartree - other.total_energy_eh) < 2e-7
    assert result["reference_orthonormality_error"] < 1e-8
    assert result["candidate_orthonormality_error"] < 1e-8
    assert result["max_projector_distance"] < 1e-3
    assert result["density_orthonormal_max_abs"] < 1e-4
    assert result["density_trace_difference"] < 1e-6
    assert [label["atom"] for label in labels] == np.asarray(maps["shell_map"])[maps["orbital_map"]].tolist()
    if name in {"phosphine", "dimethyl_sulfoxide", "chloromethane"}:
        assert -1 in phases  # a d-shell phase, distinct from arbitrary MO signs
        wrong = native.GetOrbitalCoefficients()[rows, :]
        s = np.asarray(maps["overlap_matrix"])
        assert np.max(abs(wrong.T @ s @ wrong - np.eye(native.numOrbitals))) > 1e-3
    assert mol.ToBinary() == before
