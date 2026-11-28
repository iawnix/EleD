"""Tests for Grimme FOD static-correlation and condensed linear response chi_AB."""
from __future__ import annotations

import json
import math
from pathlib import Path
import re
import subprocess
import tempfile

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (
    CalcFODFromResult,
    CalcLinearResponseMatrix,
    RunFOD,
    XTBCalculator,
    XTBOptions,
    XTBResults,
)
from eled.fod import BOLTZMANN_EV_PER_K
from eled.input import _input_metadata, _xyz_from_mol
from eled.options import _validate_calculation_settings
from eled.runtime import _find_xtb_library


def _make_mol(smiles: str, seed: int = 11):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _make_h2(dist_a: float):
    rw = Chem.RWMol()
    for _ in range(2):
        a = Chem.Atom(1)
        a.SetNoImplicit(True)
        rw.AddAtom(a)
    rw.AddBond(0, 1, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    mol.UpdatePropertyCache()
    conf = Chem.Conformer(2)
    conf.Set3D(True)
    conf.SetAtomPosition(0, (0.0, 0.0, 0.0))
    conf.SetAtomPosition(1, (dist_a, 0.0, 0.0))
    mol.AddConformer(conf, assignId=True)
    return mol


def _rotate_translate_permute(mol, perm):
    theta = 0.73
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([1.2, -0.7, 0.5])
    reordered = Chem.RenumberAtoms(mol, list(perm))
    conf = reordered.GetConformer()
    orig = mol.GetConformer()
    for new_i, old_i in enumerate(perm):
        p = orig.GetAtomPosition(old_i)
        conf.SetAtomPosition(new_i, (rot @ np.array([p.x, p.y, p.z]) + shift).tolist())
    return reordered


def _synthetic_h2_result(
    *,
    e_occ: float = -10.0,
    e_vir: float = 2.0,
    occ=(2.0, 0.0),
    uhf: int | None = 0,
    etemp: float = 5000.0,
    with_overlap: bool = True,
    s12: float = 0.4,
):
    mol = _make_h2(0.74)
    meta = _input_metadata(mol, -1, 0, uhf or 0, _xyz_from_mol(mol, -1))
    meta['calculation_settings'] = _validate_calculation_settings(electronicTemperature=etemp)
    meta['uhf'] = uhf
    s = np.array([[1.0, s12], [s12, 1.0]])
    c0 = np.array([1.0, 1.0]) / math.sqrt(2.0 * (1.0 + s12))
    c1 = np.array([1.0, -1.0]) / math.sqrt(2.0 * (1.0 - s12))
    c = np.column_stack((c0, c1))
    extra = {}
    if with_overlap:
        extra = dict(
            overlapMatrix=tuple(map(tuple, s.tolist())),
            orbitalCoefficients=tuple(map(tuple, c.tolist())),
            aoAtomMap=(0, 1),
            aoShellMap=(0, 1),
            shellAtomMap=(0, 1),
            shellAngularMomenta=(0, 0),
        )
    return XTBResults(
        method='GFN2-xTB',
        totalEnergyHartree=-1.0,
        numOrbitals=2,
        numElectrons=int(round(sum(occ))),
        numUnpairedElectrons=uhf,
        orbitalEnergiesEV=(e_occ, e_vir),
        orbitalOccupancies=tuple(occ),
        atomicCharges=(0.0, 0.0),
        metadata=meta,
        **extra,
    )


# ---------------------------------------------------------------------------
# Pure-Python unit tests: FOD
# ---------------------------------------------------------------------------


def test_fod_analytical_two_level_and_temperature_limits():
    # Symmetric 2-level closed-shell: E_F = midpoint, f_occ = 1/(1 + exp(-dE/(2 kT))),
    # f_vir = 1/(1 + exp(+dE/(2 kT))). Per spin: (1 - f_occ) + f_vir = 2 * f_vir.
    # Total over 2 spins: N_FOD = 4 / (1 + exp(dE / (2 kT))).
    e_occ, e_vir = -2.0, 0.0
    temp = 5000.0
    res = _synthetic_h2_result(e_occ=e_occ, e_vir=e_vir, etemp=temp)
    out = CalcFODFromResult(res)
    kt = BOLTZMANN_EV_PER_K * temp
    expected_nfod = 4.0 / (1.0 + math.exp((e_vir - e_occ) / (2.0 * kt)))
    assert out['n_fod_e'] == pytest.approx(expected_nfod, rel=1e-10)
    assert out['fod_fermi_alpha_ev'] == pytest.approx(-1.0, abs=1e-10)
    assert out['fod_fermi_beta_ev'] == pytest.approx(-1.0, abs=1e-10)
    assert out['self_consistent_temperature'] is True
    assert out['atomic_fod_available'] is True
    assert out['shell_fod_available'] is True
    assert sum(a['fod_pop_e'] for a in out['atoms']) == pytest.approx(out['n_fod_e'], abs=1e-12)
    assert out['atoms'][0]['fod_s_e'] == pytest.approx(out['atoms'][0]['fod_pop_e'], abs=1e-12)
    assert out['atoms'][0]['fod_p_e'] == 0.0
    assert out['atoms'][0]['fod_d_e'] == 0.0

    # Zero electronic temperature -> zero FOD.
    out_zero = CalcFODFromResult(res, electronicTemperature=0.0)
    assert out_zero['n_fod_e'] == 0.0
    assert out_zero['self_consistent_temperature'] is False

    # Higher temperature increases N_FOD monotonically.
    out_hot = CalcFODFromResult(res, electronicTemperature=10000.0)
    assert out_hot['n_fod_e'] > out['n_fod_e']
    json.dumps(out, allow_nan=False)


def test_fod_without_overlap_marks_atoms_missing_and_validates_inputs():
    res = _synthetic_h2_result(with_overlap=False)
    out = CalcFODFromResult(res)
    assert out['atomic_fod_available'] is False
    assert out['atoms'] == []
    assert out['missing']['atoms'] == 'overlap_basis_not_retained'

    with pytest.raises(ValueError, match='scheme'):
        CalcFODFromResult(res, scheme='hirshfeld')
    with pytest.raises(ValueError, match='electronicTemperature'):
        CalcFODFromResult(res, electronicTemperature=-10.0)
    with pytest.raises(ValueError, match='interface'):
        RunFOD(_make_h2(0.74), interface='cli')


# ---------------------------------------------------------------------------
# Pure-Python unit tests: Condensed Linear Response chi_AB
# ---------------------------------------------------------------------------


def test_linear_response_analytical_h2_and_sum_rule():
    # For symmetric H2 (1 occ, 1 vir) with overlap s12:
    # Lowdin:   Q^0_{0,1} = +1/2                  -> chi_00 = -1 / dE
    # Mulliken: Q^0_{0,1} = +1 / (2*sqrt(1-s12^2)) -> chi_00 = -1 / (dE * (1 - s12^2))
    e_occ, e_vir, s12 = -10.0, 2.0, 0.4
    de = e_vir - e_occ
    res = _synthetic_h2_result(e_occ=e_occ, e_vir=e_vir, s12=s12)
    expected_diag = {
        'lowdin': -1.0 / de,
        'mulliken': -1.0 / (de * (1.0 - s12 * s12)),
    }
    for scheme, diag in expected_diag.items():
        out = CalcLinearResponseMatrix(res, scheme=scheme)
        assert out['available'] is True
        chi = np.asarray(out['chi_e_per_ev'])
        expected = np.array([[diag, -diag], [-diag, diag]])
        np.testing.assert_allclose(chi, expected, atol=1e-12)
        assert out['sum_rule_max_error_e_per_ev'] < 1e-14
        assert out['atoms'][0]['chi_diag_e_per_ev'] == pytest.approx(diag)
        assert out['atoms'][0]['neg_chi_diag_e_per_ev'] == pytest.approx(-diag)
        assert out['atoms'][0]['nonlocal_coupling_e_per_ev'] == pytest.approx(-diag)
        json.dumps(out, allow_nan=False)


@pytest.mark.parametrize(
    'kwargs,expected_reason',
    [
        ({'uhf': 1, 'occ': (1.0, 0.0)}, 'open_shell_not_supported'),
        ({'uhf': None, 'occ': (2.0, 0.0)}, 'unknown_spin_state'),
        ({'uhf': 0, 'occ': (1.5, 0.5)}, 'fractional_occupation_not_supported'),
        ({'uhf': 0, 'e_occ': -1.0, 'e_vir': -1.0}, 'degenerate_homo_lumo_gap'),
    ],
)
def test_linear_response_withheld_explicitly_without_zeros(kwargs, expected_reason):
    res = _synthetic_h2_result(**kwargs)
    out = CalcLinearResponseMatrix(res)
    assert out['available'] is False
    assert out['chi_e_per_ev'] is None
    assert out['chi_e_per_ha'] is None
    assert out['atoms'] == []
    assert out['missing']['chi_e_per_ev'] == expected_reason


# ---------------------------------------------------------------------------
# Integration tests (require installed libxtb + rdXTBTools)
# ---------------------------------------------------------------------------


@pytest.fixture(scope='module')
def library():
    return _find_xtb_library(None)


@pytest.mark.integration
def test_fod_stretched_h2_static_correlation_and_native_compiled_parity(library):
    h2_eq = _make_h2(0.74)
    h2_stretched = _make_h2(2.60)
    opts = XTBOptions(method='GFN2-xTB', xtbLibrary=library)

    fod_eq = RunFOD(h2_eq, electronicTemperature=5000.0, options=opts, interface='native')
    fod_str_native = RunFOD(h2_stretched, electronicTemperature=5000.0, options=opts, interface='native')
    fod_str_compiled = RunFOD(h2_stretched, electronicTemperature=5000.0, options=opts, interface='compiled')

    assert fod_str_native['n_fod_e'] > fod_eq['n_fod_e'] + 0.1
    assert fod_str_native['n_fod_e'] == pytest.approx(fod_str_compiled['n_fod_e'], abs=1e-10)
    np.testing.assert_allclose(
        [a['fod_pop_e'] for a in fod_str_native['atoms']],
        [a['fod_pop_e'] for a in fod_str_compiled['atoms']],
        atol=1e-10,
    )
    assert abs(fod_str_native['fod_sum_error_e']) < 1e-10
    assert fod_str_native['atoms'][0]['fod_pop_e'] == pytest.approx(
        fod_str_native['atoms'][1]['fod_pop_e'], abs=1e-10
    )


@pytest.mark.integration
def test_fod_matches_xtb_cli_fod_on_stretched_h2_and_water(library):
    xtb_bin = Path(library).resolve().parent.parent / 'bin/xtb'
    if not xtb_bin.exists():
        pytest.skip('xtb CLI binary not present alongside libxtb')

    mol = _make_h2(2.40)
    conf = mol.GetConformer()
    xyz_lines = ['2', 'stretched H2']
    for i in range(2):
        p = conf.GetAtomPosition(i)
        xyz_lines.append(f'H {p.x:.10f} {p.y:.10f} {p.z:.10f}')

    with tempfile.TemporaryDirectory() as tmp:
        xyz_path = Path(tmp) / 'h2.xyz'
        xyz_path.write_text('\n'.join(xyz_lines) + '\n')
        proc = subprocess.run(
            [str(xtb_bin), str(xyz_path), '--gfn', '2', '--fod', '--etemp', '5000'],
            cwd=tmp,
            capture_output=True,
            text=True,
            check=True,
        )
    match = re.search(r'NFOD\s*:\s*([0-9.]+)', proc.stdout)
    assert match is not None, f'Could not parse NFOD from xtb output:\n{proc.stdout}'
    xtb_nfod = float(match.group(1))
    pop_matches = re.findall(r'^\s*\d+H\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)', proc.stdout, re.M)
    assert len(pop_matches) == 2

    eled_fod = RunFOD(mol, electronicTemperature=5000.0, options=XTBOptions(xtbLibrary=library))
    assert eled_fod['n_fod_e'] == pytest.approx(xtb_nfod, abs=2e-4)
    for idx, (tot_s, s_s, _, _) in enumerate(pop_matches):
        assert eled_fod['atoms'][idx]['fod_pop_e'] == pytest.approx(float(tot_s), abs=2e-4)
        assert eled_fod['atoms'][idx]['fod_s_e'] == pytest.approx(float(s_s), abs=1e-3)


@pytest.mark.integration
def test_linear_response_integration_invariance_and_parity(library):
    from eled import rdXTBTools

    mol = _make_mol('O', seed=17)
    calc = XTBCalculator(method='GFN2-xTB', keepOverlap=True, electronicTemperature=0.0, xtbLibrary=library)
    res_native = calc.Calculate(mol)
    ok, res_compiled = rdXTBTools.RunMol(
        mol,
        method='GFN2-xTB',
        keepOverlap=True,
        keepOrbitalCoefficients=True,
        electronicTemperature=0.0,
        xtbLibrary=library,
    )
    assert ok

    lr_native = CalcLinearResponseMatrix(res_native, scheme='lowdin')
    lr_compiled = CalcLinearResponseMatrix(res_compiled, scheme='lowdin')
    assert lr_native['available'] is True
    chi_n = np.asarray(lr_native['chi_e_per_ev'])
    chi_c = np.asarray(lr_compiled['chi_e_per_ev'])
    np.testing.assert_allclose(chi_n, chi_c, atol=1e-10)
    np.testing.assert_allclose(chi_n, chi_n.T, atol=1e-12)
    assert lr_native['sum_rule_max_error_e_per_ev'] < 1e-12
    assert np.all(np.diag(chi_n) < 0.0)
    eigvals = np.linalg.eigvalsh(chi_n)
    assert np.max(eigvals) < 1e-10  # negative semi-definite

    # Rigid rotation + atom permutation invariance.
    perm = [2, 0, 1]
    mol_rot = _rotate_translate_permute(mol, perm)
    lr_rot = CalcLinearResponseMatrix(calc.Calculate(mol_rot), scheme='lowdin')
    chi_rot = np.asarray(lr_rot['chi_e_per_ev'])
    np.testing.assert_allclose(chi_rot, chi_n[np.ix_(perm, perm)], atol=1e-8)

    # Open-shell radical is explicitly withheld.
    oh = _make_mol('[OH]', seed=5)
    res_oh = calc.Calculate(oh)
    lr_oh = CalcLinearResponseMatrix(res_oh)
    assert lr_oh['available'] is False
    assert lr_oh['missing']['chi_e_per_ev'] == 'open_shell_not_supported'
