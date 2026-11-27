"""Tests for Tier 2 extended libxtb properties ABI, C++/Python parity, and descriptors."""
from __future__ import annotations

import ctypes
from dataclasses import replace
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from eled import (
    CalcESP,
    DescriptorUnavailableError,
    XTBCalculator,
    XTBResults,
)
from eled.properties import (
    CalcEnergyDecomposition,
    CalcFiniteFieldPolarizability,
    CalcMultipoleDescriptors,
    CalcNuclearESP,
    CalcPolarizabilityDescriptors,
)
from eled.runtime import _PROPERTIES_EXPORT_PATCH, _find_xtb_library, _native_properties_build


def _make_mol(smiles: str = 'O', seed: int = 17):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol, randomSeed=seed) == 0
    return mol


def _rotate_translate_permute(mol, perm):
    theta = 0.73
    c, s = np.cos(theta), np.sin(theta)
    rot = np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])
    shift = np.array([1.7, -2.3, 0.9])
    reordered = Chem.RenumberAtoms(mol, list(perm))
    conf = reordered.GetConformer()
    orig_conf = mol.GetConformer()
    for new_idx, old_idx in enumerate(perm):
        pt = orig_conf.GetAtomPosition(old_idx)
        pos = np.array([pt.x, pt.y, pt.z], dtype=float)
        moved = rot @ pos + shift
        conf.SetAtomPosition(new_idx, moved.tolist())
    return reordered


# ---------------------------------------------------------------------------
# Pure-Python contract & analytical tests (run in both .runtime-env and .eled)
# ---------------------------------------------------------------------------

def test_properties_patch_lock_consistency():
    root = Path(__file__).resolve().parents[1]
    lock = json.loads((root / 'tools/xtb-source.json').read_text())['properties_extension']
    patch_bytes = (root / 'tools' / lock['patch']).read_bytes()
    patch_hash = hashlib.sha256(patch_bytes).hexdigest()
    assert patch_hash == lock['patch_sha256'] == _PROPERTIES_EXPORT_PATCH
    assert _PROPERTIES_EXPORT_PATCH in (root / 'cpp/rdXTBTools/XTBTools.cpp').read_text()
    assert b'/home/' not in patch_bytes
    assert b'diff --git a/include/xtb.h b/include/xtb.h' in patch_bytes


def test_synthetic_energy_decomp_and_missing_semantics():
    gfn0 = XTBResults(
        method='GFN0-xTB',
        totalEnergyHartree=-5.0,
        energyDecompositionFlags=1 | 2 | 4 | 8 | 64 | 128,
        electronicEnergyHartree=-4.0,
        repulsionEnergyHartree=0.5,
        electrostaticEnergyHartree=-1.0,
        dispersionEnergyHartree=-0.2,
        halogenBondEnergyHartree=-0.3,
        metadata={'calculation_status': 'success'},
    )
    out = CalcEnergyDecomposition(gfn0)
    assert out['available'] is True
    assert out['e_aes_ha'] is None and out['missing']['e_aes_ha'] == 'not_applicable_to_gfn0'
    assert out['e_axc_ha'] is None and out['missing']['e_axc_ha'] == 'not_applicable_to_gfn0'
    assert out['g_solv_ha'] is None and out['missing']['g_solv_ha'] == 'gas_phase_no_solvation'
    assert out['fermi_alpha_ha'] is None and out['missing']['fermi_alpha_ha'] == 'fermi_smearing_disabled_low_etemp'
    assert out['xb_term_definition'] == 'short_range_bond_correction'
    json.dumps(out, allow_nan=False)

    unavailable = XTBResults(method='GFN2-xTB', totalEnergyHartree=-5.0, metadata={'calculation_status': 'success'})
    out_u = CalcEnergyDecomposition(unavailable)
    assert out_u['available'] is False
    assert out_u['e_elec_ha'] is None
    assert out_u['missing']['e_elec_ha'] == 'properties_extension_unavailable'


def test_synthetic_multipole_rotational_invariants_and_gfn0_missing():
    # Diagonal traceless quadrupole Q = diag(1, -2, 1) -> frob = sqrt(6), anisotropy = sqrt(9) = 3.
    q_diag = (1.0, 0.0, -2.0, 0.0, 0.0, 1.0)
    d_vec = (0.3, -0.4, 0.0)  # norm 0.5
    res = XTBResults(
        method='GFN2-xTB',
        totalEnergyHartree=-1.0,
        atomicDipolesAU=(d_vec,),
        atomicQuadrupolesAU=(q_diag,),
        molecularQuadrupoleAU=q_diag,
        metadata={'calculation_status': 'success', 'atom_count': 1},
    )
    out = CalcMultipoleDescriptors(res)
    assert out['available'] is True
    assert out['mol_quadrupole_frobenius_ebohr2'] == pytest.approx(np.sqrt(6.0), abs=1e-12)
    assert out['mol_quadrupole_anisotropy_ebohr2'] == pytest.approx(3.0, abs=1e-12)
    np.testing.assert_allclose(out['mol_quadrupole_eigenvalues_ebohr2'], [-2.0, 1.0, 1.0], atol=1e-12)
    assert out['atoms'][0]['dip_norm_ebohr'] == pytest.approx(0.5, abs=1e-12)

    # Rotate tensor Q -> R Q R^T and vector d -> R d; invariants must match identically.
    theta = 0.9
    r = np.array([[np.cos(theta), -np.sin(theta), 0.0], [np.sin(theta), np.cos(theta), 0.0], [0.0, 0.0, 1.0]])
    q_rot_mat = r @ np.diag([1.0, -2.0, 1.0]) @ r.T
    q_rot = (q_rot_mat[0, 0], q_rot_mat[0, 1], q_rot_mat[1, 1], q_rot_mat[0, 2], q_rot_mat[1, 2], q_rot_mat[2, 2])
    d_rot = tuple((r @ np.array(d_vec)).tolist())
    res_rot = replace(res, atomicDipolesAU=(d_rot,), atomicQuadrupolesAU=(q_rot,), molecularQuadrupoleAU=q_rot)
    out_rot = CalcMultipoleDescriptors(res_rot)
    assert out_rot['mol_quadrupole_frobenius_ebohr2'] == pytest.approx(out['mol_quadrupole_frobenius_ebohr2'], abs=1e-12)
    assert out_rot['mol_quadrupole_anisotropy_ebohr2'] == pytest.approx(out['mol_quadrupole_anisotropy_ebohr2'], abs=1e-12)
    np.testing.assert_allclose(out_rot['mol_quadrupole_eigenvalues_ebohr2'], out['mol_quadrupole_eigenvalues_ebohr2'], atol=1e-12)
    assert out_rot['atoms'][0]['dip_norm_ebohr'] == pytest.approx(out['atoms'][0]['dip_norm_ebohr'], abs=1e-12)

    # GFN0 missing semantics and schema/units parity with available=True
    gfn0 = XTBResults(method='GFN0-xTB', totalEnergyHartree=-1.0, metadata={'calculation_status': 'success'})
    out0 = CalcMultipoleDescriptors(gfn0)
    assert out0['available'] is False
    assert out0['mol_quadrupole_ebohr2'] is None
    assert out0['missing']['mol_quadrupole_ebohr2'] == 'not_computed_by_gfn0'
    assert out0['missing']['multipole_model'] == 'not_computed_by_gfn0'
    assert out0['units'] == out['units']
    json.dumps(out0, allow_nan=False)

    # Nuclear ESP available=False shares full units schema with available=True
    epn_ok = CalcNuclearESP(XTBResults(method='GFN2-xTB', totalEnergyHartree=-1.0, nuclearESPAU=(-1.0,), metadata={'calculation_status': 'success'}))
    epn_miss = CalcNuclearESP(gfn0)
    assert epn_miss['available'] is False
    assert epn_miss['units'] == epn_ok['units']
    json.dumps(epn_miss, allow_nan=False)


def test_synthetic_polarizability_gfn1_missing_and_finite_field_error():
    gfn1 = XTBResults(method='GFN1-xTB', totalEnergyHartree=-1.0, metadata={'calculation_status': 'success'})
    out = CalcPolarizabilityDescriptors(gfn1)
    assert out['available'] is False
    assert out['mol_pol_bohr3'] is None
    assert out['missing']['mol_pol_bohr3'] == 'not_computed_by_gfn1_d3_dispersion'
    assert out['finite_field_supported'] is False
    json.dumps(out, allow_nan=False)

    with pytest.raises(DescriptorUnavailableError, match='Finite-field polarizability'):
        CalcFiniteFieldPolarizability(_make_mol('O'))


def test_synthetic_nuclear_esp_and_result_roundtrip():
    res = XTBResults(
        method='GFN2-xTB',
        totalEnergyHartree=-5.0,
        energyDecompositionFlags=79,
        electronicEnergyHartree=-5.5,
        repulsionEnergyHartree=0.5,
        molecularPolarizabilityAU=9.2,
        molecularC6AU=45.0,
        molecularC8AU=300.0,
        atomicPolarizabilitiesAU=(5.2, 2.0, 2.0),
        atomicC6AU=(25.0, 3.0, 3.0),
        atomicDipolesAU=((0.0, 0.0, 0.1), (0.0, 0.05, -0.05), (0.0, -0.05, -0.05)),
        atomicQuadrupolesAU=((0.1, 0.0, -0.2, 0.0, 0.0, 0.1),) * 3,
        molecularQuadrupoleAU=(0.3, 0.0, -0.6, 0.0, 0.0, 0.3),
        nuclearESPAU=(-22.3, -1.1, -1.1),
        metadata={'calculation_status': 'success'},
    )
    restored = XTBResults.FromDict(json.loads(json.dumps(res.ToDict(), allow_nan=False)))
    assert restored == res
    epn = CalcNuclearESP(restored)
    assert epn['available'] is True
    assert len(epn['atoms']) == 3
    assert epn['epn_min_au'] == pytest.approx(-22.3)
    assert epn['epn_max_au'] == pytest.approx(-1.1)
    json.dumps(epn, allow_nan=False)

    missing_epn = CalcNuclearESP(XTBResults(method='GFN2-xTB', totalEnergyHartree=-1.0, metadata={'calculation_status': 'success'}))
    assert missing_epn['available'] is False
    assert missing_epn['epn_min_au'] is None
    assert missing_epn['missing']['atoms'] == 'properties_extension_unavailable'


# ---------------------------------------------------------------------------
# Integration tests against rebuilt libxtb + compiled rdXTBTools
# ---------------------------------------------------------------------------

@pytest.fixture(scope='module')
def library():
    path = _find_xtb_library(None)
    if not _native_properties_build(path):
        pytest.skip('verified extended properties libxtb build required')
    return path


@pytest.mark.integration
@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
@pytest.mark.parametrize('smiles', ['O', 'CCO', '[NH4+]', '[OH]'])
def test_native_compiled_properties_parity(library, method, smiles):
    from eled import rdXTBTools
    mol = _make_mol(smiles)
    native = XTBCalculator(method=method, xtbLibrary=library).Calculate(mol)
    ok, compiled = rdXTBTools.RunMol(mol, method=method, xtbLibrary=library)
    assert ok, compiled.errorMessage
    assert compiled.propertiesExportVerified

    # Energy decomposition parity & sum identities
    ed_nat = CalcEnergyDecomposition(native)
    ed_cmp = CalcEnergyDecomposition(compiled)
    assert ed_nat['available'] is True and ed_cmp['available'] is True
    assert ed_nat['missing'] == ed_cmp['missing']
    for key in (
        'e_total_ha', 'e_elec_ha', 'e_rep_ha', 'e_es_ha', 'e_aes_ha', 'e_axc_ha',
        'e_disp_ha', 'e_xb_ha', 'g_solv_ha', 'fermi_alpha_ha', 'fermi_beta_ha',
        'fermi_alpha_ev', 'fermi_beta_ev',
    ):
        if ed_nat[key] is None:
            assert ed_cmp[key] is None
        else:
            assert ed_cmp[key] == pytest.approx(ed_nat[key], abs=1e-9, rel=1e-9)

    if method == 'GFN2-xTB':
        assert ed_nat['e_elec_ha'] + ed_nat['e_rep_ha'] == pytest.approx(ed_nat['e_total_ha'], abs=1e-9)
        assert ed_nat['e_xb_ha'] is None and ed_nat['missing']['e_xb_ha'] == 'not_applicable_to_gfn2'
    elif method == 'GFN1-xTB':
        assert (
            ed_nat['e_elec_ha'] + ed_nat['e_rep_ha'] + ed_nat['e_disp_ha'] + ed_nat['e_xb_ha']
            == pytest.approx(ed_nat['e_total_ha'], abs=1e-9)
        )
        assert ed_nat['e_aes_ha'] is None and ed_nat['e_axc_ha'] is None
    else:
        assert (
            ed_nat['e_elec_ha'] + ed_nat['e_rep_ha'] + ed_nat['e_es_ha'] + ed_nat['e_disp_ha'] + ed_nat['e_xb_ha']
            == pytest.approx(ed_nat['e_total_ha'], abs=1e-9)
        )

    # Polarizability parity
    pol_nat = CalcPolarizabilityDescriptors(native, mol=mol)
    pol_cmp = CalcPolarizabilityDescriptors(compiled, mol=mol)
    assert pol_nat['available'] == pol_cmp['available'] == (method in ('GFN0-xTB', 'GFN2-xTB'))
    if pol_nat['available']:
        assert pol_cmp['mol_pol_bohr3'] == pytest.approx(pol_nat['mol_pol_bohr3'], abs=1e-9)
        assert pol_cmp['mol_c6_ha_bohr6'] == pytest.approx(pol_nat['mol_c6_ha_bohr6'], abs=1e-8)
        assert pol_cmp['mol_c8_ha_bohr8'] == pytest.approx(pol_nat['mol_c8_ha_bohr8'], abs=1e-7)
        np.testing.assert_allclose(compiled.GetAtomicPolarizabilities(), native.GetAtomicPolarizabilities(), atol=1e-9)
        np.testing.assert_allclose(compiled.GetAtomicC6(), native.GetAtomicC6(), atol=1e-8)
        assert sum(r['alpha_bohr3'] for r in pol_nat['atoms']) == pytest.approx(pol_nat['mol_pol_bohr3'], abs=1e-9)
    else:
        assert pol_nat['missing']['mol_pol_bohr3'] == 'not_computed_by_gfn1_d3_dispersion'
        with pytest.raises(ValueError, match='unavailable'):
            compiled.GetAtomicPolarizabilities()
        with pytest.raises(ValueError, match='unavailable'):
            native.GetAtomicPolarizabilities()

    # Multipoles parity
    mp_nat = CalcMultipoleDescriptors(native, mol=mol)
    mp_cmp = CalcMultipoleDescriptors(compiled, mol=mol)
    assert mp_nat['available'] == mp_cmp['available'] == (method in ('GFN1-xTB', 'GFN2-xTB'))
    if mp_nat['available']:
        np.testing.assert_allclose(compiled.GetAtomicDipoles(), native.GetAtomicDipoles(), atol=1e-8)
        np.testing.assert_allclose(compiled.GetAtomicQuadrupoles(), native.GetAtomicQuadrupoles(), atol=1e-8)
        np.testing.assert_allclose(compiled.GetMolecularQuadrupole(), native.GetMolecularQuadrupole(), atol=1e-8)
        assert abs(mp_nat['mol_quadrupole_trace_ebohr2']) < 1e-11
        for row in mp_nat['atoms']:
            assert abs(row['quad_trace_ebohr2']) < 1e-11
    else:
        assert mp_nat['missing']['mol_quadrupole_ebohr2'] == 'not_computed_by_gfn0'
        with pytest.raises(ValueError, match='unavailable'):
            compiled.GetAtomicDipoles()

    # Nuclear ESP parity & owned copy check
    epn_nat = CalcNuclearESP(native, mol=mol)
    epn_cmp = CalcNuclearESP(compiled, mol=mol)
    assert epn_nat['available'] is True and epn_cmp['available'] is True
    np.testing.assert_allclose(compiled.GetNuclearESP(), native.GetNuclearESP(), atol=1e-9)
    owned = compiled.GetNuclearESP()
    assert owned.flags.owndata
    orig = owned.copy()
    owned[:] = 999.0
    np.testing.assert_array_equal(compiled.GetNuclearESP(), orig)

    # All outputs JSON-safe
    for obj in (ed_nat, ed_cmp, pol_nat, pol_cmp, mp_nat, mp_cmp, epn_nat, epn_cmp):
        json.dumps(obj, allow_nan=False)


@pytest.mark.integration
@pytest.mark.parametrize('method', ['GFN0-xTB', 'GFN1-xTB', 'GFN2-xTB'])
def test_properties_rigid_motion_and_permutation_invariance(library, method):
    mol = _make_mol('CCO', seed=5)
    perm = list(reversed(range(mol.GetNumAtoms())))
    moved = _rotate_translate_permute(mol, perm)

    r1 = XTBCalculator(method=method, xtbLibrary=library).Calculate(mol)
    r2 = XTBCalculator(method=method, xtbLibrary=library).Calculate(moved)

    ed1, ed2 = CalcEnergyDecomposition(r1), CalcEnergyDecomposition(r2)
    for key in ('e_total_ha', 'e_elec_ha', 'e_rep_ha', 'e_es_ha', 'e_disp_ha', 'fermi_alpha_ha'):
        assert ed1[key] == pytest.approx(ed2[key], abs=1e-7)

    epn1 = r1.GetNuclearESP()
    epn2 = r2.GetNuclearESP()
    np.testing.assert_allclose(epn1[perm], epn2, atol=1e-6)

    if r1.polarizabilityAvailable:
        assert r1.molecularPolarizabilityAU == pytest.approx(r2.molecularPolarizabilityAU, abs=1e-7)
        np.testing.assert_allclose(r1.GetAtomicPolarizabilities()[perm], r2.GetAtomicPolarizabilities(), atol=1e-7)
        np.testing.assert_allclose(r1.GetAtomicC6()[perm], r2.GetAtomicC6(), atol=1e-6)

    if r1.multipolesAvailable:
        mp1, mp2 = CalcMultipoleDescriptors(r1), CalcMultipoleDescriptors(r2)
        for new_idx, old_idx in enumerate(perm):
            assert mp1['atoms'][old_idx]['dip_norm_ebohr'] == pytest.approx(
                mp2['atoms'][new_idx]['dip_norm_ebohr'], abs=1e-6
            )
            assert mp1['atoms'][old_idx]['quad_frobenius_ebohr2'] == pytest.approx(
                mp2['atoms'][new_idx]['quad_frobenius_ebohr2'], abs=1e-6
            )
            assert mp1['atoms'][old_idx]['quad_anisotropy_ebohr2'] == pytest.approx(
                mp2['atoms'][new_idx]['quad_anisotropy_ebohr2'], abs=1e-6
            )
            np.testing.assert_allclose(
                mp1['atoms'][old_idx]['quad_eigenvalues_ebohr2'],
                mp2['atoms'][new_idx]['quad_eigenvalues_ebohr2'],
                atol=1e-6,
            )


@pytest.mark.integration
def test_zero_electronic_temperature_omits_fermi_and_nuclear_esp_matches_esp_limit(library):
    mol = _make_mol('O', seed=17)
    res_zero = XTBCalculator(method='GFN2-xTB', electronicTemperature=0.0, xtbLibrary=library).Calculate(mol)
    ed = CalcEnergyDecomposition(res_zero)
    assert ed['fermi_alpha_ha'] is None and ed['fermi_beta_ha'] is None
    assert ed['missing']['fermi_alpha_ha'] == 'fermi_smearing_disabled_low_etemp'

    # Check nuclear ESP matches spherical limit of CalcESP(r) - Z_O / |r - R_O|
    coords = np.asarray(mol.GetConformer().GetPositions(), dtype=float)
    ox_pos = coords[0]
    delta_bohr = 1e-4
    delta_ang = delta_bohr / 1.8897261254578281
    probes = [
        tuple((ox_pos + np.array([s * delta_ang, 0.0, 0.0])).tolist()) for s in (-1, 1)
    ] + [
        tuple((ox_pos + np.array([0.0, s * delta_ang, 0.0])).tolist()) for s in (-1, 1)
    ] + [
        tuple((ox_pos + np.array([0.0, 0.0, s * delta_ang])).tolist()) for s in (-1, 1)
    ]
    res_esp = XTBCalculator(method='GFN2-xTB', xtbLibrary=library).Calculate(mol, espPointsAngstrom=probes)
    esp_vals = np.asarray(CalcESP(res_esp), dtype=float)
    z_valence_oxygen = 6.0
    limit_est = float(np.mean(esp_vals - z_valence_oxygen / delta_bohr))
    epn_ox = float(res_esp.GetNuclearESP()[0])
    assert abs(limit_est - epn_ox) < 1e-3


@pytest.mark.integration
def test_c_api_properties_copy_and_invalidation(library):
    from eled.backends.native import _NativeXTB
    api = _NativeXTB(library)
    api.bind_properties_api()
    env = api.new_environment()
    api.set_verbosity(env, 0)
    n = ctypes.c_int(3)
    q = ctypes.c_double(0.0)
    uhf = ctypes.c_int(0)
    at = (ctypes.c_int * 3)(8, 1, 1)
    xyz = (ctypes.c_double * 9)(0.0, 0.0, 0.0, 0.0, 0.0, 1.8, 1.7, 0.0, -0.6)
    mol = api.new_molecule(env, ctypes.byref(n), at, xyz, ctypes.byref(q), ctypes.byref(uhf), None, None)
    calc = api.new_calculator()
    res = api.new_results()
    copied = None
    try:
        api.load_gfn2xtb(env, mol, calc, None)
        api.singlepoint(env, mol, calc, res)
        api.check(env, 'singlepoint')
        copied = api.copy_results(res)

        out1 = (ctypes.c_double * 3)()
        assert api.get_nuclear_esp(env, copied, 3, out1) == 1
        # Force SCF failure on res -> invalidates res but NOT copied
        api.set_max_iterations(env, calc, 1)
        api.singlepoint(env, mol, calc, res)
        assert api.check_environment(env)
        api.error(env)  # clear error stack

        sentinel = (ctypes.c_double * 3)(999.0, 999.0, 999.0)
        assert api.get_nuclear_esp(env, res, 3, sentinel) == 0
        assert list(sentinel) == [999.0, 999.0, 999.0]
        api.error(env)

        out2 = (ctypes.c_double * 3)()
        assert api.get_nuclear_esp(env, copied, 3, out2) == 1
        np.testing.assert_array_equal(list(out1), list(out2))
    finally:
        for fn, handle in (
            (api.del_results, copied),
            (api.del_results, res),
            (api.del_calculator, calc),
            (api.del_molecule, mol),
            (api.del_environment, env),
        ):
            api._delete(fn, handle)
