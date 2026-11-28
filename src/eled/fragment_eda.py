"""Fragment interaction energy decomposition (EDA) and charge-transfer descriptors.

Evaluates fixed-geometry intermolecular or fragment-pair electronic interaction
descriptors across three matched xTB single points (``complex``, ``fragment_a``,
and ``fragment_b``) with retained overlap basis (``keepOverlap=True``):

- **Energy decomposition** (``delta_X = X(AB) - X(A) - X(B)``, in Hartree and
  kcal/mol): ``delta_e_total``, ``delta_e_elec``, ``delta_e_rep``,
  ``delta_e_es``, ``delta_e_aes``, ``delta_e_axc``, ``delta_e_disp``,
  ``delta_e_xb``, and ``delta_g_solv``. Method-inapplicable terms remain
  ``None`` with explicit ``missing`` reasons.
- **Charge transfer** (electrons transferred from fragment A to fragment B;
  positive means A donates electron density to B): ``ct_scc_a_to_b_e``,
  ``ct_cm5_a_to_b_e`` (GFN1/2 only), ``ct_lowdin_a_to_b_e``, and
  ``ct_mulliken_a_to_b_e``.
- **Inter-fragment complex bonding & response**: ``wbo_inter``,
  ``coop_inter_e`` (positive = net bonding), ``ecohp_inter_ev`` (negative =
  bonding stabilization), closed-shell ``chi_inter_e_per_ev`` (inter-fragment
  sum of independent-particle linear response), and ``delta_dipole_debye``.
"""
from __future__ import annotations

from dataclasses import replace
import json
from typing import Any

import numpy as np

from .calculator import XTBCalculator
from .electronic_descriptors import (
    AU_DIPOLE_TO_DEBYE,
    _as_result,
    _charges,
    _cm5_charges,
    _indices,
    _wbo,
)
from .errors import XTBEnvironmentError, XTBRunError
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .orbital_descriptors import CalcLinearResponseMatrix, _overlap_basis
from .populations import CalcOrbitalPopulations
from .properties import HARTREE_TO_KCAL_MOL, CalcEnergyDecomposition

FRAGMENT_EDA_PROTOCOL = 'eled-fragment-eda-v1'

__all__ = [
    'CalcFragmentInteraction',
]


def _extract_submol(mol: Any, conf_id: int, indices: np.ndarray) -> Any:
    """Build an immutable explicit-H RDKit sub-molecule at exact fixed coordinates."""
    from rdkit import Chem

    orig_conf = mol.GetConformer(conf_id)
    index_set = set(indices.tolist())
    old_to_new = {int(old): new for new, old in enumerate(indices.tolist())}
    rw = Chem.RWMol()
    for old_idx in indices.tolist():
        atom = mol.GetAtomWithIdx(int(old_idx))
        copy_atom = Chem.Atom(int(atom.GetAtomicNum()))
        copy_atom.SetFormalCharge(int(atom.GetFormalCharge()))
        copy_atom.SetNumRadicalElectrons(int(atom.GetNumRadicalElectrons()))
        copy_atom.SetIsotope(int(atom.GetIsotope()))
        copy_atom.SetNoImplicit(True)
        copy_atom.SetNumExplicitHs(0)
        rw.AddAtom(copy_atom)
    for bond in mol.GetBonds():
        a = int(bond.GetBeginAtomIdx())
        b = int(bond.GetEndAtomIdx())
        if a in index_set and b in index_set:
            rw.AddBond(old_to_new[a], old_to_new[b], bond.GetBondType())
    sub = rw.GetMol()
    sub.UpdatePropertyCache(strict=False)
    conf = Chem.Conformer(len(indices))
    conf.Set3D(True)
    conf.SetId(0)
    for new_i, old_i in enumerate(indices.tolist()):
        pos = orig_conf.GetAtomPosition(int(old_i))
        conf.SetAtomPosition(new_i, pos)
    sub.AddConformer(conf, assignId=False)
    return sub


def _run_single(mol: Any, opts: XTBOptions, interface: str):
    if interface == 'native':
        return XTBCalculator(**opts.ToDict()).Calculate(mol)
    from . import rdXTBTools

    ok, raw = rdXTBTools.RunMol(
        mol,
        confId=opts.confId,
        method=opts.method,
        charge=opts.charge,
        numUnpairedElectrons=opts.numUnpairedElectrons,
        keepOrbitalCoefficients=True,
        xtbLibrary=opts.xtbLibrary or '',
        maxIterations=opts.maxIterations,
        accuracy=opts.accuracy,
        electronicTemperature=opts.electronicTemperature,
        keepOverlap=True,
        solvent=opts.solvent,
        solvationModel=opts.solvationModel,
        solvationGrid=opts.solvationGrid,
    )
    if not ok:
        error_type = XTBRunError if raw.failureKind == 'calculation' else XTBEnvironmentError
        raise error_type(raw.errorMessage)
    return raw


def _delta_opt(c_val: float | None, a_val: float | None, b_val: float | None) -> float | None:
    if c_val is None or a_val is None or b_val is None:
        return None
    return float(c_val - a_val - b_val)


def CalcFragmentInteraction(
    mol: Any,
    fragmentA: Any,
    fragmentB: Any,
    *,
    chargeA: int | None = None,
    chargeB: int | None = None,
    unpairedA: int | None = None,
    unpairedB: int | None = None,
    allowCovalentCut: bool = False,
    occupationThreshold: float = 1e-7,
    options: XTBOptions | None = None,
    interface: str = 'native',
) -> dict[str, Any]:
    """Compute fragment interaction energy decomposition, charge transfer, and inter-fragment bonding.

    Runs three matched fixed-geometry xTB single points (complex ``AB``,
    isolated ``A``, isolated ``B``) with ``keepOverlap=True``. By default
    (``allowCovalentCut=False``), rejects RDKit bonds crossing between
    ``fragmentA`` and ``fragmentB``; pass ``allowCovalentCut=True`` with
    explicit fragment charges/spins for covalent bond-cleavage EDA.
    """
    from rdkit import Chem

    if interface not in ('native', 'compiled'):
        raise ValueError("interface must be 'native' or 'compiled'")
    if not isinstance(allowCovalentCut, bool):
        raise ValueError('allowCovalentCut must be bool')
    if options is not None and not isinstance(options, XTBOptions):
        raise TypeError('options must be an XTBOptions instance or None')

    snapshot = Chem.Mol(mol)
    base = options or XTBOptions()
    _xyz_from_mol(snapshot, base.confId)
    natoms = snapshot.GetNumAtoms()
    conf_id = int(snapshot.GetConformer(base.confId).GetId())

    idx_a = _indices(fragmentA, natoms, 'fragmentA')
    idx_b = _indices(fragmentB, natoms, 'fragmentB')
    set_a = set(idx_a.tolist())
    set_b = set(idx_b.tolist())
    if set_a.intersection(set_b):
        raise ValueError('fragmentA and fragmentB must be disjoint')

    crossing = [
        (int(b.GetBeginAtomIdx()), int(b.GetEndAtomIdx()))
        for b in snapshot.GetBonds()
        if (int(b.GetBeginAtomIdx()) in set_a and int(b.GetEndAtomIdx()) in set_b)
        or (int(b.GetBeginAtomIdx()) in set_b and int(b.GetEndAtomIdx()) in set_a)
    ]
    if crossing and not allowCovalentCut:
        raise ValueError(
            f'fragmentA and fragmentB are linked by {len(crossing)} covalent bond(s) {crossing[:3]}; '
            'pass allowCovalentCut=True with explicit chargeA/chargeB/unpairedA/unpairedB for covalent cleavage EDA'
        )

    full_partition = (len(set_a) + len(set_b) == natoms)
    mol_a = _extract_submol(snapshot, conf_id, idx_a)
    mol_b = _extract_submol(snapshot, conf_id, idx_b)

    q_ab, u_ab = _effective_charge_uhf(snapshot, base.charge, base.numUnpairedElectrons)
    q_a, u_a = _effective_charge_uhf(mol_a, chargeA, unpairedA)
    q_b, u_b = _effective_charge_uhf(mol_b, chargeB, unpairedB)
    if full_partition and (q_a + q_b != q_ab):
        raise ValueError(
            f'fragment charges ({q_a} + {q_b}) must sum to complex charge ({q_ab}) for a full partition'
        )

    opts_ab = replace(
        base, confId=conf_id, charge=q_ab, numUnpairedElectrons=u_ab, keepOverlap=True,
        espPointsAngstrom=(), realSpacePointsAngstrom=(),
    )
    opts_a = replace(
        base, confId=0, charge=q_a, numUnpairedElectrons=u_a, keepOverlap=True,
        espPointsAngstrom=(), realSpacePointsAngstrom=(),
    )
    opts_b = replace(
        base, confId=0, charge=q_b, numUnpairedElectrons=u_b, keepOverlap=True,
        espPointsAngstrom=(), realSpacePointsAngstrom=(),
    )

    raw_ab = _run_single(snapshot, opts_ab, interface)
    raw_a = _run_single(mol_a, opts_a, interface)
    raw_b = _run_single(mol_b, opts_b, interface)

    res_ab = _as_result(raw_ab)
    res_a = _as_result(raw_a)
    res_b = _as_result(raw_b)

    # 1. Energy decomposition delta = AB - A - B.
    dec_ab = CalcEnergyDecomposition(raw_ab)
    dec_a = CalcEnergyDecomposition(raw_a)
    dec_b = CalcEnergyDecomposition(raw_b)

    energy_keys = (
        ('delta_e_total', 'e_total_ha'),
        ('delta_e_elec', 'e_elec_ha'),
        ('delta_e_rep', 'e_rep_ha'),
        ('delta_e_es', 'e_es_ha'),
        ('delta_e_aes', 'e_aes_ha'),
        ('delta_e_axc', 'e_axc_ha'),
        ('delta_e_disp', 'e_disp_ha'),
        ('delta_e_xb', 'e_xb_ha'),
        ('delta_g_solv', 'g_solv_ha'),
    )
    energies: dict[str, float | None] = {}
    missing: dict[str, str] = {}
    for out_prefix, src_key in energy_keys:
        if src_key == 'e_total_ha':
            val = _delta_opt(res_ab.totalEnergyHartree, res_a.totalEnergyHartree, res_b.totalEnergyHartree)
        else:
            val = _delta_opt(dec_ab.get(src_key), dec_a.get(src_key), dec_b.get(src_key))
        energies[f'{out_prefix}_ha'] = val
        energies[f'{out_prefix}_kcalmol'] = float(val * HARTREE_TO_KCAL_MOL) if val is not None else None
        if val is None:
            reason = (
                dec_ab['missing'].get(src_key)
                or dec_a['missing'].get(src_key)
                or dec_b['missing'].get(src_key)
                or 'energy_component_unavailable'
            )
            missing[f'{out_prefix}_ha'] = reason
            missing[f'{out_prefix}_kcalmol'] = reason

    # 2. Charge transfer from fragment A to fragment B (positive = A donates electrons to B).
    q_comp = _charges(res_ab)
    q_iso_a = _charges(res_a)
    ct_scc = float(np.sum(q_comp[idx_a]) - np.sum(q_iso_a))

    cm5_comp = _cm5_charges(res_ab)
    cm5_iso_a = _cm5_charges(res_a)
    if len(cm5_comp) and len(cm5_iso_a):
        ct_cm5: float | None = float(np.sum(cm5_comp[idx_a]) - np.sum(cm5_iso_a))
    else:
        ct_cm5 = None
        missing['ct_cm5_a_to_b_e'] = (
            'cm5_unavailable_for_gfn0' if res_ab.method == 'GFN0-xTB' else 'cm5_charges_unavailable'
        )

    low_ab = np.asarray(CalcOrbitalPopulations(res_ab, scheme='lowdin')['mo_n_e'], dtype=float)
    low_a = np.asarray(CalcOrbitalPopulations(res_a, scheme='lowdin')['mo_n_e'], dtype=float)
    mul_ab = np.asarray(CalcOrbitalPopulations(res_ab, scheme='mulliken')['mo_n_e'], dtype=float)
    mul_a = np.asarray(CalcOrbitalPopulations(res_a, scheme='mulliken')['mo_n_e'], dtype=float)
    # Electron population lost by fragment A = N_A(isolated) - N_A(complex).
    ct_lowdin = float(np.sum(low_a) - np.sum(low_ab[idx_a]))
    ct_mulliken = float(np.sum(mul_a) - np.sum(mul_ab[idx_a]))

    # 3. Inter-fragment bonding in complex (WBO, COOP, ECOHP, condensed linear response, dipole shift).
    wbo_mat = _wbo(res_ab, natoms)
    wbo_inter = float(np.sum(wbo_mat[np.ix_(idx_a, idx_b)])) if wbo_mat is not None else None
    if wbo_inter is None:
        missing['wbo_inter'] = 'wbo_unavailable'

    e_mo, occ_mo, occupied, _, c_mat, s_mat, ao_atom, ortho_err = _overlap_basis(
        res_ab, occupationThreshold=occupationThreshold
    )
    if len(occupied):
        c_occ = c_mat[:, occupied]
        n_occ = occ_mo[occupied]
        en_occ = e_mo[occupied]
        p_mat = (c_occ * n_occ[None, :]) @ c_occ.T
        w_mat = (c_occ * (n_occ * en_occ)[None, :]) @ c_occ.T
        ps = p_mat * s_mat
        ws = w_mat * s_mat
        ao_a = np.flatnonzero(np.isin(ao_atom, idx_a))
        ao_b = np.flatnonzero(np.isin(ao_atom, idx_b))
        coop_inter: float | None = float(2.0 * np.sum(ps[np.ix_(ao_a, ao_b)]))
        ecohp_inter: float | None = float(2.0 * np.sum(ws[np.ix_(ao_a, ao_b)]))
    else:
        coop_inter = ecohp_inter = None
        missing['coop_inter_e'] = 'no_occupied_orbitals'
        missing['ecohp_inter_ev'] = 'no_occupied_orbitals'

    lr = CalcLinearResponseMatrix(res_ab, scheme='lowdin', occupationThreshold=occupationThreshold)
    if lr['available']:
        chi_mat = np.asarray(lr['chi_e_per_ev'], dtype=float)
        chi_inter: float | None = float(np.sum(chi_mat[np.ix_(idx_a, idx_b)]))
    else:
        chi_inter = None
        missing['chi_inter_e_per_ev'] = lr['missing'].get('chi_e_per_ev', 'linear_response_unavailable')

    dip_ab = np.asarray(res_ab.dipoleAtomicUnits, dtype=float)
    dip_a = np.asarray(res_a.dipoleAtomicUnits, dtype=float)
    dip_b = np.asarray(res_b.dipoleAtomicUnits, dtype=float)
    if len(dip_ab) == 3 and len(dip_a) == 3 and len(dip_b) == 3:
        delta_dip_au = float(np.linalg.norm(dip_ab - dip_a - dip_b))
        delta_dip_debye: float | None = float(delta_dip_au * AU_DIPOLE_TO_DEBYE)
    else:
        delta_dip_au = None
        delta_dip_debye = None
        missing['delta_dipole_debye'] = 'dipole_unavailable'

    payload = {
        **energies,
        'ct_scc_a_to_b_e': ct_scc,
        'ct_cm5_a_to_b_e': ct_cm5,
        'ct_lowdin_a_to_b_e': ct_lowdin,
        'ct_mulliken_a_to_b_e': ct_mulliken,
        'wbo_inter': wbo_inter,
        'coop_inter_e': coop_inter,
        'ecohp_inter_ev': ecohp_inter,
        'chi_inter_e_per_ev': chi_inter,
        'delta_dipole_ebohr': delta_dip_au,
        'delta_dipole_debye': delta_dip_debye,
        'orthonormality_error': ortho_err,
        'fragments': {
            'fragment_a': {'atom_indices': idx_a.tolist(), 'charge': int(q_a), 'uhf': int(u_a)},
            'fragment_b': {'atom_indices': idx_b.tolist(), 'charge': int(q_b), 'uhf': int(u_b)},
            'complex': {'charge': int(q_ab), 'uhf': int(u_ab)},
            'full_partition': full_partition,
            'covalent_bonds_cut': len(crossing),
        },
        'method': res_ab.method,
        'backend': res_ab.backend,
        'units': {
            **{f'{p}_ha': 'Hartree' for p, _ in energy_keys},
            **{f'{p}_kcalmol': 'kcal/mol' for p, _ in energy_keys},
            'ct_scc_a_to_b_e': 'e',
            'ct_cm5_a_to_b_e': 'e',
            'ct_lowdin_a_to_b_e': 'e',
            'ct_mulliken_a_to_b_e': 'e',
            'wbo_inter': 'dimensionless',
            'coop_inter_e': 'e',
            'ecohp_inter_ev': 'eV',
            'chi_inter_e_per_ev': 'e/eV',
            'delta_dipole_ebohr': 'e*Bohr',
            'delta_dipole_debye': 'Debye',
        },
        'analysis_settings': {
            'protocol': FRAGMENT_EDA_PROTOCOL,
            'interface': interface,
            'allow_covalent_cut': allowCovalentCut,
            'occupation_threshold': float(occupationThreshold),
            'charge_transfer_sign': 'positive means electrons transferred from fragmentA to fragmentB',
            'coop_sign': 'positive for net inter-fragment bonding, negative for antibonding/Pauli repulsion',
            'ecohp_sign': 'negative for inter-fragment bonding stabilization, positive for destabilization',
            'dipole_origin_note': 'delta_dipole is origin-independent only when both fragments are neutral',
            'conformer_id': conf_id,
            'solvent': res_ab.solvent,
            'solvation_model': res_ab.solvationModel,
        },
        'invariance': (
            'interaction energy components, charge transfer, WBO, COOP, ECOHP, and chi_inter '
            'are rigid-motion, AO-permutation, and MO-phase invariant'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))
