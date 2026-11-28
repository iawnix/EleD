"""DIPRO dimer-projection electronic coupling (transfer integral J_AB) descriptors.

Computes intermolecular hole (HOMO-HOMO), electron (LUMO-LUMO), and cross
charge-transfer (HOMO_A-LUMO_B, LUMO_A-HOMO_B) electronic couplings ``J_eff``
following Valeev (2006) and Grimme's xTB DIPRO formulation (Kohn et al.,
*J. Chem. Phys.* 2023, 159, 144106):

- Reconstructs the exact complex AO Fock/Hamiltonian matrix from the converged
  complex eigenbasis (``keepOverlap=True``):
  ``H_AO = (S_AB C_AB) diag(eps_AB) (S_AB C_AB)^T``.
- Embeds the isolated fragment A and B frontier MO subspaces into the dimer AO
  basis (verifying exact AO overlap alignment ``S_AB[ao_A, ao_A] == S_A`` and
  ``S_AB[ao_B, ao_B] == S_B``).
- Performs exact symmetric Lowdin orthogonalization ``H_eff = S_sub^(-1/2) H_sub
  S_sub^(-1/2)`` over the combined frontier subspace, which reduces identically
  to ``J_eff = (H_AB - 0.5*(e_A + e_B)*S_AB) / (1 - S_AB^2)`` for nondegenerate
  orbitals and remains strictly invariant to MO phase flips and unitary
  rotations within degenerate frontier subspaces (e.g. benzene dimers).

No new C ABI is required. Supports both ``native`` and ``compiled`` backends.
"""
from __future__ import annotations

from dataclasses import replace
import json
import math
from typing import Any

import numpy as np

from .electronic_descriptors import _as_result, _indices, _nonnegative
from .fragment_eda import _extract_submol, _run_single
from .input import _effective_charge_uhf, _xyz_from_mol
from .options import XTBOptions
from .orbital_descriptors import _overlap_basis
from .populations import _frontier_group

DIPRO_PROTOCOL = 'eled-dipro-v1'

__all__ = [
    'CalcElectronicCouplingFromResults',
    'CalcElectronicCoupling',
]


def _map_fragment_aos(
    ao_atom_ab: np.ndarray,
    frag_indices: np.ndarray,
    s_ab: np.ndarray,
    s_frag: np.ndarray,
    label: str,
) -> tuple[np.ndarray, float]:
    """Map isolated fragment AOs into the dimer AO basis and verify overlap alignment."""
    chunks = [np.flatnonzero(ao_atom_ab == int(atom_idx)) for atom_idx in frag_indices.tolist()]
    ao_in_ab = np.concatenate(chunks) if chunks else np.empty(0, dtype=int)
    if len(ao_in_ab) != s_frag.shape[0]:
        raise ValueError(
            f'{label} AO count ({s_frag.shape[0]}) does not match dimer AO subset ({len(ao_in_ab)})'
        )
    sub_s = s_ab[np.ix_(ao_in_ab, ao_in_ab)]
    err = float(np.max(np.abs(sub_s - s_frag)))
    if err > 1e-4:
        raise ValueError(
            f'{label} overlap matrix disagrees with the dimer overlap block (max error {err:.3e}); '
            'ensure fragment and complex use identical fixed coordinates, method, and basis'
        )
    return ao_in_ab, err


def _dipro_channel(
    h_ao_ev: np.ndarray,
    s_ab: np.ndarray,
    c_a: np.ndarray,
    c_b: np.ndarray,
    ao_a: np.ndarray,
    ao_b: np.ndarray,
    idx_a: list[int],
    idx_b: list[int],
    missing_label: str,
) -> dict[str, Any]:
    """Evaluate Lowdin-orthogonalized DIPRO coupling between MO subsets idx_a and idx_b."""
    if not idx_a or not idx_b:
        reason = f'{missing_label}_frontier_unavailable'
        return {
            'available': False,
            'orb_idx_a': list(idx_a),
            'orb_idx_b': list(idx_b),
            'n_orb_a': len(idx_a),
            'n_orb_b': len(idx_b),
            'j_eff_ev': None,
            'j_eff_mev': None,
            'j_frob_ev': None,
            'j_frob_mev': None,
            'j_signed_ev': None,
            'site_energy_a_ev': None,
            'site_energy_b_ev': None,
            'delta_e_site_ev': None,
            'overlap_rms': None,
            'coupling_matrix_ev': None,
            'missing': {
                'j_eff_ev': reason,
                'j_eff_mev': reason,
                'j_frob_ev': reason,
                'j_frob_mev': reason,
            },
        }

    nao_ab = s_ab.shape[0]
    na = len(idx_a)
    nb = len(idx_b)
    phi = np.zeros((nao_ab, na + nb), dtype=float)
    phi[ao_a, :na] = c_a[:, idx_a]
    phi[ao_b, na:] = c_b[:, idx_b]

    s_sub = phi.T @ s_ab @ phi
    s_sub = 0.5 * (s_sub + s_sub.T)
    h_sub = phi.T @ h_ao_ev @ phi
    h_sub = 0.5 * (h_sub + h_sub.T)

    evals, evecs = np.linalg.eigh(s_sub)
    if evals.min() <= 1e-10:
        reason = 'singular_interfragment_orbital_overlap'
        return {
            'available': False,
            'orb_idx_a': list(idx_a),
            'orb_idx_b': list(idx_b),
            'n_orb_a': na,
            'n_orb_b': nb,
            'j_eff_ev': None,
            'j_eff_mev': None,
            'j_frob_ev': None,
            'j_frob_mev': None,
            'j_signed_ev': None,
            'site_energy_a_ev': None,
            'site_energy_b_ev': None,
            'delta_e_site_ev': None,
            'overlap_rms': None,
            'coupling_matrix_ev': None,
            'missing': {'j_eff_ev': reason, 'j_eff_mev': reason},
        }

    s_inv_half = (evecs * (1.0 / np.sqrt(evals))[None, :]) @ evecs.T
    h_eff = s_inv_half @ h_sub @ s_inv_half
    h_eff = 0.5 * (h_eff + h_eff.T)

    j_block = h_eff[:na, na:]
    s_block = s_sub[:na, na:]
    ea_eff = np.diag(h_eff[:na, :na])
    eb_eff = np.diag(h_eff[na:, na:])

    frob = float(np.linalg.norm(j_block, ord='fro'))
    j_eff = float(math.sqrt(float(np.sum(j_block * j_block)) / float(na * nb)))
    s_rms = float(math.sqrt(float(np.sum(s_block * s_block)) / float(na * nb)))
    ea_mean = float(np.mean(ea_eff))
    eb_mean = float(np.mean(eb_eff))

    missing: dict[str, str] = {}
    if na == 1 and nb == 1:
        j_signed: float | None = float(j_block[0, 0])
    else:
        j_signed = None
        missing['j_signed_ev'] = 'degenerate_subspace_sign_arbitrary'

    return {
        'available': True,
        'orb_idx_a': list(idx_a),
        'orb_idx_b': list(idx_b),
        'n_orb_a': na,
        'n_orb_b': nb,
        'j_eff_ev': j_eff,
        'j_eff_mev': float(j_eff * 1000.0),
        'j_frob_ev': frob,
        'j_frob_mev': float(frob * 1000.0),
        'j_signed_ev': j_signed,
        'site_energy_a_ev': ea_mean,
        'site_energy_b_ev': eb_mean,
        'delta_e_site_ev': float(eb_mean - ea_mean),
        'overlap_rms': s_rms,
        'coupling_matrix_ev': j_block.tolist(),
        'missing': missing,
    }


def CalcElectronicCouplingFromResults(
    mol: Any,
    complexResult: Any,
    fragmentAResult: Any,
    fragmentBResult: Any,
    fragmentA: Any,
    fragmentB: Any,
    *,
    degeneracyEV: float = 1e-3,
    occupationTolerance: float = 1e-6,
) -> dict[str, Any]:
    """Compute DIPRO electronic couplings from three retained ``keepOverlap=True`` results.

    Never invokes xTB. Validates exact AO overlap block agreement between the
    complex and isolated fragments.
    """
    _nonnegative(degeneracyEV, 'degeneracyEV')
    _nonnegative(occupationTolerance, 'occupationTolerance')

    natoms = int(mol.GetNumAtoms())
    idx_a = _indices(fragmentA, natoms, 'fragmentA')
    idx_b = _indices(fragmentB, natoms, 'fragmentB')
    if set(idx_a.tolist()).intersection(idx_b.tolist()):
        raise ValueError('fragmentA and fragmentB must be disjoint')

    res_ab = _as_result(complexResult)
    res_a = _as_result(fragmentAResult)
    res_b = _as_result(fragmentBResult)

    e_ab, _, _, _, c_ab, s_ab, ao_atom_ab, ortho_ab = _overlap_basis(
        res_ab
    )
    e_a, occ_a, occ_idx_a, vir_idx_a, c_a, s_a, _, ortho_a = _overlap_basis(
        res_a
    )
    e_b, occ_b, occ_idx_b, vir_idx_b, c_b, s_b, _, ortho_b = _overlap_basis(
        res_b
    )

    ao_a, err_a = _map_fragment_aos(ao_atom_ab, idx_a, s_ab, s_a, 'fragmentA')
    ao_b, err_b = _map_fragment_aos(ao_atom_ab, idx_b, s_ab, s_b, 'fragmentB')

    # Exact AO Hamiltonian reconstruction in eV from complex eigenbasis.
    sc_ab = s_ab @ c_ab
    h_ao_ev = (sc_ab * e_ab[None, :]) @ sc_ab.T
    h_ao_ev = 0.5 * (h_ao_ev + h_ao_ev.T)

    homo_a = _frontier_group(e_a, occ_a, occ_idx_a, vir_idx_a, True, degeneracyEV, occupationTolerance)
    lumo_a = _frontier_group(e_a, occ_a, occ_idx_a, vir_idx_a, False, degeneracyEV, occupationTolerance)
    homo_b = _frontier_group(e_b, occ_b, occ_idx_b, vir_idx_b, True, degeneracyEV, occupationTolerance)
    lumo_b = _frontier_group(e_b, occ_b, occ_idx_b, vir_idx_b, False, degeneracyEV, occupationTolerance)

    hole = _dipro_channel(h_ao_ev, s_ab, c_a, c_b, ao_a, ao_b, homo_a, homo_b, 'hole')
    electron = _dipro_channel(h_ao_ev, s_ab, c_a, c_b, ao_a, ao_b, lumo_a, lumo_b, 'electron')
    ct_ab = _dipro_channel(h_ao_ev, s_ab, c_a, c_b, ao_a, ao_b, homo_a, lumo_b, 'ct_a_to_b')
    ct_ba = _dipro_channel(h_ao_ev, s_ab, c_a, c_b, ao_a, ao_b, lumo_a, homo_b, 'ct_b_to_a')

    missing: dict[str, str] = {}
    for key, ch in (
        ('j_hole_ev', hole),
        ('j_electron_ev', electron),
        ('j_ct_ab_ev', ct_ab),
        ('j_ct_ba_ev', ct_ba),
    ):
        if not ch['available']:
            missing[key] = ch['missing'].get('j_eff_ev', 'channel_unavailable')

    payload = {
        'j_hole_ev': hole['j_eff_ev'],
        'j_hole_mev': hole['j_eff_mev'],
        'j_electron_ev': electron['j_eff_ev'],
        'j_electron_mev': electron['j_eff_mev'],
        'j_ct_ab_ev': ct_ab['j_eff_ev'],
        'j_ct_ab_mev': ct_ab['j_eff_mev'],
        'j_ct_ba_ev': ct_ba['j_eff_ev'],
        'j_ct_ba_mev': ct_ba['j_eff_mev'],
        'hole': hole,
        'electron': electron,
        'ct_a_to_b': ct_ab,
        'ct_b_to_a': ct_ba,
        'basis_overlap_error': max(err_a, err_b),
        'orthonormality_error': max(ortho_ab, ortho_a, ortho_b),
        'method': res_ab.method,
        'backend': res_ab.backend,
        'fragments': {
            'fragment_a': idx_a.tolist(),
            'fragment_b': idx_b.tolist(),
        },
        'units': {
            'j_hole_ev': 'eV',
            'j_hole_mev': 'meV',
            'j_electron_ev': 'eV',
            'j_electron_mev': 'meV',
            'j_ct_ab_ev': 'eV',
            'j_ct_ab_mev': 'meV',
            'j_ct_ba_ev': 'eV',
            'j_ct_ba_mev': 'meV',
            'j_eff_ev': 'eV',
            'j_eff_mev': 'meV',
            'j_frob_ev': 'eV',
            'j_frob_mev': 'meV',
            'j_signed_ev': 'eV',
            'site_energy_a_ev': 'eV',
            'site_energy_b_ev': 'eV',
            'delta_e_site_ev': 'eV',
            'overlap_rms': 'dimensionless',
        },
        'analysis_settings': {
            'protocol': DIPRO_PROTOCOL,
            'reference': 'Valeev (2006); Kohn, Gildemeister, Fazzi & Hansen, J. Chem. Phys. 2023, 159, 144106',
            'degeneracy_ev': float(degeneracyEV),
            'frontier_definition': 'aufbau electron count',
            'occupation_tolerance': float(occupationTolerance),
            'definition': 'H_eff = S_sub^(-1/2) H_sub S_sub^(-1/2); j_eff = sqrt(sum(J_ab^2) / (n_a * n_b))',
        },
        'invariance': (
            'j_eff, j_frob, and overlap_rms are invariant to rigid 3D motion, MO phase flips, '
            'degenerate frontier subspace rotations, and swapping fragmentA/fragmentB'
        ),
        'missing': missing,
    }
    return json.loads(json.dumps(payload, allow_nan=False))


def CalcElectronicCoupling(
    mol: Any,
    fragmentA: Any,
    fragmentB: Any,
    *,
    chargeA: int | None = None,
    chargeB: int | None = None,
    unpairedA: int | None = None,
    unpairedB: int | None = None,
    allowCovalentCut: bool = False,
    degeneracyEV: float = 1e-3,
    occupationTolerance: float = 1e-6,
    options: XTBOptions | None = None,
    interface: str = 'native',
) -> dict[str, Any]:
    """Run three matched xTB single points (AB, A, B) and compute DIPRO couplings J_AB."""
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
            'pass allowCovalentCut=True with explicit fragment charges/spins if intended'
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

    report = CalcElectronicCouplingFromResults(
        snapshot,
        raw_ab,
        raw_a,
        raw_b,
        idx_a,
        idx_b,
        degeneracyEV=degeneracyEV,
        occupationTolerance=occupationTolerance,
    )
    report['analysis_settings']['interface'] = interface
    report['analysis_settings']['conformer_id'] = conf_id
    return report
