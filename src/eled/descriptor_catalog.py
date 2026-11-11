"""Canonical, machine-readable descriptor definitions.

The catalog is intentionally descriptive: it does not calculate values.  A
field's ``missing`` text is part of its public contract and means callers must
preserve ``None``/the reason instead of imputing zero.
"""
from __future__ import annotations

from copy import deepcopy

_CATALOG = {
    'molecular': {
        'xtb_energy_ha': {'definition': 'xTB total electronic energy', 'unit': 'Hartree', 'invariance': 'rigid-motion invariant; method/reference dependent', 'requires': 'single point', 'missing': 'engine did not return energy'},
        'xtb_homo_ev': {'definition': 'highest orbital with occupation above threshold', 'unit': 'eV', 'invariance': 'orbital phase and rigid-motion invariant', 'requires': 'orbital energies and occupations', 'missing': 'no occupied orbital under threshold'},
        'xtb_lumo_ev': {'definition': 'lowest orbital at or below occupation threshold', 'unit': 'eV', 'invariance': 'orbital phase and rigid-motion invariant', 'requires': 'orbital energies and occupations', 'missing': 'no virtual orbital under threshold'},
        'xtb_gap_ev': {'definition': 'LUMO minus HOMO', 'unit': 'eV', 'invariance': 'same as frontier energies', 'requires': 'HOMO and LUMO', 'missing': 'HOMO or LUMO unavailable'},
        'xtb_dip_debye': {'definition': 'norm of the xTB dipole', 'unit': 'Debye', 'invariance': 'rigid rotation invariant; origin dependent for charged systems', 'requires': 'dipole vector', 'missing': 'dipole not returned'},
        'xtb_dip_*_ebohr': {'definition': 'Cartesian xTB dipole component', 'unit': 'e*Bohr', 'invariance': 'components rotate with the coordinate frame; origin dependent for charged systems', 'requires': 'dipole vector', 'missing': 'dipole not returned'},
        'homo_* / lumo_*': {'definition': 'statistics of the degeneracy-grouped HOMO or LUMO atom distribution', 'unit': 'field dependent (eV, count, dimensionless or Angstrom)', 'invariance': 'rigid-motion invariant except centroid components; invariant to phase and rotations inside the selected degenerate group', 'requires': 'retained overlap and orbital basis', 'missing': 'frontier group or basis unavailable'},
        'xtb_q_*_e': {'definition': 'summary statistics of xTB atomic partial charges', 'unit': 'e', 'invariance': 'atom reordering invariant for aggregate statistics', 'requires': 'partial charges', 'missing': 'charges not returned'},
        'eht_*': {'definition': 'YAeHMOP scalar orbital features', 'unit': 'field suffix', 'invariance': 'molecule invariant at fixed input geometry', 'requires': 'RDKit rdEHTTools result', 'missing': 'EHT calculation failed'},
        'xtb_wbo_sum': {'definition': 'sum of all off-diagonal Wiberg pairs', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'xtb_n_*': {'definition': 'reported orbital/electron/unpaired-electron count', 'unit': 'count', 'invariance': 'molecule invariant', 'requires': 'engine result', 'missing': 'engine did not return count'},
        'xtb_*_idx': {'definition': 'zero-based orbital index under the recorded occupation convention', 'unit': 'index', 'invariance': 'depends on the retained orbital ordering', 'requires': 'orbital energies and occupations', 'missing': 'frontier unavailable'},
        'xtb_occ_cut': {'definition': 'occupation threshold used to define occupied/virtual orbitals', 'unit': 'dimensionless', 'invariance': 'protocol setting', 'requires': 'descriptor extraction', 'missing': 'never missing for a valid extraction'},
    },
    'atom': {
        'q_e': {'definition': 'xTB engine atomic partial charge', 'unit': 'e', 'invariance': 'rigid-motion invariant; atom index follows input order', 'requires': 'partial charges', 'missing': 'charges not returned'},
        'wbo_sum': {'definition': 'sum of off-diagonal Wiberg orders incident on atom', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'mul_qrec_e / low_qrec_e': {'definition': 'reference valence electrons minus Mulliken/Löwdin shell population', 'unit': 'e', 'invariance': 'AO permutation and orbital phase invariant with maps', 'requires': 'same-calculation density and supported element', 'missing': 'unsupported element, GFN0 convention, or basis unavailable'},
        'mul_*_e': {'definition': 'Mulliken density population by angular momentum shell', 'unit': 'e', 'invariance': 'AO permutation and orbital phase invariant with maps', 'requires': 'same-calculation P, S and shell map', 'missing': 'overlap/basis extension unavailable'},
        'low_*_e': {'definition': 'Löwdin density population by angular momentum shell', 'unit': 'e', 'invariance': 'AO permutation and orbital phase invariant with maps', 'requires': 'same-calculation P, S and shell map', 'missing': 'overlap/basis extension unavailable'},
        'mul_ps / low_ps': {'definition': 'p-shell population divided by s-shell population', 'unit': 'dimensionless', 'invariance': 'same as shell population', 'requires': 's and p shell populations; nonzero s population', 'missing': 's population unavailable or zero'},
        'pi_n_e': {'definition': 'p-shell density projected onto local plane normal', 'unit': 'e', 'invariance': 'rigid-motion and normal-sign invariant', 'requires': 'sp2/aromatic atom, noncollinear planar geometry, one p shell and AO component directions', 'missing': 'atom not sp2/aromatic or normal/p map undefined'},
        'pi_q_e': {'definition': 'reference pi electrons minus pi population; default 1 e only for neutral nonradical sp2 carbon', 'unit': 'e', 'invariance': 'same as pi population', 'requires': 'pi population and defined pi reference', 'missing': 'undefined reference or pi projection'},
        'f_plus/f_minus/f_radical/dual': {'definition': 'condensed finite-difference Fukui indices', 'unit': 'dimensionless', 'invariance': 'fixed atom indexing and geometry; rigid-motion invariant', 'requires': 'compatible Q, Q+1 and Q-1 results', 'missing': 'required charge state failed'},
        's_loc_1ev': {'definition': 'global softness times f_plus', 'unit': '1/eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and f_plus', 'missing': 'hardness undefined or charge state failed'},
        'omega_loc_ev': {'definition': 'global electrophilicity times f_plus', 'unit': 'eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and f_plus', 'missing': 'electrophilicity undefined or charge state failed'},
        'n_loc_1ev': {'definition': 'softness-scaled f_minus proxy; no absolute nucleophilicity reference assumed', 'unit': '1/eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and f_minus', 'missing': 'softness undefined or charge state failed'},
    },
    'fragment': {
        'q_e': {'definition': 'sum of atomic partial charges in the named fragment', 'unit': 'e', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'partial charges and explicit fragment indices', 'missing': 'charges not returned'},
        'q_min_e / q_max_e / q_sd_e': {'definition': 'minimum, maximum and population standard deviation of fragment atomic charges', 'unit': 'e', 'invariance': 'atom reordering invariant within the fragment', 'requires': 'partial charges and explicit fragment indices', 'missing': 'charges not returned'},
        'wbo_intra / wbo_inter': {'definition': 'within-fragment and fragment-to-outside Wiberg sums', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix and explicit fragment indices', 'missing': 'bond orders not returned'},
    },
    'bond': {
        'wbo': {'definition': 'xTB Wiberg bond order on an RDKit bond', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'mbo': {'definition': 'closed-shell Mayer order from P S products', 'unit': 'dimensionless', 'invariance': 'AO permutation and orbital phase invariant', 'requires': 'same-calculation P/S and closed-shell density', 'missing': 'open-shell spin density unavailable or basis missing'},
        'mul_ovlp_e': {'definition': 'signed Mulliken bond overlap of the one-electron-normalized HOMO/LUMO group projector', 'unit': 'e', 'invariance': 'frontier-group rotations and orbital phase invariant', 'requires': 'frontier group and same-calculation P/S', 'missing': 'no frontier group or basis missing'},
        'pi_bo': {'definition': 'closed-shell projected Mayer (Mulliken scheme) or squared orthogonal-density (Lowdin scheme) pi bond index', 'unit': 'dimensionless', 'invariance': 'rigid-motion and normal-sign invariant', 'requires': 'closed shell and two valid local p projections', 'missing': 'open/unknown spin state or missing endpoint projection'},
    },
    'density': {
        'density_matrix': {'definition': 'spin-summed AO density P = C diag(n) C.T', 'unit': 'basis-dependent; trace(P S) is e', 'invariance': 'matrix transforms with AO basis; trace with S and MO phase are invariant', 'requires': 'same-calculation C, S and occupations', 'missing': 'overlap/basis extension unavailable'},
        'tr_ps_e': {'definition': 'electron count from trace(P S)', 'unit': 'e', 'invariance': 'AO permutation invariant with transformed S', 'requires': 'density matrix and overlap', 'missing': 'basis data unavailable'},
        'n_e': {'definition': 'reported modeled electron count; occupation sum if absent', 'unit': 'e', 'invariance': 'AO and MO representation invariant', 'requires': 'retained density result', 'missing': 'basis data unavailable'},
        'n_err_e': {'definition': 'trace(P S) minus reported modeled electron count (occupation sum if absent)', 'unit': 'e', 'invariance': 'AO permutation invariant', 'requires': 'density matrix and occupations', 'missing': 'basis data unavailable'},
        'shell_err_e': {'definition': 'maximum shell-to-atom population sum discrepancy', 'unit': 'e', 'invariance': 'AO permutation invariant with maps', 'requires': 'shell map and density', 'missing': 'shell map unavailable'},
    },
    'reactivity': {
        'ip_v_ev': {'definition': 'vertical ionization energy from E(Q+1)-E(Q)', 'unit': 'eV', 'invariance': 'rigid-motion invariant at fixed geometry', 'requires': 'successful Q and Q+1 GFN1/2 calculations', 'missing': 'required charge state failed'},
        'ea_v_ev': {'definition': 'vertical electron affinity from E(Q)-E(Q-1)', 'unit': 'eV', 'invariance': 'rigid-motion invariant at fixed geometry', 'requires': 'successful Q and Q-1 GFN1/2 calculations', 'missing': 'required charge state failed'},
        'mu_ev': {'definition': 'negative electronegativity, -(IP+EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        'chi_ev': {'definition': 'electronegativity, (IP+EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        'eta_ev': {'definition': 'conceptual-DFT hardness, (IP-EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        's_1ev': {'definition': 'global softness, 1/(IP-EA)', 'unit': '1/eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'omega_ev': {'definition': 'global electrophilicity, mu^2/(2 hardness)', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'f_plus/f_minus/f_radical/dual': {'definition': 'condensed finite-difference Fukui indices per atom', 'unit': 'dimensionless', 'invariance': 'fixed atom indexing and geometry; rigid-motion invariant', 'requires': 'compatible Q, Q+1 and Q-1 results', 'missing': 'required charge state failed'},
    },
    'diagnostic': {
        'conformer_dispersion': {'definition': 'population standard deviation across supplied conformers', 'unit': 'descriptor-specific', 'invariance': 'independent of conformer ordering', 'requires': 'at least one valid conformer report; a single conformer gives zero spread', 'missing': 'no valid conformers'},
        'gfn12_delta': {'definition': 'absolute difference between method reports', 'unit': 'descriptor-specific', 'invariance': 'same atom/bond indexing and geometry protocol', 'requires': 'both method reports', 'missing': 'one method report or value missing'},
        'et300_1000_delta': {'definition': 'absolute difference under electronic temperature perturbation', 'unit': 'descriptor-specific', 'invariance': 'same geometry and method', 'requires': 'both temperature reports', 'missing': 'one temperature report or value missing'},
    },
}

_STATUS = {
    'recommended': {
        'xtb_energy_ha', 'xtb_homo_ev', 'xtb_lumo_ev', 'xtb_gap_ev',
        'xtb_dip_debye', 'xtb_dip_*_ebohr', 'xtb_q_*_e', 'xtb_wbo_sum',
        'q_e', 'wbo_sum', 'wbo',
        'mul_*_e', 'low_*_e',
        'f_plus/f_minus/f_radical/dual', 's_loc_1ev',
        'omega_loc_ev',
        'ip_v_ev', 'ea_v_ev', 'mu_ev',
        'chi_ev', 'eta_ev', 's_1ev', 'omega_ev',
    },
    'diagnostic': {
        'xtb_n_*', 'xtb_*_idx', 'xtb_occ_cut',
        'conformer_dispersion', 'gfn12_delta',
        'et300_1000_delta',
    },
    'exploratory': {
        'mbo', 'mul_ovlp_e', 'pi_n_e',
        'pi_q_e', 'pi_bo', 'n_loc_1ev',
    },
}


def DescriptorCatalog():
    """Return a copy of the canonical descriptor definitions."""
    return deepcopy(_CATALOG)


def DescriptorStatus():
    """Return the intended validation maturity of descriptor families."""
    return {name: sorted(values) for name, values in _STATUS.items()}

