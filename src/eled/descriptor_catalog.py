"""Canonical, machine-readable descriptor definitions.

The catalog is intentionally descriptive: it does not calculate values.  A
field's ``missing`` text is part of its public contract and means callers must
preserve ``None``/the reason instead of imputing zero.
"""
from __future__ import annotations

from copy import deepcopy

_CATALOG = {
    'molecular': {
        'xtb_total_energy_hartree': {'definition': 'xTB total electronic energy', 'unit': 'Hartree', 'invariance': 'rigid-motion invariant; method/reference dependent', 'requires': 'single point', 'missing': 'engine did not return energy'},
        'xtb_homo_ev': {'definition': 'highest orbital with occupation above threshold', 'unit': 'eV', 'invariance': 'orbital phase and rigid-motion invariant', 'requires': 'orbital energies and occupations', 'missing': 'no occupied orbital under threshold'},
        'xtb_lumo_ev': {'definition': 'lowest orbital at or below occupation threshold', 'unit': 'eV', 'invariance': 'orbital phase and rigid-motion invariant', 'requires': 'orbital energies and occupations', 'missing': 'no virtual orbital under threshold'},
        'xtb_gap_ev': {'definition': 'LUMO minus HOMO', 'unit': 'eV', 'invariance': 'same as frontier energies', 'requires': 'HOMO and LUMO', 'missing': 'HOMO or LUMO unavailable'},
        'xtb_dipole_norm_debye': {'definition': 'norm of the xTB dipole', 'unit': 'Debye', 'invariance': 'rigid rotation invariant; origin dependent for charged systems', 'requires': 'dipole vector', 'missing': 'dipole not returned'},
        'xtb_dipole_*_e_bohr': {'definition': 'Cartesian xTB dipole component', 'unit': 'e*Bohr', 'invariance': 'components rotate with the coordinate frame; origin dependent for charged systems', 'requires': 'dipole vector', 'missing': 'dipole not returned'},
        'frontier_homo_* / frontier_lumo_*': {'definition': 'statistics of the degeneracy-grouped HOMO or LUMO atom distribution', 'unit': 'field dependent (eV, count, dimensionless or Angstrom)', 'invariance': 'rigid-motion invariant except centroid components; invariant to phase and rotations inside the selected degenerate group', 'requires': 'retained overlap and orbital basis', 'missing': 'frontier group or basis unavailable'},
        'xtb_charge_*_e': {'definition': 'summary statistics of xTB atomic partial charges', 'unit': 'e', 'invariance': 'atom reordering invariant for aggregate statistics', 'requires': 'partial charges', 'missing': 'charges not returned'},
        'xtb_wiberg_pair_sum': {'definition': 'sum of all off-diagonal Wiberg pairs', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'xtb_num_*': {'definition': 'reported orbital/electron/unpaired-electron count', 'unit': 'count', 'invariance': 'molecule invariant', 'requires': 'engine result', 'missing': 'engine did not return count'},
        'xtb_*_index': {'definition': 'zero-based orbital index under the recorded occupation convention', 'unit': 'index', 'invariance': 'depends on the retained orbital ordering', 'requires': 'orbital energies and occupations', 'missing': 'frontier unavailable'},
        'xtb_occupation_threshold': {'definition': 'occupation threshold used to define occupied/virtual orbitals', 'unit': 'dimensionless', 'invariance': 'protocol setting', 'requires': 'descriptor extraction', 'missing': 'never missing for a valid extraction'},
    },
    'atom': {
        'partial_charge_e': {'definition': 'xTB engine atomic partial charge', 'unit': 'e', 'invariance': 'rigid-motion invariant; atom index follows input order', 'requires': 'partial charges', 'missing': 'charges not returned'},
        'wiberg_sum': {'definition': 'sum of off-diagonal Wiberg orders incident on atom', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'mulliken_*_electron_count_e': {'definition': 'Mulliken density population by angular momentum shell', 'unit': 'e', 'invariance': 'AO permutation and orbital phase invariant with maps', 'requires': 'same-calculation P, S and shell map', 'missing': 'overlap/basis extension unavailable'},
        'lowdin_*_electron_count_e': {'definition': 'Löwdin density population by angular momentum shell', 'unit': 'e', 'invariance': 'AO permutation and orbital phase invariant with maps', 'requires': 'same-calculation P, S and shell map', 'missing': 'overlap/basis extension unavailable'},
        '*_p_to_s_ratio': {'definition': 'p-shell population divided by s-shell population', 'unit': 'dimensionless', 'invariance': 'same as shell population', 'requires': 's and p shell populations; nonzero s population', 'missing': 's population unavailable or zero'},
        'pi_electron_count_e': {'definition': 'p-shell density projected onto local plane normal', 'unit': 'e', 'invariance': 'rigid-motion and normal-sign invariant', 'requires': 'sp2/aromatic atom, defined plane, p shell map', 'missing': 'atom not sp2/aromatic or normal/p map undefined'},
        'pi_charge_e': {'definition': 'reference valence electrons minus pi population', 'unit': 'e', 'invariance': 'same as pi population', 'requires': 'pi population and supported element', 'missing': 'unsupported reference valence or pi population'},
        'f_plus/f_minus/f_radical/dual': {'definition': 'condensed finite-difference Fukui indices', 'unit': 'dimensionless', 'invariance': 'fixed atom indexing and geometry; rigid-motion invariant', 'requires': 'compatible Q, Q+1 and Q-1 results', 'missing': 'required charge state failed'},
        'local_softness_per_ev': {'definition': 'global softness times f_plus', 'unit': '1/eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and f_plus', 'missing': 'hardness undefined or charge state failed'},
        'local_electrophilicity_ev': {'definition': 'global electrophilicity times f_plus', 'unit': 'eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and f_plus', 'missing': 'electrophilicity undefined or charge state failed'},
        'local_nucleophilicity_per_ev': {'definition': 'softness-scaled f_minus proxy; no absolute nucleophilicity reference assumed', 'unit': '1/eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and f_minus', 'missing': 'softness undefined or charge state failed'},
    },
    'bond': {
        'wiberg_order': {'definition': 'xTB Wiberg bond order on an RDKit bond', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'mayer_bond_order': {'definition': 'closed-shell Mayer order from P S products', 'unit': 'dimensionless', 'invariance': 'AO permutation and orbital phase invariant', 'requires': 'same-calculation P/S and closed-shell density', 'missing': 'open-shell spin density unavailable or basis missing'},
        'mulliken_overlap_population': {'definition': 'signed HOMO/LUMO group Mulliken overlap population', 'unit': 'e', 'invariance': 'frontier-group rotations and orbital phase invariant', 'requires': 'frontier group and same-calculation P/S', 'missing': 'no frontier group or basis missing'},
        'pi_bond_order': {'definition': 'projected p-shell overlap between bonded atoms', 'unit': 'dimensionless', 'invariance': 'rigid-motion and normal-sign invariant', 'requires': 'two valid local p-shell normals', 'missing': 'one endpoint lacks a valid pi projection'},
    },
    'density': {
        'density_matrix': {'definition': 'spin-summed AO density P = C diag(n) C.T', 'unit': 'basis-dependent; trace(P S) is e', 'invariance': 'AO permutation and orbital phase invariant with maps', 'requires': 'same-calculation C, S and occupations', 'missing': 'overlap/basis extension unavailable'},
        'trace_PS_e': {'definition': 'electron count from trace(P S)', 'unit': 'e', 'invariance': 'AO permutation invariant with transformed S', 'requires': 'density matrix and overlap', 'missing': 'basis data unavailable'},
        'electron_conservation_error_e': {'definition': 'trace(P S) minus reported occupation sum', 'unit': 'e', 'invariance': 'AO permutation invariant', 'requires': 'density matrix and occupations', 'missing': 'basis data unavailable'},
        'shell_atom_conservation_error_e': {'definition': 'maximum shell-to-atom population sum discrepancy', 'unit': 'e', 'invariance': 'AO permutation invariant with maps', 'requires': 'shell map and density', 'missing': 'shell map unavailable'},
    },
    'reactivity': {
        'vertical_ip_ev': {'definition': 'vertical ionization energy from E(Q+1)-E(Q)', 'unit': 'eV', 'invariance': 'rigid-motion invariant at fixed geometry', 'requires': 'successful Q and Q+1 GFN1/2 calculations', 'missing': 'required charge state failed'},
        'vertical_ea_ev': {'definition': 'vertical electron affinity from E(Q)-E(Q-1)', 'unit': 'eV', 'invariance': 'rigid-motion invariant at fixed geometry', 'requires': 'successful Q and Q-1 GFN1/2 calculations', 'missing': 'required charge state failed'},
        'chemical_potential_ev': {'definition': 'negative electronegativity, -(IP+EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        'hardness_ev': {'definition': 'conceptual-DFT hardness, (IP-EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        'softness_per_ev': {'definition': 'global softness, 1/(IP-EA)', 'unit': '1/eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'electrophilicity_ev': {'definition': 'global electrophilicity, mu^2/(2 hardness)', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'f_plus/f_minus/f_radical/dual': {'definition': 'condensed finite-difference Fukui indices per atom', 'unit': 'dimensionless', 'invariance': 'fixed atom indexing and geometry; rigid-motion invariant', 'requires': 'compatible Q, Q+1 and Q-1 results', 'missing': 'required charge state failed'},
    },
    'diagnostic': {
        'conformer_dispersion': {'definition': 'population standard deviation across supplied conformers', 'unit': 'descriptor-specific', 'invariance': 'independent of conformer ordering', 'requires': 'at least two valid conformer reports', 'missing': 'insufficient valid conformers'},
        'gfn1_vs_gfn2_abs_difference': {'definition': 'absolute difference between method reports', 'unit': 'descriptor-specific', 'invariance': 'same atom/bond indexing and geometry protocol', 'requires': 'both method reports', 'missing': 'one method report or value missing'},
        'electronic_temperature_300_vs_1000_abs_difference': {'definition': 'absolute difference under electronic temperature perturbation', 'unit': 'descriptor-specific', 'invariance': 'same geometry and method', 'requires': 'both temperature reports', 'missing': 'one temperature report or value missing'},
    },
}

_STATUS = {
    'recommended': {
        'xtb_total_energy_hartree', 'xtb_homo_ev', 'xtb_lumo_ev', 'xtb_gap_ev',
        'xtb_dipole_norm_debye', 'xtb_dipole_*_e_bohr', 'xtb_charge_*_e', 'xtb_wiberg_pair_sum',
        'partial_charge_e', 'wiberg_sum', 'wiberg_order',
        'mulliken_*_electron_count_e', 'lowdin_*_electron_count_e',
        'f_plus/f_minus/f_radical/dual', 'local_softness_per_ev',
        'local_electrophilicity_ev',
        'vertical_ip_ev', 'vertical_ea_ev', 'chemical_potential_ev',
        'hardness_ev', 'softness_per_ev', 'electrophilicity_ev',
    },
    'diagnostic': {
        'xtb_num_*', 'xtb_*_index', 'xtb_occupation_threshold',
        'conformer_dispersion', 'gfn1_vs_gfn2_abs_difference',
        'electronic_temperature_300_vs_1000_abs_difference',
    },
    'exploratory': {
        'mayer_bond_order', 'mulliken_overlap_population', 'pi_electron_count_e',
        'pi_charge_e', 'pi_bond_order', 'local_nucleophilicity_per_ev',
    },
}


def DescriptorCatalog():
    """Return a copy of the canonical descriptor definitions."""
    return deepcopy(_CATALOG)


def DescriptorStatus():
    """Return the intended validation maturity of descriptor families."""
    return {name: sorted(values) for name, values in _STATUS.items()}
    return output

