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
        'xtb_q_*_e': {'definition': 'summary statistics of engine atomic charges; exact model is recorded as charge_model in provenance (GFN1/2 SCC Mulliken, GFN0 PEEQ/EEQ)', 'unit': 'e', 'invariance': 'atom reordering invariant for aggregate statistics', 'requires': 'partial charges and charge_model provenance', 'missing': 'charges not returned'},
        'xtb_cm5_q_*_e': {'definition': 'summary statistics of CM5 charges derived from GFN1/2 SCC Mulliken charges with the xTB CM5 geometry correction', 'unit': 'e', 'invariance': 'atom reordering invariant for aggregate statistics', 'requires': 'GFN1/2 SCC charges and calculation geometry', 'missing': 'CM5 is undefined for GFN0/EEQ or missing geometry'},
        'xtb_wbo_sum': {'definition': 'sum of all off-diagonal Wiberg pairs', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'xtb_n_*': {'definition': 'reported orbital/electron/unpaired-electron count', 'unit': 'count', 'invariance': 'molecule invariant', 'requires': 'engine result', 'missing': 'engine did not return count'},
        'xtb_*_idx': {'definition': 'zero-based orbital index under the recorded occupation convention', 'unit': 'index', 'invariance': 'depends on the retained orbital ordering', 'requires': 'orbital energies and occupations', 'missing': 'frontier unavailable'},
        'xtb_occ_cut': {'definition': 'occupation threshold used to define occupied/virtual orbitals', 'unit': 'dimensionless', 'invariance': 'protocol setting', 'requires': 'descriptor extraction', 'missing': 'never missing for a valid extraction'},
    },
    'atom': {
        'q_e': {'definition': 'engine atomic charge; exact model is result.charge_model (GFN1/2 SCC Mulliken or GFN0 PEEQ/EEQ)', 'unit': 'e', 'invariance': 'rigid-motion invariant; atom index follows input order', 'requires': 'partial charges and charge_model provenance', 'missing': 'charges not returned'},
        'cm5_q_e': {'definition': 'CM5 charge derived from an xTB GFN1/2 SCC Mulliken charge', 'unit': 'e', 'invariance': 'rigid-motion invariant; atom index follows input order', 'requires': 'GFN1/2 SCC charges and calculation geometry', 'missing': 'CM5 is undefined for GFN0/EEQ or missing geometry'},
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
        'omega_minus_loc_ev / omega_plus_loc_ev / delta_omega_loc_ev': {'definition': 'local electrodonating power (omega_minus * f_minus), local electroaccepting power (omega_plus * f_plus), and multiphilic difference (omega_plus_loc - omega_minus_loc)', 'unit': 'eV', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and condensed Fukui indices', 'missing': 'hardness nonpositive or charge state failed'},
        's2_loc_1ev2': {'definition': 'local hyper-softness S^2 * dual = s_1ev^2 * (f_plus - f_minus)', 'unit': '1/eV^2', 'invariance': 'same as Fukui index', 'requires': 'positive hardness and dual descriptor', 'missing': 'hardness nonpositive or charge state failed'},
        'fod_pop_e / fod_s_e / fod_p_e / fod_d_e': {'definition': 'Grimme atomic and s/p/d shell fractional occupation number weighted electron density (FOD) static-correlation populations', 'unit': 'e', 'invariance': 'rigid-motion, AO permutation, MO phase and degenerate subspace rotation invariant; preserves RDKit atom order', 'requires': 'retained overlap basis (keepOverlap=True) and finite electronic temperature', 'missing': 'atoms empty with overlap_basis_not_retained when keepOverlap=False'},
        'chi_diag_e_per_ev / atomic_softness_e_per_ev / nonlocal_coupling_e_per_ev': {'definition': 'diagonal independent-particle atom-condensed linear response chi_AA (<= 0), atomic softness -chi_AA (>= 0), and nonlocal coupling sum_{B!=A} |chi_AB|', 'unit': 'e/eV (also chi_diag_e_per_ha in e/Hartree)', 'invariance': 'rigid-motion, AO permutation, MO phase and degenerate subspace rotation invariant; preserves RDKit atom order', 'requires': 'closed-shell integer-occupied retained basis (keepOverlap=True) with positive HOMO-LUMO gap', 'missing': 'withheld for open-shell, unknown spin, fractional occupation or degenerate gap'},
        'f_w_plus/f_w_minus/f_w_radical/dual_w': {'definition': 'Pino-Rios / Martinez-Araya orbital-weighted condensed Fukui indices from a single retained calculation', 'unit': 'dimensionless', 'invariance': 'AO permutation, MO phase, degenerate subspace rotation and rigid-motion invariant', 'requires': 'retained overlap, MO coefficients, AO-atom map, and frontier orbitals', 'missing': 'HOMO or LUMO unavailable'},
        'occ_n_e / vir_n': {'definition': 'sum of occupied / virtual atomic orbital populations (Lowdin or Mulliken)', 'unit': 'e (occ_n_e), dimensionless (vir_n)', 'invariance': 'AO permutation, MO phase and rigid-motion invariant', 'requires': 'retained overlap, MO coefficients and AO-atom map', 'missing': 'occupied or virtual subspace unavailable'},
        'occ_e_mean_ev / vir_e_mean_ev': {'definition': 'population-weighted mean occupied (negative condensed ALIE proxy: ALIE_A = -occ_e_mean_ev) and virtual atomic orbital energies', 'unit': 'eV', 'invariance': 'AO permutation, MO phase, degenerate subspace rotation and rigid-motion invariant', 'requires': 'retained basis and nonzero occupied/virtual atom weight', 'missing': 'subspace absent or zero atomic population weight'},
        'occ_e_frontier_ev / vir_e_frontier_ev': {'definition': 'Gaussian frontier-weighted occupied and virtual atomic orbital energies', 'unit': 'eV', 'invariance': 'same as occ_e_mean_ev / vir_e_mean_ev', 'requires': 'retained basis, frontier orbital and nonzero weighted denominator', 'missing': 'frontier orbital absent or zero weighted denominator'},
        'pdos_homo_win_e / pdos_lumo_win': {'definition': 'projected density of states integrated in [HOMO-window, HOMO] (e) and [LUMO, LUMO+window] (dimensionless)', 'unit': 'e (pdos_homo_win_e), dimensionless (pdos_lumo_win)', 'invariance': 'AO permutation, MO phase and rigid-motion invariant', 'requires': 'retained basis and HOMO/LUMO frontier', 'missing': 'HOMO or LUMO unavailable'},
        'epn_au / epn_kcalmol': {'definition': 'exact valence electronic plus other-nuclei electrostatic potential at each nuclear position (EPN, V_N), omitting the singular self-nuclear term Z_A/r_iA', 'unit': 'Hartree/e (epn_au), kcal/(mol*e) (epn_kcalmol)', 'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices', 'requires': 'libxtb extended properties ABI v1', 'missing': 'properties_extension_unavailable'},
        'alpha_bohr3 / alpha_angstrom3': {'definition': 'D4 atomic static dipole polarizability', 'unit': 'Bohr^3 (alpha_bohr3), Angstrom^3 (alpha_angstrom3)', 'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices', 'requires': 'GFN2-xTB or GFN0-xTB with libxtb extended properties ABI v1', 'missing': 'not_computed_by_gfn1_d3_dispersion or properties_extension_unavailable'},
        'c6_aa_ha_bohr6': {'definition': 'D4 homoatomic dispersion coefficient C6_AA', 'unit': 'Hartree*Bohr^6', 'invariance': 'rigid-motion invariant; atom-permutation invariant with remapped indices', 'requires': 'GFN2-xTB or GFN0-xTB with libxtb extended properties ABI v1', 'missing': 'not_computed_by_gfn1_d3_dispersion or properties_extension_unavailable'},
        'dip_ebohr / dip_norm_ebohr / dip_norm_debye': {'definition': 'CAMM atomic dipole vector and rotational-invariant magnitude', 'unit': 'e*Bohr / Debye', 'invariance': 'norm is rigid-motion invariant; Cartesian vector rotates with the coordinate frame', 'requires': 'GFN2-xTB (SCC CAMM) or GFN1-xTB (post-SCF CAMM) with properties ABI v1', 'missing': 'not_computed_by_gfn0 or properties_extension_unavailable'},
        'quad_ebohr2 / quad_frobenius_ebohr2 / quad_anisotropy_ebohr2': {'definition': 'traceless CAMM atomic quadrupole tensor (xx,xy,yy,xz,yz,zz), Frobenius norm and rotational anisotropy', 'unit': 'e*Bohr^2', 'invariance': 'eigenvalues, Frobenius norm and anisotropy are rigid-motion invariant; tensor components rotate with the coordinate frame', 'requires': 'GFN2-xTB or GFN1-xTB with properties ABI v1', 'missing': 'not_computed_by_gfn0 or properties_extension_unavailable'},
    },
    'fragment': {
        'q_e': {'definition': 'sum of engine atomic charges in the named fragment; exact model is recorded in provenance', 'unit': 'e', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'partial charges, charge_model provenance and explicit fragment indices', 'missing': 'charges not returned'},
        'cm5_q_e': {'definition': 'sum of CM5 charges in the named fragment', 'unit': 'e', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'GFN1/2 CM5 charges and explicit fragment indices', 'missing': 'CM5 unavailable'},
        'q_min_e / q_max_e / q_sd_e': {'definition': 'minimum, maximum and population standard deviation of fragment atomic charges', 'unit': 'e', 'invariance': 'atom reordering invariant within the fragment', 'requires': 'partial charges and explicit fragment indices', 'missing': 'charges not returned'},
        'wbo_intra / wbo_inter': {'definition': 'within-fragment and fragment-to-outside Wiberg sums', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix and explicit fragment indices', 'missing': 'bond orders not returned'},
    },
    'bond': {
        'wbo': {'definition': 'xTB Wiberg bond order on an RDKit bond', 'unit': 'dimensionless', 'invariance': 'atom reordering invariant with remapped indices', 'requires': 'Wiberg matrix', 'missing': 'bond orders not returned'},
        'mul_ovlp_e': {'definition': 'signed Mulliken bond overlap of the one-electron-normalized HOMO/LUMO group projector', 'unit': 'e', 'invariance': 'frontier-group rotations and orbital phase invariant', 'requires': 'frontier group and same-calculation P/S', 'missing': 'no frontier group or basis missing'},
        'pi_bo': {'definition': 'closed-shell projected Mayer (Mulliken scheme) or squared orthogonal-density (Lowdin scheme) pi bond index', 'unit': 'dimensionless', 'invariance': 'rigid-motion and normal-sign invariant', 'requires': 'closed shell and two valid local p projections', 'missing': 'open/unknown spin state or missing endpoint projection'},
        'coop_e': {'definition': 'crystal orbital overlap population integrated over occupied MOs on an RDKit bond (positive bonding, negative antibonding)', 'unit': 'e', 'invariance': 'AO permutation, MO phase, degenerate subspace rotation, bond reversal and rigid-motion invariant', 'requires': 'retained overlap, MO coefficients, AO-atom map and occupied orbitals', 'missing': 'no occupied orbitals under threshold or basis missing'},
        'ecohp_ev': {'definition': 'eigenvalue-weighted bond overlap stabilization proxy Tr(W_AB S_BA) (negative bonding, positive antibonding)', 'unit': 'eV', 'invariance': 'AO permutation, MO phase, degenerate subspace rotation, bond reversal and rigid-motion invariant', 'requires': 'retained overlap, MO coefficients, AO-atom map and occupied orbitals', 'missing': 'no occupied orbitals under threshold or basis missing'},
        'delta_q_e / signed_delta_q_e / q_prod_e2': {'definition': 'absolute bond charge difference, signed difference q(begin)-q(end), and charge product from engine atomic charges', 'unit': 'e / e / e^2', 'invariance': 'delta_q_e and q_prod_e2 are bond-reversal and rigid-motion invariant; signed_delta_q_e antisymmetrizes on bond reversal', 'requires': 'atomic charges and RDKit bond list', 'missing': 'charges not returned'},
        'delta_cm5_q_e / signed_delta_cm5_q_e / cm5_q_prod_e2': {'definition': 'bond polarity descriptors computed from CM5 charges', 'unit': 'e / e / e^2', 'invariance': 'same as engine-charge bond polarity descriptors', 'requires': 'GFN1/2 CM5 charges', 'missing': 'cm5_unavailable_for_gfn0 or missing CM5 charges'},
        'e_parallel_au / e_parallel_mv_cm / e_perp_au / e_norm_au / e_vec_au': {'definition': 'bond-midpoint local electric field E = -grad V from 6-point central finite differences of libxtb ESP, projected onto the RDKit bond vector (1 a.u. = 5142.2067476326 MV/cm)', 'unit': 'Hartree/(e*Bohr) and MV/cm', 'invariance': 'e_norm and e_perp are rigid-motion and bond-reversal invariant; e_parallel antisymmetrizes on bond reversal; preserves RDKit bond order', 'requires': 'ESP extension ABI v1 and 3D conformer with bonds', 'missing': 'bonds empty with molecule_has_no_bonds when molecule has 0 bonds'},
    },
    'density': {
        'density_matrix': {'definition': 'spin-summed AO density P = C diag(n) C.T', 'unit': 'basis-dependent; trace(P S) is e', 'invariance': 'matrix transforms with AO basis; trace with S and MO phase are invariant', 'requires': 'same-calculation C, S and occupations', 'missing': 'overlap/basis extension unavailable'},
        'tr_ps_e': {'definition': 'electron count from trace(P S)', 'unit': 'e', 'invariance': 'AO permutation invariant with transformed S', 'requires': 'density matrix and overlap', 'missing': 'basis data unavailable'},
        'n_e': {'definition': 'reported modeled electron count; occupation sum if absent', 'unit': 'e', 'invariance': 'AO and MO representation invariant', 'requires': 'retained density result', 'missing': 'basis data unavailable'},
        'n_err_e': {'definition': 'trace(P S) minus reported modeled electron count (occupation sum if absent)', 'unit': 'e', 'invariance': 'AO permutation invariant', 'requires': 'density matrix and occupations', 'missing': 'basis data unavailable'},
        'shell_err_e': {'definition': 'maximum shell-to-atom population sum discrepancy', 'unit': 'e', 'invariance': 'AO permutation invariant with maps', 'requires': 'shell map and density', 'missing': 'shell map unavailable'},
        'mci': {'definition': 'closed-shell Giambiagi-Bultinck multicenter index (and i_ring) computed by AO-block matrix multiplication and trace over fixed-first cyclic permutations of an ordered ring', 'unit': 'dimensionless', 'invariance': 'rigid-motion, AO orthogonal/permutation/sign, MO phase and ring start/direction invariant', 'requires': 'same-calculation closed-shell density, overlap and ordered bonded ring (3..9 atoms)', 'missing': 'open-shell/unknown spin state, basis data unavailable, or ring invalid'},
    },
    'reactivity': {
        'ip_v_ev': {'definition': 'vertical ionization energy from E(Q+1)-E(Q)', 'unit': 'eV', 'invariance': 'rigid-motion invariant at fixed geometry', 'requires': 'successful Q and Q+1 GFN1/2 calculations', 'missing': 'required charge state failed'},
        'ea_v_ev': {'definition': 'vertical electron affinity from E(Q)-E(Q-1)', 'unit': 'eV', 'invariance': 'rigid-motion invariant at fixed geometry', 'requires': 'successful Q and Q-1 GFN1/2 calculations', 'missing': 'required charge state failed'},
        'mu_ev': {'definition': 'negative electronegativity, -(IP+EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        'chi_ev': {'definition': 'electronegativity, (IP+EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        'eta_ev': {'definition': 'conceptual-DFT hardness, (IP-EA)/2', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'vertical IP and EA', 'missing': 'required charge state failed'},
        's_1ev': {'definition': 'global softness, 1/(IP-EA)', 'unit': '1/eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'omega_ev': {'definition': 'global electrophilicity, mu^2/(2 hardness)', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'omega_minus_ev / omega_plus_ev': {'definition': 'Gazquez electrodonating power (3*IP+EA)^2/(16*(IP-EA)) and electroaccepting power (IP+3*EA)^2/(16*(IP-EA))', 'unit': 'eV', 'invariance': 'same as vertical IP/EA', 'requires': 'positive hardness', 'missing': 'hardness nonpositive or unavailable'},
        'f_plus/f_minus/f_radical/dual': {'definition': 'condensed finite-difference Fukui indices per atom', 'unit': 'dimensionless', 'invariance': 'fixed atom indexing and geometry; rigid-motion invariant', 'requires': 'compatible Q, Q+1 and Q-1 results', 'missing': 'required charge state failed'},
        'chi_e_per_ev / chi_e_per_ha': {'definition': 'independent-particle atom-condensed linear response matrix chi_AB = -4 sum_{i in occ, a in vir} Q^A_{ia} Q^B_{ia} / (eps_a - eps_i) with exact zero row/column sums', 'unit': 'e/eV (chi_e_per_ev), e/Hartree (chi_e_per_ha)', 'invariance': 'rigid-motion, AO permutation, MO phase and degenerate subspace rotation invariant', 'requires': 'closed-shell integer-occupied retained basis (keepOverlap=True) with positive HOMO-LUMO gap', 'missing': 'withheld with explicit missing reason for open-shell, unknown spin, fractional occupation or degenerate gap'},
    },
    'fod': {
        'n_fod_e': {'definition': 'Grimme total fractional occupation number weighted electron density static-correlation index N_FOD = sum_i f_FOD,i (matching xtb --fod)', 'unit': 'e', 'invariance': 'rigid-motion and orbital-phase invariant', 'requires': 'orbital energies, electron/spin counts and finite electronic temperature (default 5000 K in RunFOD)', 'missing': 'never missing on valid completion; 0.0 at T = 0 K'},
        'fermi_alpha_ev / fermi_beta_ev': {'definition': 'spin-resolved Fermi-Dirac chemical potentials solved at the FOD electronic temperature', 'unit': 'eV', 'invariance': 'rigid-motion and orbital-phase invariant', 'requires': 'orbital energies and spin-channel electron counts', 'missing': 'never missing on valid completion'},
        'orbital_fod_e': {'definition': 'per-MO spin-summed Grimme FOD weights (1 - f for eps <= E_F, f for eps > E_F)', 'unit': 'e', 'invariance': 'rigid-motion and orbital-phase invariant', 'requires': 'orbital energies and finite electronic temperature', 'missing': 'never missing on valid completion'},
    },
    'realspace': {
        'density_au / vir_density_au': {'definition': 'exact Cartesian GTO valence electron density rho(r) and virtual orbital density evaluated via libxtb realspace ABI v1', 'unit': 'e/Bohr^3 (density_au), 1/Bohr^3 (vir_density_au)', 'invariance': 'rigid-motion invariant when points co-transform, MO phase and degenerate subspace invariant', 'requires': 'libxtb realspace extension ABI v1 and realSpacePointsAngstrom', 'missing': 'DescriptorUnavailableError if realSpacePointsAngstrom was not supplied'},
        'alie_ev / lea_ev': {'definition': 'exact real-space Average Local Ionization Energy ALIE(r) and Local Electron Affinity LEA(r) at sampled 3D points', 'unit': 'eV', 'invariance': 'same as density_au', 'requires': 'libxtb realspace extension ABI v1 and density > densityCutoffAU', 'missing': 'None when density (or virtual density) <= densityCutoffAU; never zero-filled'},
        'alie_min_ev / alie_max_ev / alie_mean_ev / alie_sigma2_ev2 / alie_pi_ev': {'definition': 'area-weighted molecular and per-atom surface ALIE statistics on a vdW or 0.001 a.u. density isosurface grid', 'unit': 'eV (and eV^2 for sigma2)', 'invariance': 'rigid-motion and point-permutation invariant when surface grid co-transforms; preserves RDKit atom order', 'requires': 'sampled real-space properties on a surface grid', 'missing': 'None with density_below_cutoff or no_surface_points_assigned'},
        'lea_min_ev / lea_max_ev / lea_mean_ev / lea_sigma2_ev2 / lea_pi_ev': {'definition': 'area-weighted molecular and per-atom surface LEA statistics on a vdW or 0.001 a.u. density isosurface grid', 'unit': 'eV (and eV^2 for sigma2)', 'invariance': 'same as surface ALIE statistics', 'requires': 'sampled real-space properties with virtual orbitals and virtual density > densityCutoffAU', 'missing': 'None with virtual_density_below_cutoff or no_surface_points_assigned'},
        'volume_a3 / volume_bohr3 / area_a2 / area_bohr2 / sphericity': {'definition': 'molecular volume, surface area and sphericity of the exact electron-density isosurface (default 0.001 e/Bohr^3) computed via Voronoi-partitioned batched radial bisection in one xTB SCF', 'unit': 'Angstrom^3 / Bohr^3 / Angstrom^2 / Bohr^2 / dimensionless', 'invariance': 'translation invariant and rotation invariant up to finite Fibonacci angular mesh discretization', 'requires': 'libxtb realspace extension ABI v1 (and ESP ABI v1 when includeESP=True)', 'missing': 'raises ValueError if isovalue cannot be bracketed'},
    },
    'local_geometry': {
        'l_a / b1_a / b5_a': {'definition': 'Sterimol axial extent plus recorded correction and minimum/maximum transverse sphere-envelope support', 'unit': 'angstrom', 'invariance': 'rigid-motion invariant up to the reported B1 angular error bound', 'requires': '3D conformer, directed bonded axis, selected atoms and radii', 'missing': 'invalid axis or radii raise ValueError'},
        'vbur_percent': {'definition': 'fraction of centered sphere occupied by the union of selected atom spheres, times 100', 'unit': 'percent', 'invariance': 'translation invariant; rotation dependence decreases with Cartesian grid spacing', 'requires': 'center atom excluded from selected atoms, sphere radius and atom radii', 'missing': 'invalid selection or grid raises ValueError'},
        'plane_rms_a / plane_max_a': {'definition': 'RMS and maximum absolute distance to an unweighted least-squares plane', 'unit': 'angstrom', 'invariance': 'rigid-motion invariant', 'requires': 'at least three selected atoms', 'missing': 'None when the plane normal is not unique'},
        'plane_angle_deg': {'definition': 'acute angle between two least-squares planes', 'unit': 'degree', 'invariance': 'rigid-motion and normal-sign invariant', 'requires': 'two atom selections with unique plane normals', 'missing': 'None when either plane normal is not unique'},
        'centroid_distance_a': {'definition': 'distance between unweighted atom centroids', 'unit': 'angstrom', 'invariance': 'rigid-motion invariant', 'requires': 'two selected atom sets', 'missing': 'second selection not requested'},
        'homa': {'definition': '1 - mean(alpha_AB * (R_AB - Ropt_AB)^2) on the ordered ring', 'unit': 'dimensionless', 'invariance': 'rigid-motion and cyclic/reverse ring-order invariant', 'requires': 'bonded ring path and element-pair parameters', 'missing': 'None if any element pair lacks parameters'},
        'q_total_a / q_alternating_a / mode_q_a': {'definition': 'Cremer-Pople total, signed alternating and paired Fourier-mode amplitudes', 'unit': 'angstrom', 'invariance': 'proper rigid-motion invariant; alternating sign depends on ring order', 'requires': 'ordered bonded ring with at least four atoms', 'missing': 'None when the Cremer-Pople ring normal is undefined'},
        'theta_deg / mode_phi_deg': {'definition': 'Cremer-Pople six-membered polar angle and mode phases', 'unit': 'degree', 'invariance': 'proper rigid-motion invariant; phases depend on ring origin and direction', 'requires': 'defined ring normal and nonzero relevant amplitude', 'missing': 'None for undefined normal or zero relevant amplitude'},
    },
    'esp_surface': {
        'v_min_au / v_max_au / v_mean_au': {'definition': 'minimum, maximum and weighted mean molecular electrostatic potential on an explicit scaled RDKit vdW Fibonacci shell (not a 0.001 a.u. density isosurface)', 'unit': 'Hartree/e (optional *_kcalmol in kcal/(mol*e))', 'invariance': 'rigid-motion, point-permutation and atom-renumbering invariant when grid co-transforms', 'requires': 'sampled ESP values on an explicit surface grid', 'missing': 'all surface weights zero raises ValueError'},
        'v_pos_mean_au / v_neg_mean_au': {'definition': 'weighted mean positive (V > 0) and negative (V < 0) surface electrostatic potentials', 'unit': 'Hartree/e (optional *_kcalmol)', 'invariance': 'same as v_mean_au', 'requires': 'positive or negative surface points', 'missing': 'None with no_positive_surface_points / no_negative_surface_points when single-sign'},
        'v_pos_frac / v_neg_frac': {'definition': 'weighted fraction of active surface with V > 0 and V < 0', 'unit': 'dimensionless', 'invariance': 'same as v_mean_au', 'requires': 'sampled surface ESP', 'missing': 'never missing on valid active grid'},
        'pi_esp_au': {'definition': 'Politzer-Murray average local polarity (weighted mean absolute deviation of V from v_mean)', 'unit': 'Hartree/e (optional pi_esp_kcalmol)', 'invariance': 'same as v_mean_au', 'requires': 'sampled surface ESP', 'missing': 'never missing on valid active grid'},
        'sigma2_pos_au2 / sigma2_neg_au2 / sigma2_tot_au2': {'definition': 'Politzer-Murray weighted population variances (ddof=0) of positive, negative and combined surface ESP', 'unit': '(Hartree/e)^2 (optional *_kcalmol2)', 'invariance': 'same as v_mean_au', 'requires': 'sampled surface ESP', 'missing': 'absent-sign variance is None when surface is single-sign'},
        'nu / nu_sigma2_au2': {'definition': 'Politzer-Murray electrostatic balance parameter nu = (sigma2_pos * sigma2_neg) / sigma2_tot^2 (in [0, 0.25]) and product nu * sigma2_tot', 'unit': 'dimensionless (nu), (Hartree/e)^2 (nu_sigma2_au2)', 'invariance': 'same as v_mean_au', 'requires': 'both positive and negative surface points and positive sigma2_tot', 'missing': 'None when single-sign or sigma2_tot == 0'},
        'sigma2_au2': {'definition': 'per-atom weighted population variance (ddof=0) of surface ESP around atom v_mean', 'unit': '(Hartree/e)^2 (optional sigma2_kcalmol2)', 'invariance': 'same as v_mean_au; preserves RDKit atom ordering', 'requires': 'sampled surface ESP with assigned atom points', 'missing': 'None with no_surface_points_assigned or zero_surface_weight_assigned for buried atoms'},
    },
    'localized_orbitals': {
        'localization_functional / initial_functional': {'definition': 'Pipek-Mezey localization functional sum_i sum_A (Q^A_{ii})^2 after and before 2x2 Jacobi unitary rotation', 'unit': 'dimensionless', 'invariance': 'AO permutation, MO phase, degenerate occupied rotation and rigid-motion invariant', 'requires': 'closed-shell integer-occupied retained basis (keepOverlap=True)', 'missing': 'None when withheld for open-shell, unknown spin, fractional occupation or empty occupied subspace'},
        'lmo_energy_ev / n_centers / kind': {'definition': 'rotated Fock diagonal energy U^T diag(eps) U, inverse participation center count 1/sum_A q_A^2, and lone_pair/bond/delocalized classification per LMO', 'unit': 'eV / dimensionless / category', 'invariance': 'same as localization_functional', 'requires': 'converged or completed Pipek-Mezey localization', 'missing': 'empty orbitals list when withheld'},
        'n_lone_pairs / n_localized_bonds / max_lone_pair_energy_ev': {'definition': 'per-atom Pipek-Mezey lone-pair count, localized-bond count and highest lone-pair energy', 'unit': 'count / count / eV', 'invariance': 'rigid-motion invariant; preserves RDKit atom ordering', 'requires': 'completed Pipek-Mezey localization', 'missing': 'max_lone_pair_energy_ev is None with no_lone_pairs_on_atom when n_lone_pairs == 0'},
    },
    'properties': {
        'e_elec_ha / e_rep_ha / e_es_ha / e_aes_ha / e_axc_ha / e_disp_ha / e_xb_ha': {'definition': 'xTB energy decomposition components from libxtb extended properties ABI v1', 'unit': 'Hartree', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'libxtb extended properties ABI v1', 'missing': 'None with method-specific reason (e.g. e_aes/e_axc unavailable for GFN0/1; e_xb unavailable for GFN2)'},
        'g_solv_ha / g_born_ha / g_sasa_ha / g_hb_ha / g_shift_ha': {'definition': 'implicit solvation energy and Born/SASA/HB/shift sub-components reported by libxtb', 'unit': 'Hartree', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'implicit solvation (solvent != None) and properties ABI v1', 'missing': 'None with gas_phase_no_solvation in gas phase or not_separately_returned_by_gfn0 for GFN0 sub-terms'},
        'fermi_alpha_ha / fermi_beta_ha / fermi_alpha_ev / fermi_beta_ev': {'definition': 'electronic Fermi levels from Fermi-Dirac smearing', 'unit': 'Hartree / eV', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'properties ABI v1 and electronicTemperature > 0.1 K', 'missing': 'None with fermi_smearing_disabled_low_etemp when electronicTemperature <= 0.1 K'},
        'mol_pol_bohr3 / mol_pol_angstrom3 / mol_c6_ha_bohr6 / mol_c8_ha_bohr8': {'definition': 'analytical D4 molecular isotropic static polarizability and dispersion C6/C8 coefficients', 'unit': 'Bohr^3 / Angstrom^3 / Hartree*Bohr^6 / Hartree*Bohr^8', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'GFN2-xTB or GFN0-xTB with properties ABI v1', 'missing': 'None with not_computed_by_gfn1_d3_dispersion for GFN1-xTB'},
        'mol_quadrupole_ebohr2 / mol_quadrupole_frobenius_ebohr2 / mol_quadrupole_anisotropy_ebohr2': {'definition': 'traceless molecular quadrupole tensor (xx,xy,yy,xz,yz,zz) and its Frobenius norm and anisotropy', 'unit': 'e*Bohr^2', 'invariance': 'eigenvalues, Frobenius norm and anisotropy are rotation invariant; origin dependent when net charge or dipole is nonzero', 'requires': 'GFN2-xTB or GFN1-xTB with properties ABI v1', 'missing': 'None with not_computed_by_gfn0 for GFN0-xTB'},
    },
    'solvation': {
        'delta_g_solv_ha / delta_g_solv_kcalmol': {'definition': 'vertical implicit solvation free energy E(solvent) - E(gas) at fixed geometry via ALPB or GBSA', 'unit': 'Hartree / kcal/mol', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'matched gas-phase and solvated single points', 'missing': 'calculation failure raises XTBError'},
        'g_solv_direct_ha / solute_polarization_energy_ha': {'definition': 'direct libxtb solvation energy and solute electronic polarization contribution delta_g_solv - g_solv_direct', 'unit': 'Hartree (and *_kcalmol)', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'solvated single point with properties ABI v1', 'missing': 'None with properties_extension_unavailable'},
        'delta_dipole_debye / delta_homo_ev / delta_lumo_ev / delta_gap_ev': {'definition': 'solvent-induced shifts (solvated minus gas phase) in dipole magnitude and frontier orbital energies', 'unit': 'Debye / eV', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'matched gas-phase and solvated single points', 'missing': 'frontier shifts are None when HOMO or LUMO is unavailable'},
        'delta_q_e / delta_cm5_q_e': {'definition': 'per-atom charge redistribution q(solvent) - q(gas) for engine charges and CM5 charges', 'unit': 'e', 'invariance': 'rigid-motion invariant; preserves RDKit atom ordering', 'requires': 'matched gas-phase and solvated single points', 'missing': 'delta_cm5_q_e is None with cm5_undefined_for_gfn0_eeq for GFN0'},
    },
    'thermochemistry': {
        'frequencies_cm1 / ir_intensities_km_mol / rotational_constants_cm1': {'definition': 'harmonic vibrational frequencies (negative for imaginary modes), double-harmonic IR intensities from dipole derivatives, and principal rotational constants', 'unit': 'cm^-1 / km/mol / cm^-1', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'central finite-difference Hessian from analytical GetGradient()', 'missing': 'empty lists for monoatomic systems (0 vibrational modes)'},
        'zpe_ha / thermal_energy_correction_ha / thermal_enthalpy_correction_ha / entropy_j_mol_k / gibbs_correction_ha / total_gibbs_free_energy_ha': {'definition': 'ideal-gas RRHO or Grimme qRRHO zero-point energy, thermal corrections, molar entropy and total Gibbs free energy (with kcal/mol counterparts)', 'unit': 'Hartree / J/(mol*K) / kcal/mol', 'invariance': 'rigid-motion and atom-permutation invariant', 'requires': 'harmonic frequency analysis; non-stationary geometries are explicitly flagged via is_stationary_point and warnings', 'missing': 'never missing on valid completion; imaginary modes excluded from vibrational sums and flagged'},
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
        'omega_minus_ev / omega_plus_ev',
        'omega_minus_loc_ev / omega_plus_loc_ev / delta_omega_loc_ev',
        's2_loc_1ev2',
        'n_fod_e', 'fod_pop_e / fod_s_e / fod_p_e / fod_d_e',
        'density_au / vir_density_au', 'alie_ev / lea_ev',
        'alie_min_ev / alie_max_ev / alie_mean_ev / alie_sigma2_ev2 / alie_pi_ev',
        'lea_min_ev / lea_max_ev / lea_mean_ev / lea_sigma2_ev2 / lea_pi_ev',
        'volume_a3 / volume_bohr3 / area_a2 / area_bohr2 / sphericity',
        'e_parallel_au / e_parallel_mv_cm / e_perp_au / e_norm_au / e_vec_au',
        'v_min_au / v_max_au / v_mean_au', 'v_pos_mean_au / v_neg_mean_au',
        'v_pos_frac / v_neg_frac', 'pi_esp_au',
        'sigma2_pos_au2 / sigma2_neg_au2 / sigma2_tot_au2', 'nu / nu_sigma2_au2', 'sigma2_au2',
        'epn_au / epn_kcalmol',
        'mol_pol_bohr3 / mol_pol_angstrom3 / mol_c6_ha_bohr6 / mol_c8_ha_bohr8',
        'alpha_bohr3 / alpha_angstrom3', 'c6_aa_ha_bohr6',
        'e_elec_ha / e_rep_ha / e_es_ha / e_aes_ha / e_axc_ha / e_disp_ha / e_xb_ha',
        'g_solv_ha / g_born_ha / g_sasa_ha / g_hb_ha / g_shift_ha',
        'delta_g_solv_ha / delta_g_solv_kcalmol',
        'delta_dipole_debye / delta_homo_ev / delta_lumo_ev / delta_gap_ev',
        'delta_q_e / delta_cm5_q_e',
        'frequencies_cm1 / ir_intensities_km_mol / rotational_constants_cm1',
        'zpe_ha / thermal_energy_correction_ha / thermal_enthalpy_correction_ha / entropy_j_mol_k / gibbs_correction_ha / total_gibbs_free_energy_ha',
    },
    'diagnostic': {
        'xtb_n_*', 'xtb_*_idx', 'xtb_occ_cut',
        'fermi_alpha_ha / fermi_beta_ha / fermi_alpha_ev / fermi_beta_ev',
        'fermi_alpha_ev / fermi_beta_ev', 'orbital_fod_e',
        'g_solv_direct_ha / solute_polarization_energy_ha',
    },
    'exploratory': {
        'mul_ovlp_e', 'pi_n_e',
        'pi_q_e', 'pi_bo', 'mci', 'n_loc_1ev',
        'f_w_plus/f_w_minus/f_w_radical/dual_w',
        'occ_n_e / vir_n', 'occ_e_mean_ev / vir_e_mean_ev',
        'occ_e_frontier_ev / vir_e_frontier_ev',
        'pdos_homo_win_e / pdos_lumo_win',
        'coop_e', 'ecohp_ev',
        'delta_q_e / signed_delta_q_e / q_prod_e2',
        'delta_cm5_q_e / signed_delta_cm5_q_e / cm5_q_prod_e2',
        'localization_functional / initial_functional',
        'lmo_energy_ev / n_centers / kind',
        'n_lone_pairs / n_localized_bonds / max_lone_pair_energy_ev',
        'dip_ebohr / dip_norm_ebohr / dip_norm_debye',
        'quad_ebohr2 / quad_frobenius_ebohr2 / quad_anisotropy_ebohr2',
        'mol_quadrupole_ebohr2 / mol_quadrupole_frobenius_ebohr2 / mol_quadrupole_anisotropy_ebohr2',
        'chi_e_per_ev / chi_e_per_ha',
        'chi_diag_e_per_ev / atomic_softness_e_per_ev / nonlocal_coupling_e_per_ev',
    },
}


def _exact_catalog():
    """Expand grouped definitions into the concrete public field names."""
    catalog = deepcopy(_CATALOG)
    atom = catalog['atom']
    for prefix in ('mul', 'low'):
        grouped = atom[f'{prefix}_*_e']
        for channel in ('s', 'p', 'd'):
            atom[f'{prefix}_{channel}_e'] = deepcopy(grouped)
        atom[f'{prefix}_qrec_e'] = deepcopy(atom['mul_qrec_e / low_qrec_e'])
    density = catalog['density']
    density['density_matrix']['unit'] = 'basis-dependent matrix; trace with S is e'
    # Keep wildcard entries as compact documentation, but make exact lookups
    # deterministic for feature-table builders and schema validators.
    return catalog


def DescriptorCatalog():
    """Return a copy of the canonical descriptor definitions."""
    return _exact_catalog()


def DescriptorStatus():
    """Return the intended validation maturity of descriptor families."""
    return {name: sorted(values) for name, values in _STATUS.items()}

