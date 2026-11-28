"""Physical conversion constants shared by all modules (CODATA 2018).

``eled.fod`` deliberately keeps xtb's own Boltzmann/Hartree constants so its
Fermi smearing matches ``xtb --fod``; thermochemistry keeps its SI set.
"""

HARTREE_TO_EV = 27.211386245988
BOHR_TO_ANGSTROM = 0.529177210903
ANGSTROM_TO_BOHR = 1.0 / BOHR_TO_ANGSTROM
BOHR3_TO_ANGSTROM3 = BOHR_TO_ANGSTROM ** 3
HARTREE_TO_KCAL_MOL = 627.5094740631
HARTREE2_TO_KCAL_MOL2 = HARTREE_TO_KCAL_MOL * HARTREE_TO_KCAL_MOL
AU_DIPOLE_TO_DEBYE = 2.541746473
AU_EFIELD_TO_MV_CM = 5142.2067476326  # 1 Eh/(e*a0) = 5.1422067476326e11 V/m
