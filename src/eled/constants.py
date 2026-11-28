"""Physical conversion constants shared by all modules (CODATA 2018 except length).

``eled.fod`` deliberately keeps xtb's own Boltzmann/Hartree constants so its
Fermi smearing matches ``xtb --fod``; thermochemistry keeps its SI set.
"""

HARTREE_TO_EV = 27.211386245988
# Length: the Bohr radius used to convert RDKit coordinates for libxtb. It must be
# identical in the ctypes backend and rdXTBTools (XTBTools.cpp: angstromToBohr)
# so both interfaces feed the engine bit-identical geometries; this is the
# CODATA 2014 value both have always used. Outputs use the same constant.
ANGSTROM_TO_BOHR = 1.8897261254578281
BOHR_TO_ANGSTROM = 1.0 / ANGSTROM_TO_BOHR
BOHR3_TO_ANGSTROM3 = BOHR_TO_ANGSTROM ** 3
HARTREE_TO_KCAL_MOL = 627.5094740631
HARTREE2_TO_KCAL_MOL2 = HARTREE_TO_KCAL_MOL * HARTREE_TO_KCAL_MOL
AU_DIPOLE_TO_DEBYE = 2.541746473
AU_EFIELD_TO_MV_CM = 5142.2067476326  # 1 Eh/(e*a0) = 5.1422067476326e11 V/m
