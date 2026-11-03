"""Run after building the optional rdXTBTools extension."""
from rdkit import Chem
from rdkit.Chem import AllChem, rdEHTTools
from rdkit_eth import rdXTBTools

mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0

ok, xtb = rdXTBTools.RunMol(mol, method='GFN2-xTB')
if not ok:
    raise RuntimeError(xtb.errorMessage)
print('xTB total energy (Hartree):', xtb.totalEnergyHartree)
print('xTB orbital energies (eV):', xtb.GetOrbitalEnergies())
print('xTB partial charges (e):', xtb.GetAtomicCharges())
print('xTB Wiberg bond orders:', xtb.GetWibergBondOrders())

ok, eht = rdEHTTools.RunMol(mol)
if not ok:
    raise RuntimeError('EHT calculation failed')
print('EHT total energy (eV):', eht.totalEnergy)
print('EHT orbital energies (eV):', eht.GetOrbitalEnergies())
# These total energies have different model definitions and reference zeros.
