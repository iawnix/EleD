"""Run after installing eled in an environment containing RDKit and xTB."""
import json
from rdkit import Chem
from rdkit.Chem import AllChem
from eled import (
    XTBCache, XTBOptions, RunXTBBatch, AnnotateMol, CalcXTBDescriptors,
)

mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
if AllChem.EmbedMolecule(mol, randomSeed=1) != 0:
    raise RuntimeError('could not generate the ethanol conformer')
cache = XTBCache()
results = RunXTBBatch([mol, mol], cache=cache, options=XTBOptions(backend='native'))
annotated = AnnotateMol(mol, results[0])
print(json.dumps({
    'descriptors': CalcXTBDescriptors(results[0]),
    'engine': results[0].metadata['engine'],
    'carbon_charge': annotated.GetAtomWithIdx(0).GetDoubleProp('xtb_partial_charge'),
    'cache_entries': len(cache),
    'input_unchanged': not mol.GetAtomWithIdx(0).HasProp('xtb_partial_charge'),
}, indent=2))
