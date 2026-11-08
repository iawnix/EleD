"""Run basic native descriptors, optionally retaining native orbital basis data.

python examples/electronic_descriptors.py
python examples/electronic_descriptors.py --populations
"""
import argparse
import json

from rdkit import Chem
from rdkit.Chem import AllChem
from eled import (
    XTBCalculator, CalcXTBDescriptors, CalcXTBAtomDescriptors,
    CalcXTBBondDescriptors, CalcXTBFragmentDescriptors,
    CalcOrbitalPopulations, CalcFrontierDescriptors,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--populations', action='store_true', help='retain native overlap for orbital population analysis')
    args = parser.parse_args()
    mol = Chem.AddHs(Chem.MolFromSmiles('CCO'))
    assert AllChem.EmbedMolecule(mol, randomSeed=1) == 0
    # Include each carbon's attached H atoms explicitly in the alkyl fragment.
    alkyl = {a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6}
    alkyl.update(n.GetIdx() for i in list(alkyl) for n in mol.GetAtomWithIdx(i).GetNeighbors() if n.GetAtomicNum() == 1)
    fragments = {'alkyl': sorted(alkyl), 'hydroxyl': sorted(set(range(mol.GetNumAtoms())) - alkyl)}
    result = XTBCalculator().calculate(mol, keep_overlap=args.populations)
    report = {
        'scalars': CalcXTBDescriptors(result),
        'atoms': CalcXTBAtomDescriptors(result),
        'bonds': CalcXTBBondDescriptors(mol, result),
        'fragments': CalcXTBFragmentDescriptors(result, fragments),
        'provenance': result.metadata,
    }
    if args.populations:
        report['mulliken'] = CalcOrbitalPopulations(result, scheme='mulliken')
        report['frontier'] = CalcFrontierDescriptors(mol, result, fragments=fragments)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
