"""Sample water ESP and compare constrained ESP/single-stage RESP charges.

水分子静电势采样及带约束的 ESP / 单阶段 RESP 电荷示例。
Run using the installed EleD runtime: .eled/runtime/bin/python examples/esp_charges.py
"""
import json

from rdkit import Chem
from rdkit.Chem import AllChem
from eled import XTBCalculator, GenerateESPGrid, FitESPCharges, FitRESPCharges


def main():
    mol = Chem.AddHs(Chem.MolFromSmiles('O'))
    if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
        raise RuntimeError('Could not generate water coordinates')
    grid_settings = dict(scales=(1.4, 1.6, 1.8, 2.0), pointsPerAtom=100)
    points = GenerateESPGrid(mol, **grid_settings)
    result = XTBCalculator().Calculate(mol, espPointsAngstrom=points)
    fits = [FitESPCharges(result, equivalenceGroups=[(1, 2)]),
            FitRESPCharges(result, equivalenceGroups=[(1, 2)])]
    print(json.dumps({'grid_settings': grid_settings, 'points_angstrom': result.espPointsAngstrom,
                      'potential_au': result.espValuesAU,
                      'fits': [fit.ToDict() for fit in fits]}, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
