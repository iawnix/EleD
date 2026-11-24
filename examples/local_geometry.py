"""Selected substituent sterics and ring geometry using owned conformers."""
import json

from eled import (
    CalcBuriedVolume, CalcHOMA, CalcPlaneDescriptors, CalcRingPuckering,
    CalcSterimol, PrepareConformers,
)


def prepared(smiles):
    result = PrepareConformers(smiles)
    if not result.success:
        raise RuntimeError(result.report)
    return result.molecule


def main():
    anisole = prepared('COc1ccccc1')
    # SMILES order gives methyl C=0, O=1, and the attached ring C=2.
    methyl_hydrogens = [a.GetIdx() for a in anisole.GetAtomWithIdx(0).GetNeighbors()
                       if a.GetAtomicNum() == 1]
    substituent = [0, 1, *methyl_hydrogens]
    sterimol = CalcSterimol(anisole, substituent, axis=(2, 1))
    occupied_atoms = [a.GetIdx() for a in anisole.GetAtoms() if a.GetIdx() != 2]
    volume = CalcBuriedVolume(anisole, occupied_atoms, center=2)

    naphthalene = prepared('c1ccc2ccccc2c1')
    rings = list(naphthalene.GetRingInfo().AtomRings())
    homa = [CalcHOMA(naphthalene, ring) for ring in rings]
    planes = CalcPlaneDescriptors(naphthalene, rings[0], other_atom_indices=rings[1])
    cyclohexane = prepared('C1CCCCC1')
    puckering = CalcRingPuckering(cyclohexane, cyclohexane.GetRingInfo().AtomRings()[0])
    print(json.dumps(dict(sterimol=sterimol['values'], buried_volume=volume['values'],
                          homa=[r['values'] for r in homa], planes=planes['values'],
                          puckering=puckering['values']), indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
