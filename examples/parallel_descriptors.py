"""Run with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 in the installed runtime."""
from eled import prepare_conformers, run_descriptor_batch


def main():
    molecules = []
    for smiles in ('O', 'CCO', 'c1ccccc1'):
        prepared = prepare_conformers(smiles)
        if not prepared.success: raise RuntimeError(prepared.report['error'])
        molecules.append(prepared.molecule)
    reports = run_descriptor_batch(molecules, record_ids=['water', 'ethanol', 'benzene'],
                                   max_workers=2, timeout=30)
    for report in reports:
        if report['error']: raise RuntimeError(report['error'])
        print(report['record_id'], report['values']['xtb_gap_ev'], 'eV',
              report['quality']['scf_details'])


if __name__ == '__main__': main()
