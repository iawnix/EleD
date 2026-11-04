"""Persistent batches: rerun with the same output directory to resume."""
import argparse
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit_eth import (ReactivityJob, ReactivityCheckpoint, XTBDiskCache,
                       run_reactivity_batch, export_descriptor_reports)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--retry-failed',action='store_true')
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    jobs=[]
    for name,smi,spins in [('water','O',{}),('ethanol','CCO',{}),
                          ('methyl','[CH3]',{'ionized_uhf':0,'attached_uhf':0})]:
        mol=Chem.AddHs(Chem.MolFromSmiles(smi))
        if AllChem.EmbedMolecule(mol,randomSeed=1):raise RuntimeError('embedding failed')
        jobs.append(ReactivityJob(name,mol,**spins))
    jobs.append(ReactivityJob('failed-upstream-parse',None))
    checkpoint=ReactivityCheckpoint(args.output/'dataset.sqlite')
    rows=run_reactivity_batch(jobs,cache=XTBDiskCache(args.output/'cache.sqlite'),
                              checkpoint=checkpoint,retry_failed=args.retry_failed)
    for kind in ('csv','jsonl'):export_descriptor_reports(rows,args.output/f'descriptors.{kind}')
    for row in rows:
        print(row['record_id'],row['quality']['status'],row['execution']['source'])


if __name__=='__main__':main()
