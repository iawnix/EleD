#!/usr/bin/env python3
"""Frozen-corpus reactivity, checkpoint recovery and interface comparison.

Retains failures and initialization differences. Protocol regression, not an
IP/EA accuracy benchmark. Requires the compiled adapter for both paths.
"""
import argparse
from collections import Counter
from dataclasses import replace
import json
from pathlib import Path
from unittest.mock import patch

from rdkit import Chem
from eled import (ReactivityJob, XTBDiskCache, XTBOptions, RunReactivityBatch,
                       ExportDescriptorReports)


def load_jobs(path, method):
    jobs=[]
    for case in json.loads(Path(path).read_text())['cases']:
        mol=Chem.AddHs(Chem.MolFromSmiles(case['smiles']))
        assert [a.GetAtomicNum() for a in mol.GetAtoms()]==case['atomic_numbers']
        conf=Chem.Conformer(mol.GetNumAtoms());conf.Set3D(True)
        for i,xyz in enumerate(case['positions_angstrom']):conf.SetAtomPosition(i,xyz)
        mol.AddConformer(conf)
        jobs.append(ReactivityJob(case['name'],mol,
            options=XTBOptions(method=method,charge=case['charge'],uhf=case['uhf'],accuracy=.1),
            ionized_uhf=case['ionized_uhf'],attached_uhf=case['attached_uhf']))
    for name in ['ethanol','ammonium']:
        job=next(j for j in jobs if j.record_id==name)
        jobs.append(replace(job,record_id=name+'-cold',options=replace(job.options,electronic_temperature=0.)))
    job=next(j for j in jobs if j.record_id=='ethanol')
    jobs.append(replace(job,record_id='ethanol-one-iteration',options=replace(job.options,max_iterations=1)))
    return jobs


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cases',type=Path,default=Path(__file__).resolve().parents[1]/'tests/data/reactivity_cases.json')
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();args.output.mkdir(parents=True,exist_ok=True)
    cache=XTBDiskCache(args.output/'cache.sqlite');all_rows={};summaries={}
    for method in ('gfn1','gfn2'):
        jobs=load_jobs(args.cases,method)
        for initialization in ('reference','independent'):
            for interface in ('native','compiled'):
                label=f'{method}-{initialization}-{interface}'
                checkpoint=args.output/f'{label}.sqlite'
                rows=RunReactivityBatch(jobs,interface=interface,initialization=initialization,
                                           checkpoint=checkpoint,cache=cache)
                with patch('eled.reactivity_batch.RunReactivity',side_effect=AssertionError('resume recalculated')):
                    resumed=RunReactivityBatch(jobs,interface=interface,initialization=initialization,
                                                checkpoint=checkpoint,cache=cache)
                assert all(r['execution']['source']=='checkpoint' for r in resumed)
                for kind in ('csv','jsonl'):ExportDescriptorReports(rows,args.output/f'{label}.{kind}')
                summaries[label]=dict(Counter('failed' if r['error'] else 'completed' for r in rows))
                all_rows[label]=rows
    differences=[]
    for method in ('gfn1','gfn2'):
        for initialization in ('reference','independent'):
            a=all_rows[f'{method}-{initialization}-native'];b=all_rows[f'{method}-{initialization}-compiled']
            for ra,rb in zip(a,b):
                row={'method':method,'initialization':initialization,'id':ra['record_id']}
                if ra['error'] or rb['error']:
                    row['both_failed']=bool(ra['error'] and rb['error'])
                    row['failed_states']=[(r['error'] or {}).get('state') for r in (ra,rb)]
                    row['matched']=row['both_failed'] and len(set(row['failed_states']))==1
                else:
                    row['max_abs_energy_difference_hartree']=max(abs(ra['states'][k]['energy_hartree']-
                        rb['states'][k]['energy_hartree']) for k in ra['states'])
                    row['max_abs_fukui_difference']=max(abs(x[k]-y[k]) for x,y in zip(ra['atoms'],rb['atoms'])
                        for k in ('f_plus','f_minus','f_radical','dual'))
                    row['matched']=row['max_abs_energy_difference_hartree']<1e-8 and row['max_abs_fukui_difference']<1e-7
                differences.append(row)
    initialization_differences=[]
    for method in ('gfn1','gfn2'):
        for a,b in zip(all_rows[f'{method}-reference-native'],all_rows[f'{method}-independent-native']):
            row={'method':method,'id':a['record_id'],'failed_states':[(r['error'] or {}).get('state') for r in (a,b)]}
            if a['error'] is None and b['error'] is None:
                row['max_abs_fukui_difference']=max(abs(x[k]-y[k]) for x,y in zip(a['atoms'],b['atoms'])
                    for k in ('f_plus','f_minus','f_radical','dual'))
            initialization_differences.append(row)
    summary={'protocols':summaries,'interface_comparisons':differences,
             'initialization_comparisons':initialization_differences,
             'all_interface_checks_passed':all(d['matched'] for d in differences),
             'scope':'frozen unoptimized structures; failure accounting and protocol agreement only'}
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    print(json.dumps({'protocols':summaries,'all_interface_checks_passed':summary['all_interface_checks_passed']},indent=2))
    return int(not summary['all_interface_checks_passed'])


if __name__=='__main__':raise SystemExit(main())
