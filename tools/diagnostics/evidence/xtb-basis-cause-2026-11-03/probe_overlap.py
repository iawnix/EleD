import json,sys,numpy as np
from pathlib import Path
sys.path.insert(0,str(Path.cwd()/'tools'))
from compare_orbitals import molecule
from eled import rdXTBTools
from eled.backends.tblite import run_xtb_tblite
from eled.orbital_alignment import libxtb_to_tblite_rows
fixtures=json.load(open('tests/data/orbital_alignment_cases.json'))['cases']
records=[]
for name in ['water','bromomethane','iodomethane']:
 case=next(c for c in fixtures if c['name']==name);m=molecule(case)
 ok,n=rdXTBTools.RunMol(m,method='GFN1-xTB',xtbLibrary=sys.argv[1],accuracy=.01,keepOrbitalCoefficients=True)
 assert ok,n.errorMessage
 t=run_xtb_tblite(m,method='GFN1-xTB',accuracy=.01,save_integrals=True);b=t.raw_json
 p,sg,labels=libxtb_to_tblite_rows(case['atomic_numbers'],n.parameterFile,b['orbital_map'],b['shell_map'],b['angular_momenta'],n.numOrbitals)
 c=n.GetOrbitalCoefficients()[p,:]*sg[:,None];ic=np.linalg.inv(c);inferred=ic.T@ic;s=np.array(b['overlap_matrix']);a=np.array([l['atom'] for l in labels]);same=a[:,None]==a[None,:]
 d=abs(inferred-s);rec=dict(name=name,library=n.libraryPath,library_sha256=n.librarySHA256,coefficient_condition_number=float(np.linalg.cond(c)),same_atom_overlap_difference=float(d[same].max()),different_atom_overlap_difference=float(d[~same].max()),tblite_diagonal=s.diagonal().tolist(),inferred_xtb_diagonal=inferred.diagonal().tolist(),labels=labels)
 records.append(rec)
 print(name,rec['same_atom_overlap_difference'],rec['different_atom_overlap_difference'])
Path(sys.argv[2]).write_text(json.dumps(dict(note='S inferred as inv(C).T @ inv(C) from full real eigenvectors; not a directly exposed C API overlap and not an independent orthonormality check.',cases=records),indent=2)+'\n')
