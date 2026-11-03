"""Fresh-handle C API water gradients versus central energy finite differences."""
import ctypes as ct,hashlib,json,sys
from pathlib import Path
import numpy as np
from rdkit_eth.backends.native import _NativeXTB
library=Path(sys.argv[1]).resolve();params=Path(sys.argv[2]).resolve();output=Path(sys.argv[3])
api=_NativeXTB(str(library));h=ct.c_void_p;dptr=ct.POINTER(ct.c_double)
get_gradient=api._bind('xtb_getGradient',None,[h,h,dptr]);accuracy=api._bind('xtb_setAccuracy',None,[h,h,ct.c_double])
coords=np.array([[0,0,-.73578586109551],[1.44183152868459,0,.36789293054775],[-1.44183152868459,0,.36789293054775]],dtype=float)
def evaluate(xyz,method):
 env=api.new_environment();molecule=calc=result=None
 try:
  api.set_verbosity(env,0);n=ct.c_int(3);z=(ct.c_int*3)(8,1,1);pos=(ct.c_double*9)(*xyz.ravel());q=ct.c_double(0);uhf=ct.c_int(0)
  molecule=api.new_molecule(env,ct.byref(n),z,pos,ct.byref(q),ct.byref(uhf),None,None);api.check(env,'new molecule')
  calc=api.new_calculator();getattr(api,'load_gfn'+str(method)+'xtb')(env,molecule,calc,str(params/f'param_gfn{method}-xtb.txt').encode());api.check(env,'load method')
  accuracy(env,calc,.0001);result=api.new_results();api.singlepoint(env,molecule,calc,result);api.check(env,'single point')
  energy=ct.c_double();gradient=(ct.c_double*9)();api.get_energy(env,result,ct.byref(energy));get_gradient(env,result,gradient);api.check(env,'results')
  return energy.value,np.array(gradient).reshape(3,3)
 finally:
  api._delete(api.del_results,result);api._delete(api.del_calculator,calc);api._delete(api.del_molecule,molecule);api._delete(api.del_environment,env)
records=[]
for method in [1,2]:
 energy,gradient=evaluate(coords,method);step=1e-4;fd=np.empty((3,3))
 for atom in range(3):
  for component in range(3):
   plus=coords.copy();minus=coords.copy();plus[atom,component]+=step;minus[atom,component]-=step
   fd[atom,component]=(evaluate(plus,method)[0]-evaluate(minus,method)[0])/(2*step)
 records.append(dict(method=f'GFN{method}-xTB',energy_hartree=energy,analytic_gradient_eh_bohr=gradient.tolist(),finite_difference_eh_bohr=fd.tolist(),max_abs_error=float(abs(gradient-fd).max()),analytic_gradient_norm=float(np.linalg.norm(gradient)),finite_difference_gradient_norm=float(np.linalg.norm(fd))))
 print(method,records[-1]['max_abs_error'],flush=True)
output.write_text(json.dumps(dict(library=str(library),sha256=hashlib.sha256(library.read_bytes()).hexdigest(),accuracy=.0001,step_bohr=step,positions_bohr=coords.tolist(),cases=records),indent=2)+'\n')
