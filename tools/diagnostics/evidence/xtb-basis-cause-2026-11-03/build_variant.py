import json,os,subprocess,sys,hashlib
from pathlib import Path
root=Path.cwd(); stage=root/'build/xtb-orbital-investigation'
name=sys.argv[1];src=stage/(name+'-src');prefix=stage/(name+'-install');build=stage/(name+'-build')
reference=json.loads((root/'.xtb-rdkit/share/xtb/eled/build.json').read_text())
cache={}
for line in (root/'.xtb-rdkit/build/CMakeCache.txt').read_text().splitlines():
 if line and not line.startswith(('#','//')) and '=' in line:
  k,v=line.split('=',1);cache[k.split(':')[0]]=v
args=reference['configure_command'].copy()
args[args.index('-S')+1]=str(src);args[args.index('-B')+1]=str(build)
args=[('-DCMAKE_INSTALL_PREFIX='+str(prefix)) if s.startswith('-DCMAKE_INSTALL_PREFIX=') else s for s in args]
for key in ['CMAKE_C_COMPILER','CMAKE_Fortran_COMPILER','CMAKE_C_FLAGS','CMAKE_Fortran_FLAGS','CMAKE_SHARED_LINKER_FLAGS']:
 args.append('-D'+key+'='+cache[key])
env=os.environ.copy();env['PATH']=str(root/'.xtb-build-env/bin')+os.pathsep+env['PATH']
for k in ['LD_LIBRARY_PATH','LD_PRELOAD','XTBPATH','XTBHOME']:env.pop(k,None)
env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
commands=[args,['cmake','--build',str(build),'--parallel','6'],['ctest','--test-dir',str(build),'--output-on-failure','--parallel','6'],['cmake','--install',str(build)]]
for cmd in commands:
 print('COMMAND',cmd,flush=True);subprocess.run(cmd,env=env,check=True)
record=dict(source_commit=subprocess.check_output(['git','-C',str(src),'rev-parse','HEAD'],text=True).strip(),source_diff=subprocess.check_output(['git','-C',str(src),'diff'],text=True),commands=commands,library_sha256=hashlib.sha256((prefix/'lib/libxtb.so').read_bytes()).hexdigest())
(prefix/'diagnostic-build.json').write_text(json.dumps(record,indent=2)+'\n')
