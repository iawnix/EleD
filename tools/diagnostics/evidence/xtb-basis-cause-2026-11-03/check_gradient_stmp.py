import os,json,subprocess,hashlib,shutil
from pathlib import Path
root=Path.cwd();stage=root/'build/xtb-orbital-investigation';src=stage/'release-src';build=stage/'release-build';out=stage/'gradient-stmp';out.mkdir(exist_ok=True)
p=src/'src/xtb/hamiltonian.f90';original=p.read_text()
start=original.index('subroutine build_dSDQH0_noreset(');end=original.index('end subroutine build_dSDQH0_noreset',start)
section=original[start:end]
assert section.count('g_xyz(ixyz) = g_xyz(ixyz)+stmp+dtmp+qtmp')==1
section=section.replace('real(wp) stmp,ral(3,3)','real(wp) ral(3,3)').replace('tmp,stmp,dtmp,qtmp','tmp,dtmp,qtmp').replace('g_xyz(ixyz) = g_xyz(ixyz)+stmp+dtmp+qtmp','g_xyz(ixyz) = g_xyz(ixyz)+dtmp+qtmp')
patched=original[:start]+section+original[end:]
env=os.environ.copy();env['PATH']=str(root/'.xtb-build-env/bin')+os.pathsep+env['PATH'];env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
for key in ['XTBPATH','XTBHOME','LD_LIBRARY_PATH','LD_PRELOAD']:env.pop(key,None)
try:
 p.write_text(patched)
 (out/'source.patch').write_text(subprocess.check_output(['git','-C',str(src),'diff'],text=True))
 with (out/'build.log').open('w') as log:subprocess.run(['cmake','--build',str(build),'--parallel','6'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 shutil.copyfile(build/'libxtb.so',out/'libxtb.so')
 with (out/'ctest.log').open('w') as log:r=subprocess.run(['ctest','--test-dir',str(build),'--output-on-failure','--parallel','6'],env=env,stdout=log,stderr=subprocess.STDOUT)
 (out/'record.json').write_text(json.dumps(dict(ctest_exit_code=r.returncode,library_sha256=hashlib.sha256((out/'libxtb.so').read_bytes()).hexdigest()),indent=2)+'\n')
 print('ctest exit',r.returncode,flush=True)
finally:p.write_text(original)
