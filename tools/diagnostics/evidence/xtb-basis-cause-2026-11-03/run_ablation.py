from pathlib import Path
import hashlib,json,os,subprocess,shutil
root=Path.cwd();stage=root/'build/xtb-orbital-investigation';source=stage/'release-src';build=stage/'release-build'
env=os.environ.copy();env['PATH']=str(root/'.xtb-build-env/bin')+os.pathsep+env['PATH'];env.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',XTBPATH=str(source))
for k in ['LD_LIBRARY_PATH','LD_PRELOAD','XTBHOME']:env.pop(k,None)
repo=stage/'source'
def git(ref,path):return subprocess.check_output(['git','-C',str(repo),'show',ref+':'+path],text=True)
original_h=git('v6.7.1','src/xtb/hamiltonian.f90'); original_s=git('v6.7.1','src/slater.f90');new_s=git('a9040512','src/slater.f90')
# Experimental subset of upstream #1475, for single-point causal isolation.
# Not a deployable backport: it deliberately leaves other upstream edits out.
onsite_h=original_h.replace('            sint(i,i) = 1.0_wp + sint(i,i)\n            H0(ii) = H0(ii) + selfEnergy(ish, iat)\n','')
onsite_h=onsite_h.replace('            !call dtrf2(ss,ishtyp,jshtyp)','            call dtrf2(ss,ishtyp,jshtyp)')
needle='                  dpint(1:3, iao, jao) = dpint(1:3, iao, jao) + dd(1:3, jj, ii)'
assert original_h.count(needle)==1
onsite_h=onsite_h.replace(needle,'''                  sint(iao, jao) = sint(iao, jao) + ss(jj, ii)
                  if (iao /= jao) sint(jao, iao) = sint(jao, iao) + ss(jj, ii)
                  H0(lin(iao, jao)) = H0(lin(iao, jao)) &
                     & + 0.5_wp * (selfEnergy(ish, iat) + selfEnergy(jsh, iat)) * ss(jj, ii)
'''+needle)
for name,h,s in [('tag-release',original_h,original_s),('onsite-only',onsite_h,original_s),('onsite-and-5sp',onsite_h,new_s),('5sp-only',original_h,new_s)]:
 (source/'src/xtb/hamiltonian.f90').write_text(h);(source/'src/slater.f90').write_text(s)
 out=stage/name;out.mkdir(exist_ok=True)
 with (out/'build.log').open('w') as log:subprocess.run(['cmake','--build',str(build),'--parallel','6'],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 # Freeze each library before the next incremental rebuild. The RPATH retains
 # identical dependency libraries; each comparison starts a fresh process.
 library=out/'libxtb.so';shutil.copyfile(build/'libxtb.so',library)
 diff=subprocess.check_output(['git','-C',str(source),'diff'],text=True)
 (out/'source.patch').write_text(diff)
 record=dict(base_commit=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip(),library_sha256=hashlib.sha256(library.read_bytes()).hexdigest(),source_patch_sha256=hashlib.sha256(diff.encode()).hexdigest(),diagnostic_only=True)
 cmd=[str(root/'.rdxtb-build-env/bin/python'),str(root/'tools/compare_orbitals.py'),'--library',str(library),'--accuracy','0.01','--case','water','--case','bromomethane','--case','iodomethane','--output',str(out/'comparison')]
 with (out/'comparison.log').open('w') as log:
  result=subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
 record['comparison_command']=cmd;record['comparison_exit_code']=result.returncode
 (out/'build-record.json').write_text(json.dumps(record,indent=2)+'\n')
 print(name,result.returncode,flush=True)
 r=json.loads((out/'comparison/summary.json').read_text())
 for c in r['cases']:
  print(c['name'],c['method'],c.get('differences',{}).get('energy_eh'),c.get('differences',{}).get('candidate_orthonormality_error'),flush=True)
# Leave the research checkout at its original source state, retaining variants.
(source/'src/xtb/hamiltonian.f90').write_text(original_h);(source/'src/slater.f90').write_text(original_s)
