#!/usr/bin/env python3
"""Standard-library-only installed SCC ABI smoke test; emits JSON evidence."""
import ctypes as c
import json
from pathlib import Path
import sys


def verify(prefix):
    lib=c.CDLL(str(prefix/'lib/libxtb.so')); h=c.c_void_p
    def bind(name, result, args):
        f=getattr(lib,'xtb_'+name);f.restype=result;f.argtypes=args;return f
    version=bind('rdkit_getSCFAPIVersion',c.c_int,[])
    if version()!=1: raise RuntimeError('unexpected SCC ABI')
    newenv=bind('newEnvironment',h,[]);env=newenv()
    bind('setVerbosity',None,[h,c.c_int])(env,0)
    check=bind('checkEnvironment',c.c_int,[h])
    create=bind('newMolecule',h,[h,c.POINTER(c.c_int),c.POINTER(c.c_int),c.POINTER(c.c_double),
                                c.POINTER(c.c_double),c.POINTER(c.c_int),h,h])
    nat=c.c_int(3);at=(c.c_int*3)(8,1,1);xyz=(c.c_double*9)(0,0,0,0,0,1.8,1.7,0,-.6)
    q=c.c_double(0);spin=c.c_int(0)
    mol=create(env,c.byref(nat),at,xyz,c.byref(q),c.byref(spin),None,None)
    newcalc=bind('newCalculator',h,[]);newres=bind('newResults',h,[])
    sp=bind('singlepoint',None,[h,h,h,h]);copy=bind('copyResults',h,[h])
    get=bind('rdkit_getSCFDiagnostics',c.c_int,[h,h,c.POINTER(c.c_int),c.POINTER(c.c_double)])
    def delete(kind,value):
        if value: bind('del'+kind,None,[c.POINTER(h)])(c.byref(h(value)))
    def read(result):
        n=c.c_int(17);data=(c.c_double*4)(*([19.]*4))
        available=get(env,result,c.byref(n),data)
        return available,n.value,list(data)
    rows=[]
    try:
        for method in (0,1,2):
            calc=newcalc();res=newres();copied=None
            try:
                assert read(res)==(0,17,[19.]*4)
                loader=bind('loadGFN'+str(method)+'xTB',None,[h,h,h,c.c_char_p])
                loader(env,mol,calc,str(prefix/f'share/xtb/param_gfn{method}-xtb.txt').encode())
                sp(env,mol,calc,res)
                if check(env): raise RuntimeError('singlepoint failed')
                available,n,data=read(res)
                if method:
                    assert available==1 and n>=1 and 0<=data[0]<data[2] and 0<=data[1]<data[3]
                    copied=copy(res);assert read(copied)==read(res)
                    bind('setMaxIter',None,[h,h,c.c_int])(env,calc,1)
                    sp(env,mol,calc,res);assert check(env)
                    assert read(res)==(0,17,[19.]*4)
                    assert read(copied)==(available,n,data)
                    buf=c.create_string_buffer(4096);size=c.c_int(4096)
                    bind('getError',None,[h,c.c_char_p,c.POINTER(c.c_int)])(env,buf,c.byref(size))
                else: assert (available,n,data)==(0,17,[19.]*4)
                rows.append({'method':method,'available':bool(available),'iterations':n if available else None,
                             'values':data if available else None,'copy_and_failure_invalidation_checked':bool(method)})
            finally:
                delete('Results',copied);delete('Results',res);delete('Calculator',calc)
    finally:
        delete('Molecule',mol);delete('Environment',env)
    return rows


if __name__=='__main__': print(json.dumps(verify(Path(sys.argv[1]).resolve())))
