/* Verify the project-local basis ABI against the installed library. */
#include <xtb.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define REQUIRE(x) do { if (!(x)) { fprintf(stderr,"basis verification failed: %s at %d\n",#x,__LINE__); exit(1); } } while(0)
static void check(xtb_TEnvironment env) {
  if (xtb_checkEnvironment(env)) { char text[4096]; int n=4096; xtb_getError(env,text,&n); fprintf(stderr,"%s\n",text); exit(1); }
}
static void expected_error(xtb_TEnvironment *env) {
  REQUIRE(xtb_checkEnvironment(*env)); xtb_delEnvironment(env);
  *env=xtb_newEnvironment(); xtb_setVerbosity(*env,0);
}
static double verify(xtb_TEnvironment env, xtb_TResults res, int method) {
  int sizes[3]={0}; xtb_rdkit_getBasisSize(env,res,sizes); check(env);
  int n=sizes[2], nsh=sizes[1]; REQUIRE(sizes[0]==3 && n>0 && nsh>0);
  double *s=calloc(n*n,sizeof(double)), *c=calloc(n*n,sizeof(double));
  double *sc=calloc(n*n,sizeof(double)), *occ=calloc(n,sizeof(double));
  int *aa=calloc(n,sizeof(int)), *as=calloc(n,sizeof(int));
  int *sa=calloc(nsh,sizeof(int)), *sl=calloc(nsh,sizeof(int));
  REQUIRE(s && c && sc && occ && aa && as && sa && sl);
  xtb_rdkit_getOverlap(env,res,n*n,s);
  xtb_rdkit_getBasisMap(env,res,n,nsh,aa,as,sa,sl);
  xtb_getOrbitalCoefficients(env,res,c); xtb_getOrbitalOccupations(env,res,occ); check(env);
  for(int i=0;i<n;i++) {
    REQUIRE(aa[i]>=0 && aa[i]<3 && as[i]>=0 && as[i]<nsh && sa[as[i]]==aa[i]);
    for(int j=0;j<n;j++) { REQUIRE(isfinite(s[i+j*n]) && fabs(s[i+j*n]-s[j+i*n])<1e-12);
      for(int k=0;k<n;k++) sc[i+j*n]+=s[i+k*n]*c[k+j*n]; }
  }
  for(int sh=0;sh<nsh;sh++) {
    int count=0; for(int i=0;i<n;i++) count+=(as[i]==sh);
    REQUIRE(sl[sh]>=0 && sl[sh]<=4 && count==2*sl[sh]+1);
  }
  double error=0., population[3]={0}, q[3];
  for(int i=0;i<n;i++) for(int j=0;j<n;j++) {
    double v=0.; for(int k=0;k<n;k++) v+=c[k+i*n]*sc[k+j*n];
    error=fmax(error,fabs(v-(i==j)));
  }
  REQUIRE(error<1e-8);
  for(int a=0;a<n;a++) for(int i=0;i<n;i++) population[aa[a]]+=c[a+i*n]*sc[a+i*n]*occ[i];
  REQUIRE(fabs(population[0]+population[1]+population[2]-8.)<1e-7);
  if(method>0) {
    xtb_getCharges(env,res,q); check(env);
    for(int a=0;a<3;a++) REQUIRE(fabs((a==0?6.:1.)-population[a]-q[a])<1e-7);
  }
  free(s);free(c);free(sc);free(occ);free(aa);free(as);free(sa);free(sl);
  return error;
}
int main(int argc,char **argv) {
  REQUIRE(argc==2 && xtb_rdkit_getBasisAPIVersion()==1);
  int n=3, at[]={8,1,1}, uhf=0; double charge=0.;
  double xyz[]={0.,0.,0.,1.80884585,0.,0.,-.453984,1.752338,0.};
  xtb_TEnvironment env=xtb_newEnvironment(); xtb_setVerbosity(env,0);
  xtb_TMolecule mol=xtb_newMolecule(env,&n,at,xyz,&charge,&uhf,NULL,NULL); check(env);
  xtb_TCalculator calc=xtb_newCalculator(); xtb_TResults res=xtb_newResults();
  int sizes[3]={-1,-1,-1};
  xtb_rdkit_getBasisSize(env,res,sizes); REQUIRE(sizes[0]==-1); expected_error(&env);
  double worst=0.;
  for(int method=0;method<3;method++) {
    char param[4096]; snprintf(param,sizeof(param),"%s/param_gfn%d-xtb.txt",argv[1],method);
    if(method==0) xtb_loadGFN0xTB(env,mol,calc,param);
    else if(method==1) xtb_loadGFN1xTB(env,mol,calc,param);
    else xtb_loadGFN2xTB(env,mol,calc,param);
    check(env);
    xtb_rdkit_setKeepOverlap(env,res,1); xtb_singlepoint(env,mol,calc,res); check(env);
    worst=fmax(worst,verify(env,res,method));
    xtb_TResults copy=xtb_copyResults(res);
    double sentinel=123.; xtb_rdkit_getOverlap(env,res,1,&sentinel);
    REQUIRE(sentinel==123.); expected_error(&env);
    xtb_rdkit_getOverlap(env,res,0,NULL); expected_error(&env);
    xtb_rdkit_setKeepOverlap(env,res,0);
    xtb_rdkit_getBasisSize(env,res,sizes); expected_error(&env);
    verify(env,copy,method); xtb_delResults(&copy);
    xtb_singlepoint(env,mol,calc,res); check(env);
    xtb_rdkit_getBasisSize(env,res,sizes); expected_error(&env);
  }
  /* A failed SCF on a reused result must not expose the previous snapshot. */
  xtb_rdkit_setKeepOverlap(env,res,1); xtb_singlepoint(env,mol,calc,res); check(env);
  xyz[3]+=0.2; xtb_updateMolecule(env,mol,xyz,NULL); check(env);
  xtb_setMaxIter(env,calc,1); xtb_singlepoint(env,mol,calc,res); expected_error(&env);
  xtb_rdkit_getBasisSize(env,res,sizes); expected_error(&env);
  xtb_setMaxIter(env,calc,250); xtb_singlepoint(env,mol,calc,res); check(env); verify(env,res,2);
  /* Finite-difference workflows must invalidate the single-point snapshot. */
  double hessian[81];
  xtb_hessian(env,mol,calc,res,hessian,NULL,NULL,NULL,NULL); check(env);
  xtb_rdkit_getBasisSize(env,res,sizes); expected_error(&env);
  xtb_singlepoint(env,mol,calc,res); check(env); verify(env,res,2);
  xtb_delResults(&res); xtb_delCalculator(&calc); xtb_delMolecule(&mol); xtb_delEnvironment(&env);
  printf("{\"basis_api_version\":1,\"methods\":3,\"max_orthonormality_error\":%.15g,\"lifecycle_checks\":true}\n",worst);
  return 0;
}
