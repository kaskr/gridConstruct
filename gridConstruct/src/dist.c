#include <R.h>
#include <Rdefines.h>
#include <R_ext/Error.h>

double dist(double lon1, double lat1, double lon2, double lat2){
  double r=6378.1;
  double G2R=M_PI/180;
  double cosd;
  cosd = sin(lat1 * G2R) * sin(lat2 * G2R) + 
    cos(lat1 * G2R) * cos(lat2 * G2R) * cos((lon1 - lon2) * G2R);
  if(cosd < -1)cosd = -1;
  if(cosd > 1)cosd = 1;
  return r * acos(cosd);
}
int whichmin_dist(double lon1, double lat1, double* lon2, double* lat2, int n){
  int index=0;
  double m=dist(lon1,lat1,lon2[0],lat2[0]);
  double tmp;
  for(int i=0;i<n;i++){
    tmp=dist(lon1,lat1,lon2[i],lat2[i]);  
    if(tmp<m){index=i;m=tmp;}
  }
  return index+1;
}
SEXP distkm ( SEXP lon1, SEXP lat1, SEXP lon2, SEXP lat2 ) {
  double *plon1=REAL(lon1);
  double *plon2=REAL(lon2);
  double *plat1=REAL(lat1);
  double *plat2=REAL(lat2);
  int n=LENGTH(lon1);
  if(n!=LENGTH(lon2))error("ups!");
  if(n!=LENGTH(lat1))error("ups!");
  if(n!=LENGTH(lat2))error("ups!");
  SEXP ans;
  PROTECT(ans=allocVector(REALSXP,n));
  double *pans=REAL(ans);
#pragma omp parallel for schedule(static)
  for(int i=0;i<n;i++){
    pans[i]=dist(plon1[i],plat1[i],plon2[i],plat2[i]);
  }
  UNPROTECT(1);
  return ans;
}
SEXP nearestkm ( SEXP lon1, SEXP lat1, SEXP lon2, SEXP lat2 ) {
  if(LENGTH(lon1)!=LENGTH(lat1))error("ups!");
  if(LENGTH(lon2)!=LENGTH(lat2))error("ups!");
  double *plon1=REAL(lon1);
  double *plon2=REAL(lon2);
  double *plat1=REAL(lat1);
  double *plat2=REAL(lat2);
  int n=LENGTH(lon1);
  int m=LENGTH(lon2);
  SEXP ans;
  PROTECT(ans=allocVector(INTSXP,n));
  int *pans=INTEGER(ans);
#pragma omp parallel for schedule(static)
  for(int i=0;i<n;i++){
    pans[i]=whichmin_dist(plon1[i],plat1[i],plon2,plat2,m);
  }
  UNPROTECT(1);
  return ans;
}

