#include<math.h>

/* Translated from FORTRAN code written by Fengjuan Xuan */
void mudepth (int *nn, int *tt, int *dd, double *mpt, 
              double *x, int *count, double *sdep) {
  int n=*nn, t=*tt, d=*dd;                
  int i, j, k, l;
  double d1, d2, d3, d5, xik, xjk, mptlk;
  
  for (l=0; l<t; l++) {
    count[l] = 0;
    sdep[l] = 0.0;
    for (i=0; i<(n-1); i++) {
      for (j=i+1; j<n; j++) {
        d1 = d2 = d3 = 0.0;
        for (k=0; k<d; k++) {
          xik = x[i+k*n];
          xjk = x[j+k*n];
          mptlk = mpt[l+k*t];
          d1 = d1 + (xik-mptlk)*(xik-mptlk);
          d2 = d2 + (xjk-mptlk)*(xjk-mptlk);
          d3 = d3 + (xik-xjk)*(xik-xjk);
        }
        if (d1 + d2 - d3 <= 0.0) count[l]++;
      }
    }
    d5 = sqrt((double) n*(n-1)/8); 
    sdep[l]=(count[l]-n*(n-1)/4)/d5;   
  }
}

