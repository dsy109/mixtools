#include <Rmath.h>
#include <R.h>
#include <Rinternals.h>

void newz (int *n, int *k, double *V, double *W, double *z) {
  int i, j, l, ind1, ind2, nn=*n, kk=*k;
  double sum;
  
  for(i=0; i< nn; i++) {
    for(j=0; j< kk; j++) {
      sum=1.0;
      ind1 = i + nn * j;
      for(l = 0; l < kk; l++) {
        if (l != j) {
          ind2 = i + nn * l;
          sum += V[ind2]/V[ind1]*exp(W[ind1]-W[ind2]);
        }
      }
      z[ind1] = 1.0/sum;  
    }
  } 
}

