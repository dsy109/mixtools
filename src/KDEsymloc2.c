#include <R.h>
#include <Rmath.h>

/* Slightly modified version of KDEsymloc.c that can handle an nxm matrix of
   current mean estimates (instead of merely an m-vector).  This is useful
   for the regression case, where the (i,j) mean depends not only on the
   jth component but also on the ith predictor value */
void KDEsymloc2(
     int *n, /* Sample size */
     int *m, /* Number of components */
     double *mu, /* nn*mm vector of current mean estimates */
     double *x, /* data:  vector of length n */
     double *h, /* bandwidth */
     double *z, /* nn*mm vector of normalized posteriors (or indicators in
                   stochastic case), normalized by "column" */
     double *f  /* KDE evaluated at n*m matrix of points, namely,
                   x_i - mu_ij for 1<=i<=n and 1<=j<=m   */ 
) {
  int nn=*n, mm=*m, i, j, a, b;
  double sum, u1, u2, tmp1, tmp2, hh=*h;
  double const1 = -1.0 / (2.0 * hh * hh);
  double const2 = 0.39894228040143267794/(2.0*hh*(double)nn); /* .3989...=1/(sqrt(2*pi)) */
  
  /* Must loop n^2*m^2 times; each f entry requires n*m calculations */
  for(a=0; a<nn; a++) {
    for(b=0; b<mm; b++) {
      sum = 0.0;
      u1 = (x[a]-mu[a + b*nn]);
      for(i=0; i<nn; i++) {
        for(j=0; j<mm; j++) {
          u2 = (x[i]-mu[i + j*nn]);
          tmp1 = u1-u2;
          tmp2 = -u1-u2;
          /* Use normal kernel */
          sum += z[i + j*nn] * (exp(tmp1 * tmp1 * const1) + 
                                exp(tmp2 * tmp2 * const1));
        }
      }
      f[a + b*nn] = sum * const2; /* (a,b) entry of f matrix */
    }
  }
}



