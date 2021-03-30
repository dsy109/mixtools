#include <R.h>
#include <Rmath.h>

/* Implement symmetric kernel density-estimation step for location 
   mixture model, equation (20) in Benaglia et al (2008) */
void KDEsymloc(
     int *n, /* Sample size */
     int *m, /* Number of components */
     double *mu, /* m-vector of current mean estimates */
     double *x, /* data:  vector of length n */
     double *h, /* bandwidth */
     double *z, /* nn*mm vector of normalized posteriors (or indicators in
                   stochastic case), normalized by "column" */
     double *f  /* KDE evaluated at n*m matrix of points, namely,
                   x_i - mu_j for 1<=i<=n and 1<=j<=m   */ 
) {
  int nn=*n, mm=*m, i, j, a, b;
  double sum, u1, u2, tmp1, tmp2, hh=*h;
  double const1 = -1.0 / (2.0 * hh * hh);
  double const2 = 0.39894228040143267794/(2.0*hh*(double)nn); /* .3989...=1/(sqrt(2*pi)) */
  
  /* Must loop n^2*m^2 times; each f entry requires n*m calculations */
  for(a=0; a<nn; a++) {
    for(b=0; b<mm; b++) {
      sum = 0.0;
      u1 = (x[a]-mu[b]);
      for(i=0; i<nn; i++) {
        for(j=0; j<mm; j++) {
          u2 = (x[i]-mu[j]);
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



