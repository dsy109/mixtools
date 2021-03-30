#include <R.h>
#include <Rmath.h>

/* Implement symmetric kernel density-estimation step for location 
   mixture model with 1 component known , as in 
   Bordes Delmas & Vandekerkhove (2006)
   
   normalization before symmetrization must be 1/(n\lambda_2)h
*/
void KDEsymloc1comp(
     int *n, /* Sample size */
     double *mean, /* component 2 mean estimate, scalar */
     double *lambda, /* component 2 weigtht, scalar */
     double *x, /* data:  vector of length n */
     double *h, /* bandwidth */
     double *z, /* nn*mm vector of normalized posteriors (or indicators in
                   stochastic case), normalized by "column" */
     double *f  /* KDE evaluated at n vector of points, namely,
                   x_i - mu  for 1<=i<=n  */ 
) {
  int nn=*n, i, j, a;
  double sum, u1, u2, tmp1, tmp2, hh=*h, mu=*mean, lbd=*lambda;
  double const1 = -1.0 / (2.0 * hh * hh);
  double const2 = 0.39894228040143267794/(2.0*hh*(double)nn*lbd); 
  				/* .3989...=1/(sqrt(2*pi)) */
  
  /* loop over each f entry evaluated at x_a - mu */
  for(a=0; a<nn; a++) {
      sum = 0.0;
      u1 = x[a]-mu;
      for(i=0; i<nn; i++) {
 	  j = 1;  /* only component 2 = j index 1 is considered */
          u2 = x[i] - mu;
          tmp1 = u1 - u2;
          tmp2 = - u1 - u2;
          /* Use normal kernel and symmetrization */
          
          sum += z[i + j*nn] * (exp(tmp1 * tmp1 * const1) + 
                                exp(tmp2 * tmp2 * const1));
      }
      f[a] = sum * const2; /* a entry of f vector at x_a-mu */
  }
}



