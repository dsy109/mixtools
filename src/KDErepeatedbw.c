#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>

/* simultaneously calculate m different products of KDEs, 1 for each component,
   as in equation (8) of Benaglia et al for a fixed value of \ell.  
   If r is the number of coordinates in block \ell, then each
   final value is the r-fold product of KDEs */
void KDErepeatedbw(
     int *nn, /* Sample size */
     int *mm, /* Number of components */
     int *rr, /* size of current block */
     double *x, /* data:  vector of length nn*rr */
     double *hh, /* m-vector of bandwidths (compare to KDErepeated) */
     double *z, /* nn*mm vector of normalized posteriors (or indicators in
                   stochastic case), normalized by "column" */
     double *f  /* nxm matrix of KDE products */
) {
  int n=*nn, i, ii;
  int mn = *mm*n, rn=*rr*n, jn, kn, kkn;
  double sum1, sum2, tmp, xik;
  double const2 = 0.39894228040143267794/((double)(*rr)); /* .3989...=1/(sqrt(2*pi)) */
  double const1;

  for(jn=0; jn<mn; jn+=n, hh++) { /* jn is component index times n */
    /* at each iteration, *hh changes to the current bandwidth value */
    const1 = -0.5 / (*hh * *hh);
    for(i=0; i<n; i++) {
      f[i + jn] = 1.0;
      for(kn=0; kn<rn; kn+=n) { /* kn is coordinate index within block times n */
        sum1 = 0.0;
        xik = x[i + kn];
        for(ii=0; ii<n; ii++) { /* ii is the i index of equation (8) */
          sum2 = 0.0;
          for(kkn=0; kkn<rn; kkn+=n) { /* kk is the k index of equation (8) */
            tmp = xik - x[ii + kkn];
            sum2 +=  exp(tmp * tmp * const1); /* Using normal kernel */
          }
          sum1 += z[ii + jn] * sum2;
        }
        f[i + jn] *= sum1 * const2 / *hh;
      }
    }
  }
}



