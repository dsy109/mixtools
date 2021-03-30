#include <R.h>
#include <Rmath.h>

/* Implement kernel density-estimation step for location-scale
   mixture model (17) in Benaglia et al (2008) where each  
   component and block has the same shape but may differ from
   others by a location and scale, implementing equation (18) */
void KDElocscale(
     int *nn, /* Sample size */
     int *mm, /* Number of components */
     int *rr, /* Number of coordinates */
     int *blockid, /* r-vector of block numbers */
     double *mu, /* m by max(blockid) matrix of current mean estimates */
     double *sigma, /* m by max(blockid) matrix of current stdev estimates */
     double *x, /* n by r data matrix */
     double *hh, /* scalar bandwidth */
     double *z, /* n by m vector of normalized posteriors (or indicators in
                   stochastic case), normalized by "column" */
     double *f  /* n by m matrix of KDE products */
) {
  int n=*nn, m=*mm, r=*rr;
  int i, j, k, ii, kk, ell, ell2;
  double sum1, sum2, tmp, h=*hh, u;
  double const1 = -0.5 / (h * h);
  double const2;

  for(j=0; j<m; j++) {
    /* In const2 expression, .3989... equals 1/(sqrt(2*pi)) */
    const2 = 0.39894228040143267794/(h*sigma[j]*(double)r); 
    for (i=0; i<n; i++) {
      f[i + j*n] = 1.0;
      for (k=0; k<r; k++) {
        ell = blockid[k]-1;
        u = (x[i + k*n] - mu[j + m*ell])/sigma[j + m*ell];
        sum1 = 0.0;
        for (ii=0; ii<n; ii++) {
          sum2 = 0.0;
          for (kk=0; kk<r; kk++) {
            ell2 = blockid[kk]-1;
            tmp = (u - x[ii+kk*n] + mu[j + m*ell2])/sigma[j + m*ell2];
            sum2 += exp(tmp * tmp * const1);
          }
          sum1 += z[ii + j*n] * sum2;
        }
        f[i + j*n] *= sum1 * const2;
      }
    }
  }
}
      


//THE FOLLOWING IS NOT FINISHED!
//
///* Implement kernel density-estimation step for location-scale
//   mixture model with a different shape for each component, 
//   equation (11) in Benaglia et al (2008) */
//void KDElocscale(
//     int *nn, /* Sample size */
//     int *mm, /* Number of components */
//     int *rr, /* size of current block */
//     double *mu, /* m-vector of current mean estimates */
//     double *sigma, /* m-vector of current stdev estimates */
//     double *x, /* data:  vector of length n */
//     double *hh, /* scalar bandwidth */
//     double *z, /* nn*mm vector of normalized posteriors (or indicators in
//                   stochastic case), normalized by "column" */
//     double *f  /* nxm matrix of KDE products */
//) {
//  int n=*nn, i, ii;
//  int j, ell;//k, kk;
//  double sum1, sum2, tmp, h=*hh, xik, u;
//  double const1 = -0.5 / (h * h);
//  double const2;
//
//  for(j=0; j<m; j++) {
//    /* In const2 expression, .3989... equals 1/(sqrt(2*pi)) */
//    const2 = 0.39894228040143267794/(h*sigma[j]*(double)(*rr)); 
//    for (i=0; i<n; i++) {
//      for (k=0; k<r; k++) {
//        u = (x[i + k*n] - mu[j])/sigma[j];
//        f[i + j*n] = 1.0;
//        sum1 = 0.0;
//        for (ii=0; ii<n; ii++) {
//          sum2 = 0.0;
//          for (kk=0; kk<r; kk++) {
//            tmp = (u - x[ii+kk*n] + mu[j])/sigma[j];
//            sum2 += exp(tmp * tmp * const1);
//          }
//          sum1 += z[ii + j*n] * sum2;
//        }
//        f[i + j*n] *= sum1 * const2;
//      }
//    }
//  }
//}
//
