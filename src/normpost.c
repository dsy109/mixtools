#include <R.h>
#include <Rmath.h>

/* Compute the matrix of "posterior" probabilities in a finite mixture of
   univariate normal densities.  The algorithm used is fairly safe from
   a numerical perspective; it avoids over- or under-flow as long as the
   values of sigma are not too large or small.  */
void normpost(
    int *nn, /* sample size */
    int *mm, /* number of components */
    double *data,  /* vector of observations */
    double *mu, /* current vector of component means */
    double *sigma, /* current vector of component stdevs */
    double *lambda, /* current vector of mixing parameters */
    double *res2, /* n by m matrix of squared residuals */
    double *work, /* 3*m-vector of workspace, which will be broken into 3 parts */
    double *post, /* n by m matrix of posterior probabilities */      
    double *loglik /* scalar loglikelihood value */ 
    ) {
  int n=*nn, m=*mm, i, j, minj=0;
  double x, r, rowsum, min=0.0;
  double *LamSigRatio = work+m; /* Second 1/3 of workspace, for frequently used constants */
  double *logLamSigRatio = work+2*m; /* Third 1/3 of workspace, for frequently used constants */

  *loglik = -(double)n * 0.91893853320467274178; /* n/2 times log(2pi) */
  for (j=0; j<m; j++){ /* store some often-used values to save time later */
    LamSigRatio[j] = lambda[j] / sigma[j];
    logLamSigRatio[j] = log(LamSigRatio[j]);
  }
  for (i=0; i<n; i++){
    x=data[i];
    for (j=0; j<m; j++) {
      r = x-mu[j];
      r = r*r;
      res2[i + n*j] = r;
      work[j] = r = r / (2.0 * sigma[j] * sigma[j]);
      /* Keep track of the smallest standardized squared residual.  
         By dividing everything by the component density with the
         smallest such residual, the denominator of the posterior
         is guaranteed to be at least one and cannot be infinite unless
         the values of lambda or sigma are very large or small.  This helps
         prevent numerical problems when calculating the posteriors.*/
      if (j==0 || r < min) {        
        minj = j;
        min = r;
      }
    }
    /* At this stage, work contains the squared st'dized resids over 2 */
    rowsum = 1.0;
    for (j=0; j<m; j++) {
      if (j==minj) 
        work[j] = 1.0;
      else {
        work[j] = (LamSigRatio[j] / LamSigRatio[minj]) * exp(min - work[j]);
        rowsum += work[j];
      }
    }
    /* At this stage, work contains the normal density at data[i]
       divided by the normal density with the largest st'dized resid 
       Thus, dividing by rowsum gives the posteriors: */
    for (j=0; j<m; j++) {
      post[i + n*j] = work[j] / rowsum;
    }
    /* Finally, adjust the loglikelihood correctly */
    *loglik += log(rowsum) - min + logLamSigRatio[minj];
  }
}

void oldnormpost(
    int *nn, /* sample size */
    int *mm, /* number of components */
    double *data,  /* vector of observations */
    double *mu, /* current vector of component means */
    double *sigma, /* current vector of component stdevs */
    double *lambda, /* current vector of mixing parameters */
    double *res2, /* n by m matrix of squared residuals */
    double *work, /* m-vector of workspace */
    double *post, /* n by m matrix of posterior probabilities */      
    double *loglik /* scalar loglikelihood value */ 
    ) {
  int n=*nn, m=*mm, i, j, minj=0;
  double x, r, rowsum, min=0.0;

  *loglik = -(double)n * 0.91893853320467274178; /* n/2 times log(2pi) */
  for (i=0; i<n; i++){
    x=data[i];
    min = 1.0e+6; /* large number */
    for (j=0; j<m; j++) {
      r = x-mu[j];
      r = r*r;
      res2[i + n*j] = r;
      work[j] = r = r / (2.0 * sigma[j] * sigma[j]);
      /* Keep track of the smallest standardized squared residual.  
         By dividing everything by the component density with the
         smallest such residual, the denominator of the posterior
         is guaranteed to be at least one and cannot be infinite unless
         the values of lambda or sigma are very large or small.  This helps
         prevent numerical problems when calculating the posteriors.*/
      if (r < min) {
        minj = j;
        min = r;
      }
    }
    /* At this stage, work contains the squared st'dized resids over 2 */
    rowsum = 1.0;
    for (j=0; j<m; j++) {
      if (j==minj) 
        work[j] = 1.0;
      else 
        rowsum+=(work[j]=lambda[j]/sigma[j]*sigma[minj]/lambda[minj]*exp(min-work[j]));
    }
    /* At this stage, work contains the normal density at data[i]
       divided by the normal density with the largest st'dized resid 
       Thus, dividing by rowsum gives the posteriors: */
    for (j=0; j<m; j++) {
      post[i + n*j] = work[j] / rowsum;
    }
    /* Finally, adjust the loglikelihood correctly */
    *loglik += log(rowsum) - min + log(lambda[minj]/sigma[minj]);
  }
}

