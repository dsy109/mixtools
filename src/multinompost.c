#include <R.h>
#include <Rmath.h>

/* Compute the matrix of "posterior" probabilities in a finite mixture of
   multinomial distributions.  The algorithm used is fairly safe from
   a numerical perspective; it avoids over- or under-flow as long as the
   values of sigma are not too large or small.  */
void multinompost(
    int *nn, /* sample size */
    int *mm, /* number of components */
    double *loglamcd, /* n by m matrix of log(lambda * component density) values */
    double *post, /* n by m matrix of posterior probabilities */
    double *loglik /* scalar loglikelihood value (input value is a constant,
                      which is modified before return). */ 
    ) {
  int n=*nn, m=*mm, i, j, maxj;
  double sum, max, *loglamcolumnptr;

  for (loglamcolumnptr=loglamcd, i=0; i<n; loglamcolumnptr+=m, i++) {
    /* In each column of loglamcd array:
       Find the largest value, subtract it from each value and then exponentiate,
       then standardize by dividing by the sum.  Add this maximum to the 
       loglikelihood along with the log of the sum. */
    max = loglamcolumnptr[maxj=0];
    for (j=1; j<m; j++) {
      if (loglamcolumnptr[j] > max) {
        max = loglamcolumnptr[j];
        maxj = j;
      }
    }
    sum = 1.0;
    for(j = 0; j<m; j++) {
      if (j != maxj) {
        sum += (post[j*n+i] = exp(loglamcolumnptr[j]-max));
      }
    }
    *loglik += log(sum) + max;
    for(j=0; j<m; j++) {
      if (j==maxj) {
        post[j*n+i] = 1.0/sum;
      } else {
        post[j*n+i] /= sum;
      }
    }
  }
}


