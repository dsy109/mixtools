#include <math.h> /* in order to make the sqrt function work */

/* First, create the function prototype (global declaration that there
exists a function called new_svalues returning void and taking these 
arguments) */
void new_svalues (double *z, double *y, double *x, double *beta,
		  int *k, int *n, int *p, double *out, double *sz, double *runsum);


/* Next comes the function itself:  */
void new_svalues (double *z, double *y, double *x, double *beta,
		  int *k, int *n, int *p, double *out, double *sz, double *runsum) {
  int i, j, l;
  double sum;
  double xbeta;
  double diff;
  double diff2;
  double zdiff;
  double zdiff2;
  double zdiff3;
  
  /*  Create the sz (column sum of z) vector */
  for(j=0; j < *k ; j=j+1) {
    sum=0.0;
    for(i=0; i < *n; i=i+1) {
      sum += z[j * (*n) + i];
    }
    sz[j]=sum;
  }

  for(j=0; j < *k ; j=j+1) {
	zdiff=0.0;
    for(i=0; i < *n; i=i+1) {
	xbeta=0.0;
      /* Calculate i,j component of x %*% beta */
	for(l=0; l < *p; l++){
	xbeta += x[l * (*n) + i] * beta[j * (*p) + l];
	}
      /* Subtract it from the i component of y */
	diff = y[i] - xbeta;
      /* square the difference */
	diff2 = pow(diff,2);
      /* multiply by the i,j component of z */
	zdiff += z[j * (*n) + i] * diff2;
    }
      /* keep track of running sum */
	runsum[j] = zdiff;
    /* divide dot product by sz[j] */
	zdiff2 = runsum[j] / sz[j];
    /* take the square root */
	zdiff3 = sqrt(zdiff2);
    /* put result in out[j] */
	out[j] = zdiff3;
  }
  
}
