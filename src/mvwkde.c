#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <stdio.h>

/*  Simultaneously calculate m different multivariate weighted KDEs (mvwkde)
    1 for each component, for a fixed value of block index l,
    and same bandwidth for all components (option samebw=TRUE)
    final value is a (n,m) matrix f where f[i,j]=\hat f_{jl}(u_i)  */

void mvwkde_samebw (
     int *nn, /* Sample size */
     int *dd, /* Number of coordinates in l-th block */
     int *mm, /* Number of components */
     double *h, /* bandwidth d-vector */
     double *x, /* data:  vector of length nn*rr */
     double *u, /* data:  vector of length nn*rr (u=xfor now, ToDo for any u!) */
     double *z, /* nn*mm vector of normalized posteriors */
     double *f  /* nxm matrix of weighted KDE - multidim */
) {
  
    int n=*nn, d=*dd, dn=*dd*n, mn=*mm*n,  iu, ix;
    int kn, jn, id;
    double tmp1, tmp2, sum1, sum2, xik, uik, hk;
    double c2, det_h;
  /*const float const1 = -0.5;*/
    double const2 = 0.39894228040143267794; /* =1/(sqrt(2*pi)) */

    /* computing the constant from bandwidth matrix and d */
    det_h = 1.0;
    for (id=0; id<d; id++) det_h *= h[id];
    c2 = exp((double)d *log(const2)) / det_h;

  for(jn=0; jn<mn; jn+=n) { /* jn is component index times n */
    for(iu=0; iu<n; iu++) { 
	sum2 = 0.0;
      	for (ix=0; ix<n; ix++){
		sum1 = 0.0;
		for(kn=0; kn<dn; kn+=n) { /* kn is coordinate index within block times n */
        	   uik = u[iu + kn]; /* notation uik is confusing, same as xik */
	           xik = x[ix + kn]; 
		   hk = h[kn/n];
		   tmp1 = ((uik-xik)/hk)*((uik-xik)/hk);
		   sum1 += tmp1;
		}
        tmp2 = z[ix + jn]*exp(sum1*(-0.5));  /* does not depends on kn */
		sum2 += tmp2;
      	}
	f[iu + jn] = c2*sum2;
    }
  }
}



void mvwkde_adaptbw (
     int *nn, /* Sample size */
     int *dd, /* Number of coordinates in l-th block */
     int *mm, /* Number of components */
     double *H, /* bandwidth m*d-matrix */
     double *x, /* data:  vector of length nn*rr */
     double *u, /* data:  vector of length nn*rr (u=xfor now, ToDo for any u!) */
     double *z, /* nn*mm vector of normalized posteriors */

     double *f  /* nxm matrix of weighted KDE - multidim */
) {
  
    int n=*nn, d=*dd, dn=*dd*n, mn=*mm*n, m=*mm, iu, ix, jh, dm=*dd*m;
    int kn, jn, id;
    double tmp1, tmp2, sum1, sum2, xik, uik, hk, tmp0;
    double c2, det_h, h[100];
  /*const float const1 = -0.5;*/
    double const2 = 0.39894228040143267794; /* =1/(sqrt(2*pi)) */ 

   for(jn=0; jn<mn; jn+=n) { /* jn is component index times n */ 
   	for (jh=0; jh<dm; jh+=m){
		tmp0 = H[jn/n + jh];
		h[jh/m] = tmp0;
	}
	det_h = 1.0;
    	
	for (id=0; id<d; id++) det_h *= h[id];
   	c2 = exp((double)d *log(const2)) / det_h;

    	for(iu=0; iu<n; iu++) { 
		sum2 = 0.0;
      		for (ix=0; ix<n; ix++){
			sum1 = 0.0;
			for(kn=0; kn<dn; kn+=n) { /* kn is coordinate index within block times n */
        	  	   uik = u[iu + kn]; /* notation uik is confusing, same as xik */
	           	   xik = x[ix + kn]; 
		 	   hk = h[kn/n];
		   	   tmp1 = ((uik-xik)/hk)*((uik-xik)/hk);
		   	   sum1 += tmp1;
			}
       		tmp2 = z[ix + jn]*exp(sum1*(-0.5));  /* does not depends on kn */
		sum2 += tmp2;
      		}
	f[iu + jn] = c2*sum2;
    	}
   }
}




