#include <R.h>
#include <Rmath.h>
#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* In the npMSL algorithm (formerly NEMS), it is necessary to store 
   the values of each 
   f_{jk} density on a univariate grid of points, since the E-step involves
   a convolution (which means that it no longer suffices to
   store f_{jk} only at the data points x_{ik}, as it did for
   the original npEM algorithm of Benaglia et al).  
   
   In the M-step, we must set each f_{jk}(u_a) equal to the sum (over i) of
   p_{ij} * phi((u_a - x_{ik}) / h) / h, where phi is the std. normal density
   which gives the normal kernel with bandwidth h.  Thus, the M-step involves
   four nested loops, one each for i, j, k, and a.
   
   In the E-step, we must set each p_{ij} equal to K * lambda_j * N(f_j)(x_i),
   where:
   -- K is a normalizing constant ensuring that sum_j p_{ij} = 1 for all i.
   -- N(f_j)(x_i) is the product prod_j N(f_jk)(x_{ik})
   -- N(f_{jk})(x_{ik}) = exp{integral of (1/h)*phi((u-x_{ik})/h)*log(f_{jk}(u))du}
   -- that last integral is approximated by the sum (over a) of
      (normalizing constant)*exp((u_a-x_{ik})^2/h^2)*log(f_{jk}(u_a))

  Current version: with underflow handling in the E-step
		   and block structure,
*/
/* **********************************************************
   Version allowing blocks of conditionaly iid coordinates: 
   f is now a ngrid by m by B array of f_{j ell} values on the grid
   note: *BB is not used here, but remains in the parameter list for
   consistency with the npMSL_Estep_bw() version which needs it
*/
void npMSL_Estep(
    int *nngrid, /* size of grid */
    int *nn, /* sample size */
    int *mm, /* number of components */
    int *rr, /* number of repeated measurements */
    int *BB, /* total nb of blocks (not used in samebw version) */
    int *blockid, /* r-vector (b_k) block to which belongs coord k */
    double *hh, /* bandwidth */
    double *data,  /* n by r vector of observations */
    double *grid, /* grid points */
    double *f, /* ngrid by m by B array of density values on grid */
    double *lambda, /* current vector of mixing parameters */
    double *post, /* n by m matrix of posterior probabilities */      
    double *loglik, /* scalar value of penalized loglikelihood */
    int *nb_udfl,	/* nb of underflows log(0) cancelled  */
    int *nb_nan		/* nb of nonzero K()*log(0) cancelled */
    ) {
  int n=*nn, m=*mm, r=*rr, ngrid=*nngrid;
  int i, j, k, ell, a;
  double sum, conv, xik, *fjl, two_h_squared =2*(*hh)*(*hh);
  double Delta = (grid[2]-grid[1]) / *hh / sqrt(2*3.14159265358979);
  double t1, expminus500=exp(-500);
  double epsi=1e-323;	/* smallest number; maybe machine-dependent ? */
  double epsi2=1e-100;	/* assumed small enough for cancelling log(0) */ 
 
  *loglik=0.0;
  for (i=0; i<n; i++) {
    sum=0.0;
    for (j=0; j<m; j++) {
      post[i + n*j] = lambda[j];
      for (k=0; k<r; k++) {
        /* select the block index of coordinate k, b_k in the paper,
	   and the appropriate f_{j b_k} array of size ngrid 
	   (b_k - 1 since ell=0 to B-1) 
	*/
	ell = blockid[k] - 1; 
	fjl = f + ngrid*(j + m*ell);
        xik=data[i + n*k];
        conv = 0.0;	
        for (a = 0; a<ngrid; a++) {  /* numerical int checking underflows */
          t1 = MAX(expminus500, exp(-(xik-grid[a])*(xik-grid[a])/two_h_squared)); /* value of kernel */

          if (fjl[a] > epsi) { /* no underflow pb */
            conv += t1 * log(fjl[a]);
          }
          else if (t1 < epsi2) { /* assume kernel cancels log(0) part */
            *nb_udfl +=1;  /* count underflow replaced by 0 */
          }
          else *nb_nan +=1; /* kernel *may* be not small enough ! */
        }
        conv *= Delta; /* Delta is the normalizing constant */
        conv = exp(conv); /* conv now = Nf_{jb_k}(xik) */
        post[i + n*j] *= conv; /* numerator = lambda_j*prod_k Nf_{jb_k}(xik) */
      }
      sum += post[i + n*j];
    }
    *loglik += log(sum);
    for(j=0; j<m; j++) {
      post[i + n*j] /= sum;
    }
  }
}





/*  
  implementing the block structure
  - the role of coordinate k is replaced by block ell  
  - the densities f_{j ell} are stored on the grid 
*/
void npMSL_Mstep(
    int *nngrid, /* size of grid */
    int *nn, /* sample size */
    int *mm, /* number of components */
    int *rr, /* number of repeated measurements */
    int *BB, /* total number of blocks */
    int *BlS, /* BB vector of blocks sizes */
    int *blockid, /* r-vector (b_k) block to which belongs coord k */
    double *hh, /* bandwidth */
    double *data,  /* n by r vector of observations */
    double *grid, /* grid points */
    double *new_f, /* ngrid by m by B array to hold new values on grid */
    double *lambda, /* vector of lambda_j 's */
    double *post /* n by m matrix of posterior probabilities */
    ) {
  int n=*nn, m=*mm, r=*rr, ngrid=*nngrid, B=*BB, i, j, k, ell, a;
  double sum, xik, *fjl, two_h_squared =2*(*hh)*(*hh);
  double pij, gridpt, expminus500=exp(-500);
  double normconst = 0.39894228040143267794/ *hh; /* 0.39...=1/(sqrt(2*pi)) */
 
  for (j=0; j<m; j++) { 	/* for each component */
    for (ell=0; ell<B; ell++) {	/* for each block */
      fjl = new_f + ngrid*(j + m*ell); /* fjl is now an array of size ngrid */
      for (a=0; a<ngrid; a++) {
        gridpt = grid[a]; /* this is the ath value of the grid */
        sum = 0.0;
	
	/* now build the vector of data for each coordinate k such that
	   k belongs to block l, ie to index ell+1
	*/
	for (k=0; k<r; k++) {
	  
	 if (blockid[k] == (ell + 1)) {  /* coordinate k is in block ell */
	   for (i=0; i<n; i++) {
	     xik = data[i + n*k];  /* n obs for coordinate k */
	     pij = post[i + n*j];
	     sum += pij * MAX(expminus500, exp(-(xik - gridpt)*(xik-gridpt)/two_h_squared)); 
           /* Numerator: Take p_{ij} 
	      times the normal density at (x_{ik}-gridpt) */
           /* denom computed from n*lambda[j]*BlS[ell] 
	     and may even be computed outside the a loop */
	   } /* for each i */
	 }
	} /* for each k belonging to block ell */
        fjl[a] = MAX(expminus500, normconst * sum / (n*lambda[j]*BlS[ell])) ; 
        /* Note that each fjl "density" should now be correctly normalized;
           normconst is the normalizing constant for a normal 
	   kernel function */           
      } /* for a in grid */
    } /* for ell */
  }
}      



/* ************************************************************** 
 adaptive bandwidth implementation versions : h replaced by h_{j ell}                  
 **************************************************************** */

void npMSL_Estep_bw(
    int *nngrid, /* size of grid */
    int *nn, /* sample size */
    int *mm, /* number of components */
    int *rr, /* number of repeated measurements */
    int *BB, /* total number of blocks */
    int *blockid, /* r-vector (b_k) block to which belongs coord k */
    double *hh, /* bandwidth matrix (B,m) */
    double *data,  /* n by r vector of observations */
    double *grid, /* grid points */
    double *f, /* ngrid by m by B array of density values on grid */
    double *lambda, /* current vector of mixing parameters */
    double *post, /* n by m matrix of posterior probabilities */      
    double *loglik, /* scalar value of penalized loglikelihood */
    int *nb_udfl,	/* nb of underflows log(0) cancelled  */
    int *nb_nan		/* nb of nonzero K()*log(0) cancelled */
    ) {
  int n=*nn, m=*mm, r=*rr, ngrid=*nngrid, B=*BB;
  int i, j, k, ell, a;
  double sum, conv, xik, *fjl, hjl, two_h_squared;
  /* two_h_squared =2*(*hh)*(*hh); */
  /* double Delta = (grid[2]-grid[1]) / *hh / sqrt(2*3.14159265358979);*/
  double gsq2pi = (grid[2]-grid[1]) / sqrt(2*3.14159265358979); 
  double t1, Delta, expminus500=exp(-500);
  double epsi=1e-323;	/* smallest number; maybe machine-dependent ? */
  double epsi2=1e-100;	/* assumed small enough for cancelling log(0) */ 
 
  *loglik=0.0;
  for (i=0; i<n; i++) {
    sum=0.0;
    for (j=0; j<m; j++) {
      post[i + n*j] = lambda[j];
      for (k=0; k<r; k++) {
        /* select the block index of coordinate k, b_k in the paper,
	   and the appropriate f_{j b_k} array of size ngrid 
	   (b_k - 1 since ell=0 to B-1) 
	*/
	ell = blockid[k] - 1; 
	fjl = f + ngrid*(j + m*ell);
        xik=data[i + n*k];
        
        hjl = hh[ell + B*j];	/* current h_{jl} bandwidth */
        two_h_squared = 2.0*hjl*hjl;
        Delta = gsq2pi / hjl;        
        conv = 0.0;	
        for (a = 0; a<ngrid; a++) {  /* numerical int checking underflows */
          t1 = MAX(expminus500, exp(-(xik - grid[a])*(xik-grid[a])/two_h_squared));

          if (fjl[a] > epsi) { /* no underflow pb */
            conv += t1 * log(fjl[a]);
          }
          else if (t1 < epsi2) { /* assume kernel cancels log(0) part */
            *nb_udfl +=1;  /* count underlow replaced by 0 */
          }
          else *nb_nan +=1; /* kernel *may* be not small enough ! */
        }
        conv *= Delta; /* Delta is the normalizing constant */
        conv = exp(conv); /* conv now = Nf_{jb_k}(xik) */
        post[i + n*j] *= conv; /* numerator = lambda_j*prod_k Nf_{jb_k}(xik) */
      }
      sum += post[i + n*j];
    }
    *loglik += log(sum);
    for(j=0; j<m; j++) {
      post[i + n*j] /= sum;
    }
  }
}





/* M-step: adaptive bandwidth version */
void npMSL_Mstep_bw(
    int *nngrid,    /* size of grid */
    int *nn, 	    /* sample size */
    int *mm, 	    /* number of components */
    int *rr, 	    /* number of repeated measurements */
    int *BB, 	    /* total number of blocks */
    int *BlS, 	    /* BB vector of blocks sizes */
    int *blockid,   /* r-vector (b_k) block to which belongs coord k */
    double *hh,	    /* bandwidth matrix B times m */
    double *data,   /* n by r vector of observations */
    double *grid,   /* grid points */
    double *new_f,  /* ngrid by m by B array to hold new values on grid */
    double *lambda, /* vector of lambda_j 's */
    double *post    /* n by m matrix of posterior probabilities */
    ) {
  int n=*nn, m=*mm, r=*rr, ngrid=*nngrid, B=*BB, i, j, k, ell, a;
  double sum, xik, *fjl; /* two_h_squared =2*(*hh)*(*hh); */
  double hjl, two_h_squared; /* adaptive bw depends on j and ell */
  double pij, gridpt, expminus500=exp(-500);
  double invsq2pi = 0.39894228040143267794; /* 0.39...=1/(sqrt(2*pi)) */
  double normconst;
  /* normconst = 0.39894228040143267794/ *hh in the samebw case */
  /* here 2*h*h and normconst are computed in j,ell loops */
 
  for (j=0; j<m; j++) { 	/* for each component */
    for (ell=0; ell<B; ell++) {	/* for each block */
      fjl = new_f + ngrid*(j + m*ell); /* fjl is now an array of size ngrid */
      hjl = hh[ell + B*j];	/* current h_{jl} bandwidth */      
      two_h_squared = 2.0*hjl*hjl;	/* using adaptive h_{j ell} */
      normconst = invsq2pi/hjl;
      
      for (a=0; a<ngrid; a++) {
        gridpt = grid[a]; /* this is the ath value of the grid */
        sum = 0.0;	
	/* now build the vector of data for each coordinate k such that
	   k belongs to block ell (ie to index ell+1)
	*/
	for (k=0; k<r; k++) {  
	 if (blockid[k] == (ell + 1)) {  /* coordinate k is in block ell */	
	   for (i=0; i<n; i++) {
	   	 xik = data[i + n*k];  /* n obs for coordinate k */
	   	 pij = post[i + n*j];
	   	 sum += pij * MAX(expminus500, exp(-(xik - gridpt)*(xik-gridpt)/two_h_squared));
           /* Numerator: Take p_{ij} 
	      times the normal density at (x_{ik}-gridpt) */
           /* denom += pij; obsolete: computed from n*lambda[j]*BlS[ell] 
	     and may even be computed outside the a loop */
	   	} /* for each i */
	  }
	} /* for each k belonging to block ell */
        fjl[a] = MAX(expminus500, normconst * sum / (n*lambda[j]*BlS[ell])); 
		
        /* Note that each fjl "density" should now be correctly normalized;
           normconst is the normalizing constant for a normal 
	   kernel function */           
      } /* for a in grid */
    } /* for ell */
  }
}      

