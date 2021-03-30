## Use an ECM algorithm (in the sense of Meng and Rubin, Biometrika 1993)
## to search for a local maximum of the likelihood surface for a 
## univariate finite mixture of normals with possible equality
## constraints on the stdev parameters.

## It is assumed here that there are three components and the three normal means
## are equal to alpha, alpha-delta, and alpha+delta for unknown parameters
## alpha and delta.
## In other words, this function implements the specific model described in 
## Thomas et al (2009), Extensions of Reliability Theory.
## It is a modified version of normalmixEM.

tauequivnormalmixEM <-
function (x, lambda = NULL, mu = NULL, sigma = NULL, k = 3, 
          mean.constr = NULL, sd.constr = NULL, gparam = NULL,
          epsilon = 1e-08, maxit = 10000, maxrestarts=20, 
          verb = FALSE, fast=FALSE, ECM = TRUE,
          arbmean = TRUE, arbvar = TRUE) {
  
  M <- A <- NULL
  if (is.null(mean.constr)) {
  	  # In this case, we will be fitting a 3-component mixture model with means
  	  # constrained to be alpha, alpha-delta, and alpha+delta for 
  	  # parameters alpha and delta.
  	  k <- 3
  	  if (length(mu) != 3) mu <- NULL
  	  if (length(sigma) != 3) sigma <- NULL
  	  M <- matrix(c(1, 1, 1, 0, -1, 1), 3, 2)
  	  
  	  # We will also constain the reciprocals of the variances to be 
  	  # gamma_1+gamma_2, gamma_1, and gamma_1 for positive
  	  # parameters gamma_1 and gamma_2.
  	  A <- matrix(c(1, 1, 1, 1, 0, 0), 3, 2)
  }
  normalmixMMlc(x, 
    lambda = lambda,
    mu = mu,
    sigma = sigma,
    k = k,
    mean.constr = mean.constr,
    mean.lincstr = M,
    var.lincstr = A,
    gparam = gparam,
    epsilon = epsilon,
    maxit = maxit, 
    maxrestarts = maxrestarts,
    verb = verb)
}

