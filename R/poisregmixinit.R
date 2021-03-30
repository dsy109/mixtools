poisregmix.init <- function(y, x, lambda = NULL, beta = NULL, k = 2){
  x <- as.matrix(x)
  n <- length(y)
  if (is.null(lambda)) {
    cond = TRUE
    while(cond){
      lambda = runif(k)
      lambda = lambda/sum(lambda)
      if(min(lambda)<0.05) cond=TRUE else cond=FALSE
    }
  } else k = length(lambda)
  p <- ncol(x)
  w = cbind(y, x)
  w = w[order(w[, 1]), ]
  w.bin = list()
  for (j in 1:k) {
    w.bin[[j]] <- w[max(1, floor((j - 1) * n/k)):ceiling(j * 
                                                           n/k), ]
  }
  if (is.null(beta)) {
    beta.hyp = matrix(sapply(1:k, function(j) glm.fit(w.bin[[j]][, 
                                                                 2:(p + 1)], w.bin[[j]][, 1], family = poisson())$coeff), 
                      ncol = k)
    sd.hyp = apply(beta.hyp, 1, sd)
    beta = matrix(0, p, k)
    for (j in 1:k) {
      beta[, j] = rnorm(p, mean = as.vector(beta.hyp[, 
                                                     j]), sd = sd.hyp)
    }
  }
  else k = ncol(beta)
  list(lambda = lambda, beta = beta, k = k)
}
