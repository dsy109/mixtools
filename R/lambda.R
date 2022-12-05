lambda <- function (z, x, xi, h = NULL, kernel = c("Gaussian", "Beta", 
                                                 "Triangle", "Cosinus", "Optcosinus"), g = 0) 
{
  if (is.null(h)) {
    cat("WARNING! BANDWIDTH MUST BE SPECIFIED!", "\n")
  }
  n <- nrow(xi)
  p <- ncol(xi)
  if(length(h)==1) h <- rep(h,p)
  if(length(h)!=p) {
    stop("Check the length of the bandwidth h.")
  } else{
    h <- matrix(rep(h,each=n),nrow=n)
  }
  x <- matrix(rep(x,each=n),nrow=n)
  X.mat=cbind(1,(xi-x))
  kernel <- match.arg(kernel)
  inwindow <- (abs((xi - x)/h) <= 1)
  if (kernel == "Gaussian") {
    W=(kern.G(x, xi, h) * inwindow)
  }
  else if (kernel == "Beta") {
    W=(kern.B(x, xi, h, g) * inwindow)
  }
  else if (kernel == "Triangle") {
    W=(kern.T(x, xi, h) * inwindow)
  }
  else if (kernel == "Cosinus") {
    W=(kern.C(x, xi, h) * inwindow)
  }
  else if (kernel == "Optcosinus") {
    W=(kern.O(x, xi, h) * inwindow)
  }
  W = diag(apply(matrix(W,ncol=ncol(X.mat)-1),1,prod))
  A=try(solve(t(X.mat)%*%(W%*%X.mat)), silent=TRUE)
  if(inherits(A, "try-error", which = TRUE)) {
    A=ginv(t(X.mat)%*%(W%*%X.mat))
  }
  beta.x=A%*%t(X.mat)%*%(W%*%cbind(z))
  beta.x
}