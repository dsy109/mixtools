\name{wkde}
\title{Weighted Univariate (Normal) Kernel Density Estimate}
\alias{wkde}
\alias{wkde.symm}
\usage{
wkde(x, u=x, w=rep(1, length(x)), bw=bw.nrd0(as.vector(x)), sym=FALSE) 
}
\description{
Evaluates a weighted kernel density estimate, using a Gaussian kernel, 
at a specified vector of points.
}
\arguments{
  \item{x}{Data}
  \item{u}{Points at which density is to be estimated}
  \item{w}{Weights (same length as \code{x})}
  \item{bw}{Bandwidth}
  \item{sym}{Logical:  Symmetrize about zero?}
}
\value{A vector of the same length as \code{u}
}
\seealso{
\code{\link{npEM}}, \code{\link{ise.npEM}}
}
\references{
   \itemize{
   \item Benaglia, T., Chauveau, D., and Hunter, D. R. (2009), An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   Journal of Computational and Graphical Statistics, 18, 505-526.
   
   \item Benaglia, T., Chauveau, D., Hunter, D. R., and Young, D. (2009),
  mixtools: An R package for analyzing finite mixture models.
	Journal of Statistical Software, 32(6):1-29.
   }
}
\examples{
# Mixture with mv gaussian model
set.seed(100)
m <- 2 # no. of components
r <- 3 # no. of repeated measures (coordinates)
lambda <- c(0.4, 0.6)
mu <- matrix(c(0, 0, 0, 4, 4, 6), m, r, byrow=TRUE) # means 
sigma <- matrix(rep(1, 6), m, r, byrow=TRUE) # stdevs
centers <- matrix(c(0, 0, 0, 4, 4, 4), 2, 3, byrow=TRUE) # initial centers for est

blockid = c(1,1,2) # block structure of coordinates
n = 100
x <- rmvnormmix(n, lambda, mu, sigma) # simulated data
a <- npEM(x, centers, blockid, eps=1e-8, verb=FALSE)

par(mfrow=c(2,2))
u <- seq(min(x), max(x), len=200)
for(j in 1:2) { 
  for(b in 1:2) {
    xx <- as.vector(x[,a$blockid==b])
    wts <- rep(a$post[,j], length.out=length(xx))
    bw <- a$bandwidth
    title <- paste("j =", j, "and b =", b)
    plot(u, wkde(xx, u, wts, bw), type="l", main=title)
  }
}

}


\keyword{file}
