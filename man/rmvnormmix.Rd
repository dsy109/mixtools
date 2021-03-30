\name{rmvnormmix}
\title{Simulate from Multivariate (repeated measures) Mixtures of Normals}
\alias{normmixrm.sim}
\alias{rmvnormmix}
\usage{
rmvnormmix(n, lambda=1, mu=0, sigma=1)
}

\description{
Simulate from a mixture of multivariate zero-correlation normal distributions
}
\arguments{
  \item{n}{Number of cases to simulate.}
  \item{lambda}{Vector of mixture probabilities with length equal to \eqn{m}, the
    desired number of components.  This is assumed to sum
    to 1; if not, it is normalized.}
  \item{mu}{Matrix of means of dimensions \eqn{m\times r}{m x r}, where
  \eqn{m} is the number of components (subpopulations) and \eqn{r} is
  the number of coordinates (repeated measurements) per case.  Note:  \code{mu}
  is automatically coerced to a matrix with \eqn{m} rows even if it is not
  given in this form, which can lead to unexpected behavior in some cases.}
  \item{sigma}{Matrix of standard deviations, same dimensions as \code{mu}.
  The coordinates within a case are independent, conditional on the mixture
  component. (There is marginal correlation among the coordinates, but
  this is due to the mixture structure only.)  Note:  \code{sigma} is
  automatically coerced to a matrix with \eqn{m} rows even if it is not
  given in this form, which can lead to unexpected behavior in some cases.}
}
\details{
It is possible to generate univariate standard normal random variables using
the default values (but why bother?).  The case of conditionally iid
coordinates is covered by the situation in which
all columns in mu and sigma are identical.
}
\value{
  \code{rmvnormmix} returns an \eqn{n\times r}{n x r} matrix in which each row is
  a sample from one of the components of a mixture of zero-correlation 
  multivariate normals.  The mixture structure
  induces nonzero correlations among the coordinates.
}
 
\seealso{
\code{\link{rnormmix}}
}
\examples{
##Generate data from a 2-component mixture of trivariate normals.

set.seed(100)
n <- 200
lambda <- rep(1, 2)/2
mu <- matrix(2*(1:6), 2, 3)
sigma <- matrix(1,2,3)
mydata<-rmvnormmix(n, lambda, mu, sigma)

## Now check to see if we can estimate mixture densities well:
title <- paste("Does this resemble N(", mu[1,], ",1) and N(", mu[2,],",1)?",
                sep="")
plot(npEM(mydata, 2), title=title)
}

\keyword{file}
