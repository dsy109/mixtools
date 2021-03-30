\name{rweibullmix}
\title{Simulate from Mixtures of Weibull distributions}
\alias{rweibullmix}
\usage{
rweibullmix(n, lambda = 1, shape = 1, scale = 1)
}

\description{
Simulate from a mixture of univariate Weibull distributions.
}
\arguments{
  \item{n}{Number of cases to simulate.}
  \item{lambda}{Vector of mixture probabilities, with length equal to \eqn{m}, 
    the desired number of components (subpopulations).  This is assumed to sum
    to 1.}
  \item{shape}{Vector of component shapes.}
  \item{scale}{Vector of component scales.}
}
\value{
  \code{rexpmix} returns an \eqn{n}-vector sampled from an \eqn{m}-component 
  mixture of univariate Weibull distributions.
}
%%\details{This function simply calls \code{\link{sample}} and \code{\link{rexp}}.} 
\seealso{
\code{\link{rnormmix}} and \code{\link{rmvnormmix}} for Gaussian mixtures,
\code{\link{rexpmix}} for mixture of exponentials.
}
\examples{
n = 500 # sample size
m = 2 # nb components
lambda=c(0.4, 0.6)
shape <- c(0.5,5); scale <- c(1,20) # model parameters
set.seed(321)
x <- rweibullmix(n, lambda, shape, scale) # iid ~ weibull mixture

## histogram of the simulated data.
hist(x, col=8)
}

\keyword{file}
