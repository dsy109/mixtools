\name{rexpmix}
\title{Simulate from Mixtures of Exponentials}
\alias{rexpmix}
\usage{
rexpmix(n, lambda = 1, rate = 1)
}

\description{
Simulate from a mixture of univariate exponential distributions.
}
\arguments{
  \item{n}{Number of cases to simulate.}
  \item{lambda}{Vector of mixture probabilities, with length equal to \eqn{m}, 
    the desired number of components (subpopulations).  This is assumed to sum
    to 1.}
  \item{rate}{Vector of component rates.}
}
\value{
  \code{rexpmix} returns an \eqn{n}-vector sampled from an \eqn{m}-component 
  mixture of univariate exponential distributions.
}
%%\details{This function simply calls \code{\link{sample}} and \code{\link{rexp}}.} 
\seealso{
\code{\link{rnormmix}}, \code{\link{rmvnormmix}} for Gaussian mixtures,
\code{\link{rweibullmix}} for mixture of Weibull distributions.
}
\examples{
## Generate data from a 2-component mixture of exponentials.
n=300 # sample size
m=2 # nb components
lambda=c(1/3, 2/3); rate = c(1,1/10) # parameters
set.seed(1234)
x <- rexpmix(n, lambda, rate) # iid ~ exp mixture

## histogram of the simulated data.
hist(x, col=8)
}

\keyword{file}
