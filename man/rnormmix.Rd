\name{rnormmix}
\title{Simulate from Mixtures of Normals}
\alias{normmix.sim}
\alias{rnormmix}
\usage{
rnormmix(n, lambda=1, mu=0, sigma=1)
}

\description{
Simulate from a mixture of univariate normal distributions.
}
\arguments{
  \item{n}{Number of cases to simulate.}
  \item{lambda}{Vector of mixture probabilities, with length equal to \eqn{m},
    the desired number of components (subpopulations).  This is assumed to sum
    to 1; if not, it is normalized.}
  \item{mu}{Vector of means.}
  \item{sigma}{Vector of standard deviations.}
}
\value{
  \code{rnormmix} returns an \eqn{n}-vector sampled from an \eqn{m}-component 
  mixture of univariate normal distributions.
}
\details{This function simply calls \code{\link{rmvnormmix}}.} 
\seealso{
\code{\link{makemultdata}}, \code{\link{rmvnormmix}}
}
\examples{
##Generate data from a 2-component mixture of normals.

set.seed(100)
n <- 500
lambda <- rep(1, 2)/2
mu <- c(0, 5)
sigma <- rep(1, 2)
mixnorm.data <- rnormmix(n, lambda, mu, sigma)

##A histogram of the simulated data.

hist(mixnorm.data)
}

\keyword{file}
