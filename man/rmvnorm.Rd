\name{rmvnorm}
\alias{rmvnorm}
\title{Simulate from a Multivariate Normal Distribution}
\description{
  Simulate from a multiviate normal distribution
}
\usage{
rmvnorm(n, mu=NULL, sigma=NULL)
}
\arguments{
  \item{n}{Number of vectors to simulate}
  \item{mu}{mean vector}
  \item{sigma}{covariance matrix, assumed symmetric and nonnegative definite}
}
\value{
  An \eqn{n \times d}{n x d} matrix in which each row is an independently
  generated realization from the desired multivariate normal distribution
}
\details{
This function uses an \code{\link{eigen}} decomposition assuming \code{sigma} is symmetric.
In particular, the upper triangle of \code{sigma} is ignored.
}
\seealso{
  \code{\link{eigen}}, \code{\link{dnorm}}, \code{\link{dmvnorm}}
}
\keyword{distribution}
