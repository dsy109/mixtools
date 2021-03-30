\name{ldmult}
\title{Log-Density for Multinomial Distribution}
\alias{ldmult}
\usage{
ldmult(y, theta)
}
\description{
Return the logarithm of the multinomial density function.
}
\arguments{
  \item{y}{A vector of multinomial counts.}
  \item{theta}{A vector of multinomial probabilities.  May have same number of
    components as or one fewer component than \code{y}.  In the latter case, 
    an extra component is appended so that theta sums to one.}
}
\value{
  \code{ldmult} returns the logarithm of the multinomial density
  with parameter \code{theta}, evaluated at \code{y}.  
}
\details{
This function is called by \code{multmixEM}.
}
\seealso{
\code{\link{multmixEM}}
}
\examples{
y <- c(2, 2, 10)
theta <- c(0.2, 0.3, 0.5)
ldmult(y, theta)

}

\keyword{internal}
