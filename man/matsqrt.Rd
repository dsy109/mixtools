\name{matsqrt}
\title{Calculates the Square Root of a Diagonalizable Matrix}
\alias{matsqrt}
\usage{
matsqrt(x)
}
\description{
Returns the square root of a diagonalizable matrix.
}
\arguments{
  \item{x}{An nxn diagonalizable matrix.}
}
\value{
  \code{matsqrt} returns the square root of \code{x}.  
}
\details{
This function is called by \code{regcr}.
}
\seealso{
\code{\link{regcr}}
}
\examples{
a <- matrix(c(1, -0.2, -0.2, 1), 2, 2)
matsqrt(a)

}

\keyword{internal}
