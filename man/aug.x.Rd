\name{aug.x}
\title{Augmented Predictor Function}
\alias{aug.x}
\usage{
aug.x(X, cp.locs, cp, delta = NULL)
}
\description{
Creates the augmented predictor matrix based on an appropriately defined changepoint structure.
}
\arguments{
  \item{X}{The raw matrix of predictor values.  Note that the raw data matrix should not include a columns of 1's.}
  \item{cp.locs}{The locations of the changepoints.  The length of this vector must be equal to the sum of the entries of \code{cp}.}
  \item{cp}{A vector having length equal to the number of predictors.}
  \item{delta}{A vector to accommodate discontinuities.  If NULL, then no discontinuities are included.  Otherwise, this must be a vector of the same length as \code{cp.locs}.}
}
\value{
  \code{aug.x} returns a matrix of the original matrix \code{X} with the predictor adjusted for changepoints and (optional) discontinuities.
}
\details{
This function is called by \code{segregmixEM} and the associated internal functions.
}
\seealso{
\code{\link{segregmixEM}}
}
\examples{
x <- matrix(1:30, nrow = 10)
cp <- c(1, 3, 0)
cp.locs <- c(3, 12, 14, 16)
d <- rep(0, 4)
x1 <- aug.x(x, cp.locs, cp, delta = NULL)
x1
x2 <- aug.x(x, cp.locs, cp, delta = d)
x2
}


\keyword{internal}
