\name{lambda}
\title{Local Estimation for Lambda in Mixtures of Regressions}
\alias{lambda}
\usage{
lambda(z, x, xi, h = NULL, kernel = c("Gaussian", "Beta", 
       "Triangle", "Cosinus", "Optcosinus"), g = 0)
}

\description{
  Return local estimates of the mixing proportions from each component
  of a mixture of regressions model using output from an EM algorithm.   
}
\arguments{
  \item{z}{An nxk matrix of posterior probabilities obtained from the EM
    algorithm.}
  \item{x}{A vector of values for which the local estimation is calculated.}
  \item{xi}{An nx(p-1) matrix of the predictor values.}
  \item{h}{The bandwidth controlling the size of the window used for the
  local estimation.}
  \item{kernel}{The type of kernel to be used for the local estimation.}
  \item{g}{A shape parameter required for the symmetric beta kernel.  The default
  is \code{g} = 0 which yields the uniform kernel.  Some common values are \code{g} = 1 for the
  Epanechnikov kernel, \code{g} = 2 for the biweight kernel, and \code{g} = 3 for the triweight kernel.}
}
\value{
  \code{lambda} returns local estimates of the mixing proportions for the inputted
  \code{x} vector.
}
\seealso{
\code{\link{regmixEM.loc}}
}
\note{\code{lambda} is for use within \code{regmixEM.loc}.}

\keyword{internal}
