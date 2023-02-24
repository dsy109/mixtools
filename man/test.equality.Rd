\name{test.equality}
\title{Performs Chi-Square Tests for Scale and Location Mixtures}
\alias{test.equality}
\usage{
test.equality(y, x = NULL, arbmean = TRUE, arbvar = FALSE, 
              mu = NULL, sigma = NULL, beta = NULL, 
              lambda = NULL, ...)
}

\description{
  Performs a likelihood ratio test of a location (or scale) normal or regression mixture versus the more
  general model.  For a normal mixture, the alternative hypothesis is that each
  component has its own mean and variance, whereas the null is that all means
  (in the case of a scale mixture) or all variances (in the case of a location
  mixture) are equal.
  This test is asymptotically chi-square with degrees of freedom equal to k-1,
  where k is
  the number of components.
}
\arguments{
  \item{y}{The responses for \code{regmixEM} or the data for \code{normalmixEM}.}
  \item{x}{The predictors for \code{regmixEM}.}
  \item{arbmean}{If FALSE, then a scale mixture analysis is performed for \code{normalmixEM} or \code{regmixEM}.}
  \item{arbvar}{If FALSE, then a location mixture analysis is performed for \code{normalmixEM} or \code{regmixEM}.}
  \item{mu}{An optional vector for starting values (under the null hypothesis) for \code{mu} in \code{normalmixEM}.}
  \item{sigma}{An optional vector for starting values (under the null hypothesis) for \code{sigma} in \code{normalmixEM} or \code{regmixEM}.}
  \item{beta}{An optional matrix for starting values (under the null hypothesis) for \code{beta} in \code{regmixEM}.}
  \item{lambda}{An optional vector for starting values (under the null hypothesis) for \code{lambda} in \code{normalmixEM} or \code{regmixEM}.}
  \item{...}{Additional arguments passed to the various EM algorithms for the mixture of interest.} 
}
\value{
  \code{test.equality} returns a list with the following items:
  \item{chi.sq}{The chi-squared test statistic.}
  \item{df}{The degrees of freedom for the chi-squared test statistic.}
  \item{p.value}{The p-value corresponding to this likelihood ratio test.}
}
\seealso{
\code{\link{test.equality.mixed}}
}
\examples{
## Should a location mixture be used for the Old Faithful data?

data(faithful)
attach(faithful)
set.seed(100)
test.equality(y = waiting, arbmean = FALSE, arbvar = TRUE)

}


\keyword{file}
