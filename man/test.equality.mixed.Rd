\name{test.equality.mixed}
\title{Performs Chi-Square Test for Mixed Effects Mixtures}
\alias{test.equality.mixed}
\usage{
test.equality.mixed(y, x, w=NULL, arb.R = TRUE, 
                    arb.sigma = FALSE, lambda = NULL, 
                    mu = NULL, sigma = NULL, R = NULL, 
                    alpha = NULL, ...)
}

\description{
  Performs a likelihood ratio test of either common variance terms between the response trajectories in a mixture
  of random (or mixed) effects regressions or for common variance-covariance matrices for the random effects mixture distribution.}
\arguments{
  \item{y}{The responses for \code{regmixEM.mixed}.}
  \item{x}{The predictors for the random effects in \code{regmixEM.mixed}.}
  \item{w}{The predictors for the (optional) fixed effects in \code{regmixEM.mixed}.}
  \item{arb.R}{If FALSE, then a test for different variance-covariance matrices for the random effects mixture is performed.}
  \item{arb.sigma}{If FALSE, then a test for different variance terms between the response trajectories is performed.}
  \item{lambda}{A vector of mixing proportions (under the null hypothesis) with same purpose as outlined in \code{regmixEM.mixed}.}
  \item{mu}{A matrix of the means (under the null hypothesis) with same purpose as outlined in \code{regmixEM.mixed}.}
  \item{sigma}{A vector of standard deviations (under the null hypothesis) with same purpose as outlined in \code{regmixEM.mixed}.}
  \item{R}{A list of covariance matrices  (under the null hypothesis) with same purpose as outlined in \code{regmixEM.mixed}.}
  \item{alpha}{An optional vector of fixed effects regression coefficients  (under the null hypothesis) with same purpose as outlined
  in \code{regmixEM.mixed}.}
  \item{...}{Additional arguments passed to \code{regmixEM.mixed}.} 
}
\value{
  \code{test.equality.mixed} returns a list with the following items:
  \item{chi.sq}{The chi-squared test statistic.}
  \item{df}{The degrees of freedom for the chi-squared test statistic.}
  \item{p.value}{The p-value corresponding to this likelihood ratio test.}
}
\seealso{
\code{\link{test.equality}}
}
\examples{
##Test of equal variances in the simulated data set.

data(RanEffdata)
set.seed(100)
x<-lapply(1:length(RanEffdata), function(i) 
          matrix(RanEffdata[[i]][, 2:3], ncol = 2))
x<-x[1:15]
y<-lapply(1:length(RanEffdata), function(i) 
          matrix(RanEffdata[[i]][, 1], ncol = 1))
y<-y[1:15]

out<-test.equality.mixed(y, x, arb.R = TRUE, arb.sigma = FALSE,
                         epsilon = 1e-1,  verb = TRUE,
                         maxit = 50,
                         addintercept.random = FALSE)
out
}

\keyword{file}
