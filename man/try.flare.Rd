\name{try.flare}
\title{Mixtures of Regressions with Flare MM Algorithm}
\alias{try.flare}
\usage{
try.flare(y, x, lambda = NULL, beta = NULL, sigma = NULL, 
          alpha = NULL, nu = 1, epsilon = 1e-04, 
          maxit = 10000, verb = FALSE, restart = 50)
}

\description{
  The function which \code{flaremixEM} actually calls.  This only allows
  one barrier constant to be inputted at a time.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An n-vector of predictor values.  An intercept term will be added by default.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to 1.}
  \item{beta}{Initial value of \code{beta} parameters.  Should be a 2x2 matrix where the columns
    corresond to the component.}
  \item{sigma}{A vector of standard deviations.}
  \item{alpha}{A scalar for the exponential component's rate.}
  \item{nu}{A scalar specifying the barrier constant to use.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
  \item{restart}{The number of times to restart the algorithm in case convergence is not attained.
  The default is 50.}
}
\value{
  \code{try.flare} returns a list of class \code{mixEM} with items:
  \item{x}{The set of predictors (which includes a column of 1's).}
  \item{y}{The response values.}
  \item{posterior}{An nx2 matrix of posterior probabilities for
   observations.}
  \item{lambda}{The final mixing proportions.}
  \item{beta}{The final regression coefficients.}
  \item{sigma}{The final standard deviations.}
  \item{alpha}{The final exponential rate.}
  \item{loglik}{The final log-likelihood.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{flaremixEM}}
}

\details{
  This usually is not called by the user.  The user will likely want \code{flaremixEM}, which also
  has an example to demonstrate this algorithm.
}
\keyword{internal}
