\name{flaremixEM}
\title{EM Algorithm for Mixtures of Regressions with Flare}
\alias{flaremixEM}
\usage{
flaremixEM(y, x, lambda = NULL, beta = NULL, sigma = NULL, 
           alpha = NULL, nu = NULL, epsilon = 1e-04, 
           maxit = 10000, verb = FALSE, restart = 50)
}

\description{
  Returns output for 2-component mixture of regressions with flaring
  using an EM algorithm with one step of Newton-Raphson requiring
  an adaptive barrier for maximization of the objective function. A mixture of regressions
  with flare occurs when there appears to be a common regression relationship for the data,
  but the error terms have a mixture structure of one normal component and one exponential
  component.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An n-vector of predictor values.  An intercept term will be added by default.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to 1.}
  \item{beta}{Initial value of \code{beta} parameters.  Should be a 2x2 matrix where the columns
    correspond to the component.}
  \item{sigma}{A vector of standard deviations.}
  \item{alpha}{A scalar for the exponential component's rate.}
  \item{nu}{A vector specifying the barrier constants to use.  The first barrier constant where the algorithm
  converges is used.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
  \item{restart}{The number of times to restart the algorithm in case convergence is not attained.
  The default is 50.}
}
\value{
  \code{flaremixEM} returns a list of class \code{mixEM} with items:
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
\code{\link{regmixEM}}
}

\examples{
## Simulation output.

set.seed(100)
j=1
while(j == 1){
    x1 <- runif(30, 0, 10)
    x2 <- runif(20, 10, 20)
    x3 <- runif(30, 20, 30)
    y1 <- 3+4*x1+rnorm(30, sd = 1)
    y2 <- 3+4*x2+rexp(20, rate = .05)
    y3 <- 3+4*x3+rnorm(30, sd = 1)
    x <- c(x1, x2, x3)
    y <- c(y1, y2, y3)
    nu <- (1:30)/2

    out <- try(flaremixEM(y, x, beta = c(3, 4), nu = nu,
               lambda = c(.75, .25), sigma = 1), silent = TRUE)
    if(any(class(out) == "try-error")){
        j <- 1
    } else j <- 2
}

out[4:7]
plot(x, y, pch = 19)
abline(out$beta)



}


\keyword{file}
