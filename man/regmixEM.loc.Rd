\name{regmixEM.loc}
\title{Iterative Algorithm Using EM Algorithm for Mixtures of Regressions with 
    Local Lambda Estimates}
\alias{regmixEM.loc}
\usage{
regmixEM.loc(y, x, lambda = NULL, beta = NULL, sigma = NULL, 
             k = 2, addintercept = TRUE, kern.l = c("Gaussian",
             "Beta", "Triangle", "Cosinus", "Optcosinus"), 
             epsilon = 1e-08, maxit = 10000, kernl.g = 0, 
             kernl.h = 1, verb = FALSE) 
}

\description{
  Iterative algorithm returning EM algorithm output for mixtures of multiple regressions where the mixing proportions
  are estimated locally.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An nxp matrix of predictors.  See \code{addintercept} below.}
  \item{lambda}{An nxk matrix of initial local values of mixing proportions.  
    Entries should sum to 1.  This determines number of components.  
    If NULL, then \code{lambda} is simply one over the number of components.}
  \item{beta}{Initial global values of \code{beta} parameters.  Should be a pxk matrix,
    where p is the number of columns of x and \code{k} is number of components.
    If NULL, then \code{beta} has uniform standard normal entries.  If both
    \code{lambda} and \code{beta} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{A k-vector of initial global values of standard deviations.
    If NULL, then \eqn{1/\code{sigma}^2} has random standard exponential entries.  
    If \code{lambda}, \code{beta}, and \code{sigma} are NULL, then number of components determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless all of \code{lambda}, \code{beta},
    and \code{sigma} are NULL.}
  \item{addintercept}{If TRUE, a column of ones is appended to the x
    matrix before the value of p is calculated.}
  \item{kern.l}{The type of kernel to use in the local estimation of \code{lambda}.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{kernl.g}{A shape parameter required for the symmetric beta kernel for local estimation of \code{lambda}.  
  The default is g = 0 which yields the uniform kernel.  Some common values are g = 1 for the
  Epanechnikov kernel, g = 2 for the biweight kernel, and g = 3 for the triweight kernel.}
  \item{kernl.h}{The bandwidth controlling the size of the window used in the
  local estimation of lambda around x.}
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{regmixEM.loc} returns a list of class \code{mixEM} with items:
  \item{x}{The set of predictors (which includes a column of 1's if \code{addintercept} = TRUE).}
  \item{y}{The response values.}
  \item{lambda.x}{The final local mixing proportions.}
  \item{beta}{The final global regression coefficients.}
  \item{sigma}{The final global standard deviations.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{regmixEM.lambda}}
}
\examples{
## Compare a 2-component and 3-component fit to NOdata.

data(NOdata)
attach(NOdata)
set.seed(100)
out1 <- regmixEM.loc(Equivalence, NO, kernl.h = 2, 
                     epsilon = 1e-02, verb = TRUE)
out2 <- regmixEM.loc(Equivalence, NO, kernl.h = 2, k = 3,
                     epsilon = 1e-02, verb = TRUE)
c(out1$loglik, out2$loglik)

}

\keyword{file}
