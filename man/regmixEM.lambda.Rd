\name{regmixEM.lambda}
\title{EM Algorithm for Mixtures of Regressions with Local Lambda Estimates}
\alias{regmixEM.lambda}
\usage{
regmixEM.lambda(y, x, lambda = NULL, beta = NULL, sigma = NULL, 
                k = 2, addintercept = TRUE, arbmean = TRUE,
                arbvar = TRUE, epsilon = 1e-8, maxit = 10000,
                verb = FALSE)
}

\description{
  Returns output for one step of an EM algorithm output for mixtures of multiple regressions 
  where the mixing proportions are estimated locally.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An nxp matrix of predictors.  See \code{addintercept} below.}
  \item{lambda}{An nxk matrix of initial local values of mixing proportions.  
    Entries should sum to 1.  This determines number of components.  
    If NULL, then \code{lambda} is simply one over the number of components.}
  \item{beta}{Initial value of \code{beta} parameters.  Should be a pxk matrix,
    where p is the number of columns of x and k is number of components.
    If NULL, then \code{beta} has uniform standard normal entries.  If both
    \code{lambda} and \code{beta} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{k-vector of initial global values of standard deviations.  
    If NULL, then \eqn{1/\code{sigma}^2} has random standard exponential entries.  
    If \code{lambda}, \code{beta}, and \code{sigma} are NULL, then number of components is determined by \code{k}.}
  \item{k}{The number of components.  Ignored unless all of \code{lambda}, \code{beta},
    and \code{sigma} are NULL.}
  \item{addintercept}{If TRUE, a column of ones is appended to the x
    matrix before the value of p is calculated.}
  \item{arbmean}{If TRUE, each mixture component is assumed to have a different set of regression coefficients
  (i.e., the \code{beta}s).}
  \item{arbvar}{If TRUE, each mixture component is assumed to have a different \code{sigma}.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{regmixEM.lambda} returns a list of class \code{mixEM} with items:
  \item{x}{The set of predictors (which includes a column of 1's if \code{addintercept} = TRUE).}
  \item{y}{The response values.}
  \item{lambda}{The inputted mixing proportions.}
  \item{beta}{The final regression coefficients.}
  \item{sigma}{The final standard deviations. If \code{arbmean} = FALSE, then only the smallest standard
   deviation is returned. See \code{scale} below.}
  \item{scale}{If \code{arbmean} = FALSE, then the scale factor for the component standard deviations is returned.
   Otherwise, this is omitted from the output.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\details{
    Primarily used within \code{regmixEM.loc}.
}
\seealso{
\code{\link{regmixEM.loc}}
}
\examples{
## Compare a 2-component and 3-component fit to NOdata.

data(NOdata)
attach(NOdata)
set.seed(100)
out1 <- regmixEM.lambda(Equivalence, NO)
out2 <- regmixEM.lambda(Equivalence, NO, k = 3)
c(out1$loglik, out2$loglik)

}

\keyword{file}
