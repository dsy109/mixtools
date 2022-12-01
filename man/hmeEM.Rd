\name{hmeEM}
\title{EM Algorithm for Mixtures-of-Experts}
\alias{hmeEM}
\usage{
hmeEM(y, x, lambda = NULL, beta = NULL, sigma = NULL, w = NULL,
      k = 2, addintercept = TRUE, epsilon = 1e-08, 
      maxit = 10000, verb = FALSE)
}

\description{
  Returns EM algorithm output for a mixture-of-experts model. Currently, this code
  only handles a 2-component mixture-of-experts, but will be extended to the
  general k-component hierarchical mixture-of-experts.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An nxp matrix of predictors.  See \code{addintercept} below.}
  \item{lambda}{Initial value of mixing proportions, which are modeled as an inverse
    logit function of the predictors.  Entries should sum to 1.  
    If NULL, then \code{lambda} is taken as 1/\code{k} for each \code{x}.}
  \item{beta}{Initial value of \code{beta} parameters.  Should be a pxk matrix,
    where p is the number of columns of x and k is number of components.
    If NULL, then \code{beta} has standard normal entries according to a binning method done on the data.}
  \item{sigma}{A vector of standard deviations.  If NULL, then \eqn{1/\code{sigma}^2} has
    random standard exponential entries according to a binning method done on the data.}
  \item{w}{A p-vector of coefficients for the way the mixing proportions are modeled.  See \code{lambda}.}
  \item{k}{Number of components.  Currently, only \code{k}=2 is accepted.}
  \item{addintercept}{If TRUE, a column of ones is appended to the x
    matrix before the value of p is calculated.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{hmeEM} returns a list of class \code{mixEM} with items:
  \item{x}{The set of predictors (which includes a column of 1's if \code{addintercept} = TRUE).}
  \item{y}{The response values.}
  \item{w}{The final coefficients for the functional form of the mixing proportions.}
  \item{lambda}{An nxk matrix of the final mixing proportions.}
  \item{beta}{The final regression coefficients.}
  \item{sigma}{The final standard deviations. If \code{arbmean} = FALSE, then only the smallest standard
   deviation is returned. See \code{scale} below.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{regmixEM}}
}
\references{
  Jacobs, R. A., Jordan, M. I., Nowlan, S. J. and Hinton, G. E. (1991) Adaptive Mixtures of Local
  Experts, \emph{Neural Computation} \bold{3(1)}, 79--87.
  
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
}
\examples{
## EM output for NOdata.
 
data(NOdata)
attach(NOdata)
set.seed(100)
em.out <- regmixEM(Equivalence, NO)
hme.out <- hmeEM(Equivalence, NO, beta = em.out$beta)
hme.out[3:7]
}


\keyword{file}
