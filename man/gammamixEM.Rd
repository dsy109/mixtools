\name{gammamixEM}
\title{EM Algorithm for Mixtures of Gamma Distributions}
\alias{gammamixEM}
\usage{
gammamixEM(x, lambda = NULL, alpha = NULL, beta = NULL, k = 2,
           mom.start = TRUE, fix.alpha = FALSE, epsilon = 1e-08, 
           maxit = 1000, maxrestarts = 20, verb = FALSE)
}
\description{
  Return EM algorithm output for mixtures of gamma distributions.
}
\arguments{
  \item{x}{A vector of length n consisting of the data.}
  \item{lambda}{Initial value of mixing proportions.  If \code{NULL}, 
  then \code{lambda} is random from a uniform Dirichlet
  distribution (i.e., its entries are uniform random and then it is 
  normalized to sum to 1).}
  \item{alpha}{Starting value of vector of component shape parameters.  If non-NULL, \code{alpha} must be of length \code{k} if allowing different component shape parameters, or a single value if \code{fix.alpha = TRUE}.  If NULL, then the initial value
  is estimated by partitioning the data into \code{k} regions (with \code{lambda} determining
  the proportion of values in each region) and then calculating the method of moments estimates.}
  \item{beta}{Starting value of vector of component scale parameters.  If non-NULL and a vector,
  \code{k} is set to \code{length(beta)}.  If NULL, then the initial value
  is estimated the same method described for \code{alpha}.}
  \item{k}{Number of components.  Initial value ignored unless \code{alpha} and \code{beta}
    are both NULL.}
  \item{mom.start}{Logical to indicate if a method of moments starting value strategy should be implemented.  If \code{TRUE}, then only unspecified starting values will be generated according to this strategy.}
  \item{epsilon}{The convergence criterion.  Convergence is declared when the change in 
  the observed data log-likelihood increases by less than epsilon.}
  \item{fix.alpha}{Logical to indicate if the components should have a common shape parameter \code{alpha} estimated.  The default is \code{FALSE}.}
  \item{maxit}{The maximum number of iterations.}
  \item{maxrestarts}{The maximum number of restarts allowed in case of a problem
  with the particular starting values chosen (each restart uses randomly chosen
  starting values).}
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{gammamixEM} returns a list of class \code{mixEM} with items:
  \item{x}{The raw data.}
  \item{lambda}{The final mixing proportions.}
  \item{gamma.pars}{A 2xk matrix where each column provides the component estimates of \code{alpha} and \code{beta}.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.  This vector
  includes both the initial and the final values; thus, the number of iterations 
  is one less than its length.}
  \item{ft}{A character vector giving the name of the function.}
}
\references{
  Dempster, A. P., Laird, N. M., and Rubin, D. B. (1977) Maximum Likelihood From Incomplete Data
  Via the EM Algorithm, \emph{Journal of the Royal Statistical Society, Series B}, \bold{39(1)}, 1--38.
  
  Young, D. S., Chen, X., Hewage, D., and Nilo-Poyanco, R. (2019) Finite Mixture-of-Gamma Distributions: Estimation, Inference, and Model-Based Clustering, \emph{Advances in Data Analysis and Classification}, \bold{13(4)}, 1053--1082.
}
\examples{
##Analyzing a 3-component mixture of gammas.

set.seed(100)
x <- c(rgamma(200, shape = 0.2, scale = 14), rgamma(200, 
     shape = 32, scale = 10), rgamma(200, shape = 5, scale = 6))
out <- gammamixEM(x, lambda = c(1, 1, 1)/3, verb = TRUE)
out[2:4]
}

\keyword{file}
