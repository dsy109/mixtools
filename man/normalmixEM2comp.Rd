\name{normalmixEM2comp}
\title{Fast EM Algorithm for 2-Component Mixtures of Univariate Normals}
\alias{normalmixEM2comp}
\usage{
normalmixEM2comp(x, lambda, mu, sigsqrd, eps= 1e-8, maxit = 1000, verb=FALSE)
}
\description{
  Return EM algorithm output for mixtures of univariate normal distributions 
  for the special case of 2 components, exploiting the simple structure of the
  problem to speed up the code.
}
\arguments{
  \item{x}{A vector of length \eqn{n} consisting of the data.}
  \item{lambda}{Initial value of first-component mixing proportion.}
  \item{mu}{A 2-vector of initial values for the mean parameters.}
  \item{sigsqrd}{Either a scalar or a 2-vector with initial value(s) for
  the variance parameters.  If a scalar, the algorithm assumes that the
  two components have equal variances; if a 2-vector, it assumes that the
  two components do not have equal variances.}
  \item{eps}{The convergence criterion.  Convergence is declared when the change in 
  the observed data log-likelihood increases by less than epsilon.}
  \item{maxit}{The maximum possible number of iterations.}
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\details{
  This code is written to be very fast, sometimes more than an order of magnitude
  faster than \code{\link{normalmixEM}} for the same problem.  It is less numerically
  stable that \code{\link{normalmixEM}} in the sense that it does not safeguard 
  against underflow as carefully.
  
  Note that when the two components are assumed to have unequal variances,
  the loglikelihood is unbounded.  However, in practice this is rarely a problem
  and quite often the algorithm converges to a "nice" local maximum.  Since this 
  algorithm requires starting values for all parameters, it is recommended that the 
  user initialize this routine from different starting values to gauge the reasonableness 
  of the outputted solutions.  Cases where, say, the components are not well-separated 
  and/or the component variances are somewhat similar, may have quite different estimates 
  depending on the starting values used.  Such situations could also be accompanied by slow 
  convergence of the algorithm relative to the convergence criteria used.
}
\value{
  \code{normalmixEM2comp} returns a list of class \code{mixEM} with items:
  \item{x}{The raw data.}
  \item{lambda}{The final mixing proportions (lambda and 1-lambda).}
  \item{mu}{The final two mean parameters.}
  \item{sigma}{The final one or two standard deviations. }
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nx2 matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.  This vector
  includes both the initial and the final values; thus, the number of iterations 
  is one less than its length.}
  \item{restarts}{The number of times the algorithm restarted due to 
  unacceptable choice of initial values (always zero).}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
  \code{\link{mvnormalmixEM}}, \code{\link{normalmixEM}}
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
}
\examples{
##Analyzing the Old Faithful geyser data with a 2-component mixture of normals.

data(faithful)
attach(faithful)
set.seed(100)
system.time(out <- normalmixEM2comp(waiting, lambda=.5, 
            mu=c(50,80), sigsqrd=100))
out$all.loglik # Note:  must be monotone increasing

# Compare elapsed time with more general version
system.time(out2 <- normalmixEM(waiting, lambda=c(.5,.5), 
            mu=c(50,80), sigma=c(10,10), arbvar=FALSE))
out2$all.loglik # Values should be identical to above
}

\keyword{file}
