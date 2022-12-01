\name{normalmixEM}
\title{EM Algorithm for Mixtures of Univariate Normals}
\alias{normalmixEM}
\usage{
normalmixEM(x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, 
            mean.constr = NULL, sd.constr = NULL,
            epsilon = 1e-08, maxit = 1000, maxrestarts = 20, 
            verb = FALSE, fast = FALSE, ECM = FALSE,
            arbmean = TRUE, arbvar = TRUE) 
}
\description{
  Return EM algorithm output for mixtures of normal distributions.
}
\details{
  This is the standard EM algorithm for normal mixtures that maximizes
  the conditional expected complete-data
  log-likelihood at each M-step of the algorithm.
  If desired, the
  EM algorithm may be replaced by an ECM algorithm (see \code{ECM} argument)
  that alternates between maximizing with respect to the \code{mu}
  and \code{lambda} while holding \code{sigma} fixed, and maximizing with
  respect to \code{sigma} and \code{lambda} while holding \code{mu}
  fixed.  In the case where \code{arbmean} is \code{FALSE}
  and \code{arbvar} is \code{TRUE}, there is no closed-form EM algorithm,
  so the ECM option is forced in this case.
}
\arguments{
  \item{x}{A vector of length n consisting of the data.}

  \item{lambda}{Initial value of mixing proportions.  Automatically 
  repeated as necessary 
  to produce a vector of length \code{k}, then normalized to sum to 1.
  If \code{NULL}, then \code{lambda} is random from a uniform Dirichlet
  distribution (i.e., its entries are uniform random and then it is 
  normalized to sum to 1).}

  \item{mu}{Starting value of vector of component means.  If non-NULL and a
  scalar, \code{arbmean} is set to \code{FALSE}.  If non-NULL and a vector,
  \code{k} is set to \code{length(mu)}.  If NULL, then the initial value
  is randomly generated from a normal distribution with center(s) determined
  by binning the data.}

  \item{sigma}{Starting value of vector of component standard deviations 
  for algorithm.  If non-NULL
  and a scalar, \code{arbvar} is set to \code{FALSE}.  If non-NULL and a vector,
  \code{arbvar} is set to \code{TRUE} and \code{k} is set to \code{length(sigma)}.
  If NULL, then the initial value is the reciprocal of the square root of
  a vector of random exponential-distribution values whose means are determined
  according to a binning method done on the data.}

  \item{k}{Number of components.  Initial value ignored unless \code{mu} 
  and \code{sigma} are both NULL.}

  \item{mean.constr}{Equality constraints on the mean parameters, given as
  a vector of length \code{k}.  Each vector entry helps specify the constraints,
  if any, on the corresponding mean parameter:  If \code{NA}, the corresponding
  parameter is unconstrained.  If numeric, the corresponding
  parameter is fixed at that value.  If a character string consisting of
  a single character preceded by a coefficient, such as \code{"0.5a"}
  or \code{"-b"}, all parameters using the same single character in their
  constraints will fix these parameters equal to the coefficient times
  some the same free parameter.  For instance, if 
  \code{mean.constr = c(NA, 0, "a", "-a")}, then the first mean parameter
  is unconstrained, the second is fixed at zero, and the third and forth
  are constrained to be equal and opposite in sign.}
  
  \item{sd.constr}{Equality constraints on the standard deviation parameters.
  See \code{mean.constr}.}

  \item{epsilon}{The convergence criterion.  Convergence is declared when the change in 
  the observed data log-likelihood increases by less than epsilon.}

  \item{maxit}{The maximum number of iterations.}

  \item{maxrestarts}{The maximum number of restarts allowed in case of a problem
  with the particular starting values chosen due to one of the variance
  estimates getting too small
  (each restart uses randomly chosen
  starting values).  It is well-known that when each component of a normal
  mixture may have its own mean and variance, the likelihood has no maximizer;
  in such cases, we hope to find a "nice" local maximum with this algorithm
  instead, but occasionally the algorithm finds a "not nice" solution and
  one of the variances goes to zero, driving the likelihood to infinity.}

  \item{verb}{If TRUE, then various updates are printed during each 
  iteration of the algorithm.} 

  \item{fast}{If TRUE and k==2 and arbmean==TRUE, then use 
  \code{\link{normalmixEM2comp}}, which is a much faster version of the EM 
  algorithm for this case.
  This version is less protected against certain kinds of underflow
  that can cause numerical problems and it does not permit any restarts.  If
  k>2, \code{fast} is ignored.}
  
  \item{ECM}{logical:  Should this algorithm be an ECM algorithm in the sense
  of Meng and Rubin (1993)?  If FALSE, the algorithm is a true EM algorithm;
  if TRUE, then every half-iteration alternately updates the means conditional
  on the variances or the variances conditional on the means, with an extra
  E-step in between these updates.}

  \item{arbmean}{If TRUE, then the component densities are allowed to have different \code{mu}s. If FALSE, then
  a scale mixture will be fit.  Initial value ignored unless \code{mu} is NULL.}

  \item{arbvar}{If TRUE, then the component densities are allowed to have different \code{sigma}s. If FALSE, then
  a location mixture will be fit.  Initial value ignored unless \code{sigma} is NULL.}

}
\value{
  \code{normalmixEM} returns a list of class \code{mixEM} with items:
  \item{x}{The raw data.}
  \item{lambda}{The final mixing proportions.}
  \item{mu}{The final mean parameters.}
  \item{sigma}{The final standard deviations. If \code{arbmean} = FALSE, then only the smallest standard
   deviation is returned. See \code{scale} below.}
  \item{scale}{If \code{arbmean} = FALSE, then the scale factor for the component standard deviations is returned.
   Otherwise, this is omitted from the output.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.  This vector
  includes both the initial and the final values; thus, the number of iterations 
  is one less than its length.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
  \code{\link{mvnormalmixEM}}, \code{\link{normalmixEM2comp}},
  \code{\link{normalmixMMlc}}, \code{\link{spEMsymloc}}
}
\references{
  \itemize{
  \item McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, 
  John Wiley and Sons, Inc.
  \item Meng, X.-L. and Rubin, D. B. (1993) Maximum Likelihood Estimation
  Via the ECM Algorithm:  A General Framework, \emph{Biometrika} 80(2):
  267-278.
  
  \item Benaglia, T., Chauveau, D., Hunter, D. R., and Young, D.
	mixtools: An R package for analyzing finite mixture models.
	Journal of Statistical Software, 32(6):1-29, 2009.

  }
}
\examples{
##Analyzing the Old Faithful geyser data with a 2-component mixture of normals.

data(faithful)
attach(faithful)
set.seed(100)
system.time(out<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03))
out
system.time(out2<-normalmixEM(waiting, arbvar = FALSE, epsilon = 1e-03, fast=TRUE))
out2 # same thing but much faster
}

\keyword{file}
