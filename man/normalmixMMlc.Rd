\name{normalmixMMlc}
\title{EC-MM Algorithm for Mixtures of Univariate Normals 
	with linear constraints}
\alias{normalmixMMlc}
\usage{
normalmixMMlc(x, lambda = NULL, mu = NULL, sigma = NULL, k = 2,
              mean.constr = NULL, mean.lincstr = NULL, 
              mean.constant = NULL, var.lincstr = NULL, 
              gparam = NULL, epsilon = 1e-08, maxit = 1000, 
              maxrestarts=20, verb = FALSE) 
}
\description{
  Return EC-MM (see below) algorithm output for mixtures of normal distributions
  with linear constraints on the means and variances parameters,
  as in Chauveau and Hunter (2013).
  The linear constraint for the means is of the form
  \eqn{\mu = M \beta + C}, where \eqn{M} and \eqn{C} are matrix 
  and vector specified as parameters.
  The linear constraints for the variances are actually specified on
  the inverse variances, by \eqn{\pi = A \gamma}, where \eqn{\pi}
  is the vector of inverse variances, and \eqn{A} is a matrix
  specified as a parameter (see below).
}
\details{
  This is a specific "EC-MM" algorithm for normal mixtures 
  with linear constraints on the means and variances parameters.
  EC-MM here means that this algorithm is similar to 
  an ECM algorithm as in Meng and Rubin (1993), 
  except that it uses conditional MM 
  (Minorization-Maximization)-steps
  instead of simple M-steps. Conditional means that it 
  alternates between maximizing with respect to the \code{mu}
  and \code{lambda} while holding \code{sigma} fixed, and maximizing with
  respect to \code{sigma} and \code{lambda} while holding \code{mu}
  fixed.  This ECM generalization of EM is forced in the case of linear constraints because there is no closed-form EM algorithm.
}
\arguments{
  \item{x}{A vector of length n consisting of the data.}

  \item{lambda}{Initial value of mixing proportions.  Automatically 
  repeated as necessary 
  to produce a vector of length \code{k}, then normalized to sum to 1.
  If \code{NULL}, then \code{lambda} is random from a uniform Dirichlet
  distribution (i.e., its entries are uniform random and then it is 
  normalized to sum to 1).}

  \item{mu}{Starting value of vector of component means.
  	  If non-NULL and a vector,
  \code{k} is set to \code{length(mu)}.  If NULL, then the initial value
  is randomly generated from a normal distribution with center(s) determined
  by binning the data.}

  \item{sigma}{Starting value of vector of component standard deviations 
  for algorithm.  
  Obsolete for linear constraints on the inverse variances;
  use \code{gparam} instead to specify a starting value.}

  \item{k}{Number of components.  Initial value ignored unless \code{mu} 
  and \code{sigma} are both NULL.}

  \item{mean.constr}{First, simplest way to define 
  	  equality constraints on the mean parameters, given as
  a vector of length \code{k}, as in \code{\link{normalmixEM}}.  
  Each vector entry specifies the constraints,
  if any, on the corresponding mean parameter:  If \code{NA}, the corresponding
  parameter is unconstrained.  If numeric, the corresponding
  parameter is fixed at that value.  If a character string consisting of
  a single character preceded by a coefficient, such as \code{"0.5a"}
  or \code{"-b"}, all parameters using the same single character in their
  constraints will fix these parameters equal to the coefficient times
  some the same free parameter.  For instance, if 
  \code{mean.constr = c(NA, 0, "a", "-a")}, then the first mean parameter
  is unconstrained, the second is fixed at zero, and the third and forth
  are constrained to be equal and opposite in sign.
  Note: if there are no linear constraints for the means, it is
  more efficient to use directly \code{\link{normalmixEM}}.}
  
  \item{mean.lincstr}{Matrix \eqn{M} \eqn{(k,p)} in the linear constraint for the means
  	equation \eqn{\mu = M \beta + C}, with \eqn{p \le k}.}

  \item{mean.constant}{Vector of \eqn{k} constants \eqn{C} 
  			in the linear constraint for the means
  			equation \eqn{\mu = M \beta + C}.}

  \item{var.lincstr}{Matrix \eqn{A} \eqn{(k,q)} in the linear constraint for the
  	  inverse variances equation \eqn{\pi = A \gamma}, with \eqn{q \le k}.}
  
  	\item{gparam}{Vector of \eqn{q} starting values
  		for the \eqn{\gamma} parameter in the
  	  linear constraint for the inverse variances;
  	  see \code{var.lincstr}.  If NULL, a vector of randomly generated
    	standard exponential variables is used.}
  
  \item{epsilon}{The convergence criterion.  
  	  Convergence is declared when the change in 
 	 the observed data log-likelihood increases by less than epsilon.}

  \item{maxit}{The maximum allowed number of iterations.}

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

}
\value{
  \code{normalmixMMlc} returns a list of class \code{mixEM} with items:
  \item{x}{The raw data.}
  \item{lambda}{The final mixing proportions.}
  \item{mu}{The final mean parameters.}
  \item{sigma}{The final standard deviation(s)}
  \item{scale}{Scale factor for the component standard deviations, if applicable.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.  This vector
  includes both the initial and the final values; thus, the number of iterations 
  is one less than its length.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{beta}{The final \eqn{\beta} parameter estimate.}
  \item{gamma}{The final \eqn{\gamma} parameter estimate.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
 \code{\link{normalmixEM}}, \code{\link{mvnormalmixEM}}, 
 \code{\link{normalmixEM2comp}}, \code{\link{tauequivnormalmixEM}}
}
\references{
  \itemize{
  \item McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, 
  John Wiley & Sons, Inc.
  \item Meng, X.-L. and Rubin, D. B. (1993) Maximum Likelihood Estimation
  Via the ECM Algorithm:  A General Framework, \emph{Biometrika} 80(2):
  267-278.
  \item Chauveau, D. and Hunter, D.R. (2013)
  	ECM and MM algorithms for mixtures with constrained parameters,
  	\emph{preprint \url{https://hal.archives-ouvertes.fr/hal-00625285}}.
  \item Thomas, H., Lohaus, A., and Domsch, H. (2011) Stable Unstable Reliability
        Theory, \emph{British Journal of Mathematical and Statistical Psychology}
        65(2): 201-221.
  }
}


\author{Didier Chauveau}


\examples{
## Analyzing synthetic data as in the tau equivalent model  
## From Thomas et al (2011), see also Chauveau and Hunter (2013)
## a 3-component mixture of normals with linear constraints.
lbd <- c(0.6,0.3,0.1); m <- length(lbd)
sigma <- sig0 <- sqrt(c(1,9,9))
# means constaints mu = M beta
M <- matrix(c(1,1,1,0,-1,1), 3, 2)
beta <- c(1,5) # unknown constrained mean
mu0 <- mu <- as.vector(M \%*\% beta)
# linear constraint on the inverse variances pi = A.g
A <- matrix(c(1,1,1,0,1,0), m, 2, byrow=TRUE)
iv0 <- 1/(sig0^2)
g0 <- c(iv0[2],iv0[1] - iv0[2]) # gamma^0 init 

# simulation and EM fits
set.seed(50); n=100; x <- rnormmix(n,lbd,mu,sigma)
s <- normalmixEM(x,mu=mu0,sigma=sig0,maxit=2000) # plain EM
# EM with var and mean linear constraints
sc <- normalmixMMlc(x, lambda=lbd, mu=mu0, sigma=sig0,
					mean.lincstr=M, var.lincstr=A, gparam=g0)
# plot and compare both estimates
dnormmixt <- function(t, lam, mu, sig){
	m <- length(lam); f <- 0
	for (j in 1:m) f <- f + lam[j]*dnorm(t,mean=mu[j],sd=sig[j])
	f}
t <- seq(min(x)-2, max(x)+2, len=200)
hist(x, freq=FALSE, col="lightgrey", 
		ylim=c(0,0.3), ylab="density",main="")
lines(t, dnormmixt(t, lbd, mu, sigma), col="darkgrey", lwd=2) # true
lines(t, dnormmixt(t, s$lambda, s$mu, s$sigma), lty=2) 
lines(t, dnormmixt(t, sc$lambda, sc$mu, sc$sigma), col=1, lty=3)
legend("topleft", c("true","plain EM","constr EM"), 
	col=c("darkgrey",1,1), lty=c(1,2,3), lwd=c(2,1,1))
}

\keyword{file}
