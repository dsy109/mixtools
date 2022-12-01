\name{mvnormalmixEM}
\title{EM Algorithm for Mixtures of Multivariate Normals}
\alias{mvnormalmixEM}
\usage{
mvnormalmixEM(x, lambda = NULL, mu = NULL, sigma = NULL, k = 2,
              arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
              maxit = 10000, verb = FALSE)
}
\description{
  Return EM algorithm output for mixtures of multivariate normal distributions.
}
\arguments{
  \item{x}{A matrix of size nxp consisting of the data.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1.  This determines number of components.  If NULL, then \code{lambda} is
    random from uniform Dirichlet and number of
    components is determined by \code{mu}.}
  \item{mu}{A list of size k consisting of initial values for the p-vector mean parameters.  
    If NULL, then the vectors are generated from a normal distribution with
    mean and standard deviation according to a binning method done on the data.
    If both \code{lambda} and \code{mu} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{A list of size k consisting of initial values for the pxp variance-covariance matrices.  
    If NULL, then \code{sigma} is generated using the data.  
    If \code{lambda}, \code{mu}, and \code{sigma} are
    NULL, then number of components is determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless \code{lambda}, \code{mu}, and \code{sigma}
    are all NULL.}
  \item{arbmean}{If TRUE, then the component densities are allowed to have different \code{mu}s. If FALSE, then
   a scale mixture will be fit.}
  \item{arbvar}{If TRUE, then the component densities are allowed to have different \code{sigma}s. If FALSE, then
   a location mixture will be fit.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.}
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{normalmixEM} returns a list of class \code{mixEM} with items:
  \item{x}{The raw data.}
  \item{lambda}{The final mixing proportions.}
  \item{mu}{A list of with the final mean vectors.}
  \item{sigma}{A list with the final variance-covariance matrices.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
    observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
  \code{\link{normalmixEM}}
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
}
\examples{
##Fitting randomly generated data with a 2-component location mixture of bivariate normals.

set.seed(100)
x.1 <- rmvnorm(40, c(0, 0))
x.2 <- rmvnorm(60, c(3, 4))
X.1 <- rbind(x.1, x.2)
mu <- list(c(0, 0), c(3, 4))

out.1 <- mvnormalmixEM(X.1, arbvar = FALSE, mu = mu,
                       epsilon = 1e-02)
out.1[2:5]

##Fitting randomly generated data with a 2-component scale mixture of bivariate normals.

x.3 <- rmvnorm(40, c(0, 0), sigma = 
               matrix(c(200, 1, 1, 150), 2, 2))
x.4 <- rmvnorm(60, c(0, 0))
X.2 <- rbind(x.3, x.4)
lambda <- c(0.40, 0.60)
sigma <- list(diag(1, 2), matrix(c(200, 1, 1, 150), 2, 2))
 
out.2 <- mvnormalmixEM(X.2, arbmean = FALSE,
                       sigma = sigma, lambda = lambda,
                       epsilon = 1e-02)
out.2[2:5]
}

\keyword{file}

