\name{multmixEM}
\title{EM Algorithm for Mixtures of Multinomials}
\alias{multmixEM}
\usage{
multmixEM(y, lambda = NULL, theta = NULL, k = 2,
          maxit = 10000, epsilon = 1e-08, verb = FALSE)
}
\description{
  Return EM algorithm output for mixtures of multinomial distributions.
}
\arguments{
  \item{y}{Either An nxp matrix of data (multinomial counts), where n is the
  sample size and p is the number of multinomial bins, or the
  output of the \code{\link{makemultdata}} function.  It is not necessary
  that all of the rows contain the same number of multinomial trials (i.e.,
  the row sums of \code{y} need not be identical).}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1.  This determines number of components.  If NULL, then \code{lambda} is
    random from uniform Dirichlet and number of
    components is determined by \code{theta}.}
  \item{theta}{Initial value of \code{theta} parameters.  Should be a kxp matrix,
    where p is the number of columns of y and k is number of components.
    Each row of \code{theta} should sum to 1.
    If NULL, then each row is random from uniform Dirichlet.
    If both \code{lambda} and \code{theta} are NULL, then number of components 
    is determined by k.}
  \item{k}{Number of components.  Ignored unless \code{lambda} and \code{theta}
    are NULL.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.}
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{multmixEM} returns a list of class \code{mixEM} with items:
  \item{y}{The raw data.}
  \item{lambda}{The final mixing proportions.}
  \item{theta}{The final multinomial parameters.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
    observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.} 
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{compCDF}}, \code{\link{makemultdata}}, \code{\link{multmixmodel.sel}}
}
\references{
  \itemize{
  \item McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
  
  \item Elmore, R. T., Hettmansperger, T. P. and Xuan, F. (2004) The Sign Statistic, One-Way Layouts
  and Mixture Models, \emph{Statistical Science} \bold{19(4)}, 579--587.
  }
}
\examples{
## The sulfur content of the coal seams in Texas

set.seed(100)
A <- c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B <- c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59) 
C <- c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49) 
D <- c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32) 
E <- c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)

dis.coal <- makemultdata(A, B, C, D, E, 
                         cuts = median(c(A, B, C, D, E)))
em.out <- multmixEM(dis.coal)
em.out[1:4]
}

\keyword{file}
