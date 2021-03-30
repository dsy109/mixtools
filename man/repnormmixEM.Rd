\name{repnormmixEM}
\title{EM Algorithm for Mixtures of Normals with Repeated Measurements}
\alias{repnormmixEM}
\usage{
repnormmixEM(x, lambda = NULL, mu = NULL, sigma = NULL, k = 2, 
             arbmean = TRUE, arbvar = TRUE, epsilon = 1e-08, 
             maxit = 10000, verb = FALSE)
}

\description{
  Returns EM algorithm output for mixtures of normals with repeated
  measurements and arbitrarily many components.
}
\arguments{
  \item{x}{An mxn matrix of data.  The columns correspond to the subjects and
  the rows correspond to the repeated measurements.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1.  This determines number of components.  If NULL, then \code{lambda} is
    random from uniform Dirichlet and number of
    components is determined by \code{mu}.}
  \item{mu}{A k-vector of component means.  If NULL, then \code{mu} is determined by a
    normal distribution according to a binning method done on the data.  If both
    \code{lambda} and \code{mu} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{A vector of standard deviations.  If NULL, then \eqn{1/\code{sigma}^2} has
    random standard exponential entries according to a binning method done on the data.
    If \code{lambda}, \code{mu}, and \code{sigma} are NULL, then number of components is determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless all of \code{lambda}, \code{mu}, 
    and \code{sigma} are NULL.}
  \item{arbmean}{If TRUE, then the component densities are allowed to have different \code{mu}s. If FALSE, then
  a scale mixture will be fit.}
  \item{arbvar}{If TRUE, then the component densities are allowed to have different \code{sigma}s. If FALSE, then
  a location mixture will be fit.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
}
\value{
  \code{repnormmixEM} returns a list of class \code{mixEM} with items:
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
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{normalmixEM}}
}
\references{
  Hettmansperger, T. P. and Thomas, H.  (2000)
  Almost Nonparametric Inference for Repeated Measures in Mixture Models,
  \emph{Journal of the Royals Statistical Society, Series B} \bold{62(4)} 811--825.
}
\examples{
## EM output for the water-level task data set.

data(Waterdata)
set.seed(100)
water <- t(as.matrix(Waterdata[,3:10]))
em.out <- repnormmixEM(water, k = 2, verb = TRUE, epsilon = 1e-03)
em.out
}


\keyword{file}
