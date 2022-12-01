\name{regmixEM}
\title{EM Algorithm for Mixtures of Regressions}
\alias{regmixEM}
\usage{
regmixEM(y, x, lambda = NULL, beta = NULL, sigma = NULL, k = 2,
         addintercept = TRUE, arbmean = TRUE, arbvar = TRUE, 
         epsilon = 1e-08, maxit = 10000, verb = FALSE)
}

\description{
  Returns EM algorithm output for mixtures of multiple regressions with
  arbitrarily many components.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An nxp matrix of predictors.  See \code{addintercept} below.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1.  This determines number of components.  If NULL, then \code{lambda} is
    random from uniform Dirichlet and number of
    components is determined by \code{beta}.}
  \item{beta}{Initial value of \code{beta} parameters.  Should be a pxk matrix,
    where p is the number of columns of x and k is number of components.
    If NULL, then \code{beta} has standard normal entries according to a binning method done on the data.  If both
    \code{lambda} and \code{beta} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{A vector of standard deviations.  If NULL, then 1/\code{sigma}^2 has
    random standard exponential entries according to a binning method done on the data.
    If \code{lambda}, \code{beta}, and \code{sigma} are
    NULL, then number of components is determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless all of \code{lambda}, \code{beta}, 
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
  \code{regmixEM} returns a list of class \code{mixEM} with items:
  \item{x}{The set of predictors (which includes a column of 1's if \code{addintercept} = TRUE).}
  \item{y}{The response values.}
  \item{lambda}{The final mixing proportions.}
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
\seealso{
\code{\link{regcr}}, \code{\link{regmixMH}}
}
\references{
de Veaux, R. D. (1989), Mixtures of Linear Regressions,
  \emph{Computational Statistics and Data Analysis} 8, 227-245.

  Hurn, M., Justel, A. and Robert, C. P. (2003) Estimating Mixtures of Regressions, \emph{Journal 
  of Computational and Graphical Statistics} \bold{12(1)}, 55--79.
  
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
}
\examples{
## EM output for NOdata.
 
data(NOdata)
attach(NOdata)
set.seed(100)
em.out <- regmixEM(Equivalence, NO, verb = TRUE, epsilon = 1e-04)
em.out[3:6]
}


\keyword{file}
