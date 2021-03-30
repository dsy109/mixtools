\name{regmixMH}
\title{Metropolis-Hastings Algorithm for Mixtures of Regressions}
\alias{regmixMH}
\usage{
regmixMH(y, x, lambda = NULL, beta = NULL, s = NULL, k = 2,
         addintercept = TRUE, mu = NULL, sig = NULL, lam.hyp = NULL,
         sampsize = 1000, omega = 0.01, thin = 1)
}

\description{
  Return Metropolis-Hastings (M-H) algorithm output for mixtures of multiple regressions with
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
    If NULL, then \code{beta} has uniform standard normal entries.  If both
    \code{lambda} and \code{beta} are NULL, then number of components is determined by \code{s}.}
  \item{s}{k-vector of standard deviations.  If NULL, then \eqn{1/\code{s}^2} has
    random standard exponential entries.  If \code{lambda}, \code{beta}, and \code{s} are
    NULL, then number of components determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless all of \code{lambda}, \code{beta},
    and \code{s} are NULL.}
  \item{addintercept}{If TRUE, a column of ones is appended to the x
    matrix before the value of p is calculated.}
  \item{mu}{The prior hyperparameter of same size as \code{beta};
    the means of \code{beta} components.  If NULL,
    these are set to zero.}
  \item{sig}{The prior hyperparameter of same size as \code{beta};
    the standard deviations of \code{beta} components.  If NULL, these are 
    all set to five times the overall standard deviation of y.}
  \item{lam.hyp}{The prior hyperparameter of length \code{k} for the mixing proportions (i.e.,
	these are hyperparameters for the Dirichlet distribution).  If NULL, these are generated from a standard uniform
	distribution and then scaled to sum to 1.}
  \item{sampsize}{Size of posterior sample returned.}
  \item{omega}{Multiplier of step size to control M-H acceptance rate.
    Values closer to zero result in higher acceptance rates, generally.}
  \item{thin}{Lag between parameter vectors that will be kept.}
}
\value{
  \code{regmixMH} returns a list of class \code{mixMCMC} with items:
  \item{x}{A nxp matrix of the predictors.}
  \item{y}{A vector of the responses.}
  \item{theta}{A (\code{sampsize}/\code{thin}) x q matrix of MCMC-sampled
  q-vectors, where q is the total number of parameters in \code{beta}, \code{s}, and
  \code{lambda}.}
  \item{k}{The number of components.}
}
\seealso{
\code{\link{regcr}}
}
\references{
  Hurn, M., Justel, A. and Robert, C. P. (2003) Estimating Mixtures of Regressions, \emph{Journal 
  of Computational and Graphical Statistics} \bold{12(1)}, 55--79.
}
\examples{
## M-H algorithm for NOdata with acceptance rate about 40\%.

data(NOdata)
attach(NOdata)
set.seed(100)
beta <- matrix(c(1.3, -0.1, 0.6, 0.1), 2, 2)
sigma <- c(.02, .05)
MH.out <- regmixMH(Equivalence, NO, beta = beta, s = sigma, 
                   sampsize = 2500, omega = .0013)
MH.out$theta[2400:2499,]
}

\keyword{file}
