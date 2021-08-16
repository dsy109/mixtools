\name{post.beta}
\title{Summary of Posterior Regression Coefficients in Mixtures of Random Effects Regressions}
\alias{post.beta}
\usage{
post.beta(y, x, p.beta, p.z)
}

\description{
  Returns a 2x2 matrix of plots summarizing the posterior intercept and slope terms in a mixture of
  random effects regression with arbitrarily many components.
}
\arguments{
  \item{y}{A list of N response trajectories with (possibly) varying dimensions of
  length \eqn{n_i}.}
  \item{x}{A list of N predictor values of dimension \eqn{n_i}.  Each trajectory in y has
  its own design vector.}
  \item{p.beta}{A list of N 2xk matrices giving the posterior intercept and slope values from the output of an
  EM algorithm.}
  \item{p.z}{An Nxk matrix of posterior membership probabilities from the output of an EM algorithm.}
}
\value{
  \code{post.beta} returns a 2x2 matrix of plots giving:
  \item{(1, 1)}{The data plotted on the x-y axes with all posterior regression lines.}
  \item{(1, 2)}{The data plotted on the x-y axes with most probable posterior regression lines.}
  \item{(2, 1)}{A beta-space plot of all posterior regression coefficients.}
  \item{(1, 1)}{A beta-space plot of most probable posterior regression coefficients.}
}
\seealso{
\code{\link{regmixEM.mixed}}, \code{\link{plot.mixEM}}
}
\references{
  Young, D. S. and Hunter, D. R. (2015) Random Effects Regression Mixtures for Analyzing Infant Habituation,
  \emph{Journal of Applied Statistics}, \bold{42(7)}, 1421--1441.
}
\examples{
  \dontrun{
## EM output for simulated data from 2-component mixture of random effects.

data(RanEffdata)
set.seed(100)
x <- lapply(1:length(RanEffdata), function(i) 
            matrix(RanEffdata[[i]][, 2:3], ncol = 2))
x <- x[1:20]
y <- lapply(1:length(RanEffdata), function(i) 
            matrix(RanEffdata[[i]][, 1], ncol = 1))
y <- y[1:20]
lambda <- c(0.45, 0.55)
mu <- matrix(c(0, 4, 100, 12), 2, 2)
sigma <- 2
R <- list(diag(1, 2), diag(1, 2))
em.out <- regmixEM.mixed(y, x, sigma = sigma, arb.sigma = FALSE,
                         lambda = lambda, mu = mu, R = R,
                         addintercept.random = FALSE,
                         epsilon = 1e-02, verb = TRUE)

## Obtaining the 2x2 matrix of plots.

x.ran <- lapply(1:length(x), function(i) x[[i]][, 2])
p.beta <- em.out$posterior.beta
p.z <- em.out$posterior.z
post.beta(y, x.ran, p.beta = p.beta, p.z = p.z)
}
}
\details{
This is primarily used for within \code{plot.mixEM}.
}

\keyword{file}
