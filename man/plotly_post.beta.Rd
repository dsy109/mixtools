\name{plotly_post.beta}
\title{Visualization of Posterior Regression Coefficients in Mixtures of Random Effects Regressions using \code{plotly}}
\alias{plotly_post.beta}
\usage{
plotly_post.beta(y, x, p.beta, p.z,
                 cex = 6,lwd=1,
                 title.size = 15,
                 xlab.size = 15 , xtick.size = 15,
                 ylab.size = 15 , ytick.size = 15,
                 col.data = "#1f77b4",
                 col.comp = NULL)
}

\description{
  Returns a 2x2 matrix of plots summarizing the posterior intercept and slope terms in a mixture of random effects regression with arbitrarily many components.
}
\arguments{
  \item{y}{A list of N response trajectories with (possibly) varying dimensions of
  length \eqn{n_i}.}
  \item{x}{A list of N predictor values of dimension \eqn{n_i}.  Each trajectory in y has
  its own design vector.}
  \item{p.beta}{A list of N 2xk matrices giving the posterior intercept and slope values from the output of an
  EM algorithm.}
  \item{p.z}{An Nxk matrix of posterior membership probabilities from the output of an EM algorithm.}
  \item{cex}{Size of dots of posterior Coefficients.}
  \item{lwd}{Width of lines.}
  \item{title.size}{Size of the main title.}
  \item{xlab.size}{Size of the lable of X-axis.}
  \item{xtick.size}{Size of tick lables of X-axis.}
  \item{ylab.size}{Size of the lable of Y-axis.}
  \item{ytick.size}{Size of tick lables of Y-axis.}
  \item{col.data}{Color of original data points.}
  \item{col.comp}{Color of points and lines of components. Number of colors specified needs to be consistent with number of components.}
}
\value{
  Plots returned.
}
\seealso{
\code{\link{regmixEM.mixed}}, \code{\link{plot.mixEM}}, \code{\link{post.beta}}
}
\references{
  Young, D. S. and Hunter, D. R. (2015) Random Effects Regression Mixtures for Analyzing Infant Habituation,
  \emph{Journal of Applied Statistics}, \bold{42(7)}, 1421--1441.
}
\examples{
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

x.1 = em.out$x
n = sum(sapply(x.1, nrow))
x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][,1])))
if (x.1.sum == n) {
  x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][,-1], ncol = 1))
} else {
  x = x.1
}

plotly_post.beta(x = x, y = em.out$y, p.beta = em.out$posterior.beta, 
                 p.z = em.out$posterior.z)
}
\details{
This is primarily used for within \code{plot.mixEM}.
}

\keyword{internal}
