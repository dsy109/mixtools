\name{segregmixEM}
\title{ECM Algorithm for Mixtures of Regressions with Changepoints}
\alias{segregmixEM}
\alias{regmixEM.chgpt}
\usage{
segregmixEM(y, x, lambda = NULL, beta = NULL, sigma = NULL, 
            k = 2, seg.Z, psi, psi.locs = NULL, delta = NULL, 
            epsilon = 1e-08, maxit = 10000, verb = FALSE,
            max.restarts = 15)
}

\description{
  Returns ECM algorithm output for mixtures of multiple regressions with
  changepoints and arbitrarily many components.
}
\arguments{
  \item{y}{An n-vector of response values.}
  \item{x}{An nxp matrix of predictors.  Note that this model assumes the presence of an intercept.}
  \item{lambda}{Initial value of mixing proportions.  Entries should sum to
    1.  This determines number of components.  If NULL, then \code{lambda} is
    random from uniform Dirichlet and the number of
    components is determined by \code{beta}.}
  \item{beta}{Initial value of \code{beta} parameters.  This is a list of length \code{k} such that each element must contain a vector having length consistent with the defined changepoint structure. See \code{seg.Z}, \code{psi}, and \code{psi.loc} below. If NULL, then \code{beta} has standard normal entries according to a binning method done on the data.  If both
    \code{lambda} and \code{beta} are NULL, then number of components is determined by \code{sigma}.}
  \item{sigma}{A vector of standard deviations.  If NULL, then 1/\code{sigma}^2 has
    random standard exponential entries according to a binning method done on the data.
    If \code{lambda}, \code{beta}, and \code{sigma} are
    NULL, then number of components is determined by \code{k}.}
  \item{k}{Number of components.  Ignored unless all of \code{lambda}, \code{beta}, 
  and \code{sigma} are NULL.}
  \item{seg.Z}{A list of length \code{k} whose elements are right-hand side formulas, which are additive linear models of the predictors that have changepoints in their respective components.  See below for more details.}
  \item{psi}{A kxp matrix specifying the number of changepoints for each predictor in each component. See below for more details.}
  \item{psi.locs}{A list of length \code{k} that has initial estimates for the changepoint locations.  Each element of the list must have length equal to the number of chanegpoints specified in the corresponding row of the \code{psi} matrix.  For components with no changepoints, simply set that element equal to NULL. See below for more details.}
  \item{delta}{An optional list of values quantifying the amount of separation at each changepoint if assuming discontinuities at the changepoints.  This has the same dimensions as \code{psi.locs}.}
  \item{epsilon}{The convergence criterion.}
  \item{maxit}{The maximum number of iterations.} 
  \item{verb}{If TRUE, then various updates are printed during each iteration of the algorithm.} 
  \item{max.restarts}{The number of times to try restarting the ECM algorithm if estimation problems occur - such as choice of poor initial values or a poorly chosen changepoint structure.}
}
\value{
  \code{segregmixEM} returns a list of class \code{segregmixEM} with items:
  \item{x}{The set of predictors.}
  \item{y}{The response values.}
  \item{lambda}{The final mixing proportions.}
  \item{beta}{The final regression coefficients.}
  \item{sigma}{The final standard deviations.}
  \item{seg.Z}{The list of right-hand side formulas as defined by the user.}
	\item{psi.locs}{A list of length k with the final estimates for the changepoint locations.}
	\item{delta}{A list of the delta values that were optionally specified by the user.}
  \item{loglik}{The final log-likelihood.}
  \item{posterior}{An nxk matrix of posterior probabilities for
   observations.}
  \item{all.loglik}{A vector of each iteration's log-likelihood.}
  \item{restarts}{The number of times the algorithm restarted due to unacceptable choice of initial values.}
  \item{ft}{A character vector giving the name of the function.}
}
\details{
\code{seg.Z} is defined as a list of right-hand side linear model formulas that are used to identify which predictors have changepoints in each component.  For example, suppose you have a dataframe with three predictors: \code{V1}, \code{V2}, \code{V3}.  Suppose now that you wish to model a 3-component mixture of regressions with changepoints structure such that the first component has changepoints in V1 and V2, the second component has changepoints in \code{V3}, and the third component has no changepoints.  Then you would define \code{seg.Z = list(~V1+V2, ~V3, NULL)}.  Note that you MUST place the variables in order with respect to how they appear in the predictor matrix \code{x}.

\code{psi} is a kxp matrix specifying the number of changepoints for each predictor in each component.  For the example given above, suppose there are three changepoints for \code{V1}, two changepoints for \code{V2}, and four changepoints for \code{V3}.  Then you would define \code{psi = rbind(c(3, 2, 0), c(0, 0, 4), c(0, 0, 0))}.

\code{psi.locs} is a list of length k whose elements give the initial locations of the changepoints for each component.  Each element of the list must have length equal to the total number of changepoints for that component's regression equation.  For the example given above, in component 1, assume that the three changepoints for \code{V1} are at 3, 7, and 10 and the two changepoints for \code{V1} are at 5, 20, and 30.  In component 2, assume that the four changepoints for \code{V3} are at 2, 4, 6, and 8.  Then you would define \code{psi.locs = list(c(3, 7, 10, 5, 20, 30), c(2, 4, 6, 8), NULL)}.  Note that the order of the changepoints is determined by first sorting the predictors by how they appear in the formulas in \code{seg.Z} and then sorting in increasing order within each predictor.
}
\seealso{
\code{\link{regmixEM}}
}
\references{
  Young, D. S. (2014) Mixtures of Regressions with Changepoints, \emph{Statistics and Computing}, \bold{24(2)}, 265--281.
}
\note{
As of version 0.4.6, this more general function has replaced the now defunct \code{regmixEM.chgpt} and associated internal functions.
}
\examples{
\dontrun{
## Simulated example.

set.seed(100)
x <- 1:20
y1 <- 3 + x + rnorm(20)
y2 <- 3 - x - 5*(x - 15)*(x > 15) + rnorm(20)
y <- c(y1, y2)
x <- c(x, x)

set.seed(100)
be <- list(c(3, -1, -5), c(3, 1))
s <- c(1, 1)
psi.locs <- list(comp.1 = list(x = 15), comp.2 = NULL)
out <- segregmixEM(y, cbind(1,x), verb = TRUE, k = 2,
                   beta = be, sigma = s, lambda = c(1, 1)/2, 
                   seg.Z = list(~x, NULL), psi = rbind(1, 0),
                   psi.locs = psi.locs, epsilon = 0.9)

z <- seq(0, 21, len = 40)
plot(x, y, col = apply(out$post, 1, which.max) + 1, pch = 19, 
	   cex.lab = 1.4, cex = 1.4)
b <- out$beta
d <- out$psi.locs
lines(z, b[[1]][1] + b[[1]][2] * z + b[[1]][3] * 
      (z - d[[1]][[1]]) * (z > d[[1]][[1]]) , col = 2, lwd = 2)
lines(z, b[[2]][1] + b[[2]][2] * z, col = 3, lwd = 2)
abline(v = out$psi.locs[[1]][1], col = 2, lty = 2)
}

\dontrun{
## Example using the NOdata.
 
data(NOdata)
attach(NOdata)

set.seed(100)
be <- list(c(1.30, -0.13, 0.08), c(0.56, 0.09))
s <- c(0.02, 0.04)
psi.locs <- list(comp.1 = list(NO = 1.57), comp.2 = NULL)
out <- segregmixEM(Equivalence, cbind(NO), verb = TRUE, k = 2,
                   beta = be, sigma = s, lambda = c(1, 1)/2, 
                   seg.Z = list(~NO, NULL), psi = rbind(1, 0),
                   psi.locs = psi.locs, epsilon = 0.1)

z <- seq(0, 5, len = 1000)
plot(NOdata, col = apply(out$post, 1, which.max) + 1, pch = 19, 
	   cex.lab = 1.4, cex = 1.4, ylab = "Equivalence Ratio")
b <- out$beta
d <- out$psi.locs
lines(z, b[[1]][1] + b[[1]][2] * z + b[[1]][3] * 
      (z - d[[1]][[1]]) * (z > d[[1]][[1]]) , col = 2, lwd = 2)
lines(z, b[[2]][1] + b[[2]][2] * z, col = 3, lwd = 2)
abline(v = out$psi.locs[[1]][1], col = 2, lty = 2)

detach(NOdata)
}
}


\keyword{file}
