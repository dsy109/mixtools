\name{spregmix}
\title{EM-like Algorithm for Semiparametric Mixtures of Regressions}
\alias{spregmix}
\usage{
spregmix(lmformula, bw = NULL, constbw = FALSE,
         bwmult = 0.9, z.hat = NULL, symm = TRUE, betamethod = "LS",
         m = ifelse(is.null(z.hat), 2, ncol(z.hat)),
         epsilon = 1e-04, maxit = 1000, verbose = FALSE, 
         \dots) 
}
\description{
  Returns parameter estimates for finite mixtures of linear
  regressions with unspecified error structure.  Based on 
  Hunter and Young (2012).
}
\arguments{
  \item{lmformula}{Formula for a linear model, in the same format used by
  \code{\link{lm}}.  Additional parameters may be passed to \code{\link{lm}}
  via the \code{...} argument.}
  \item{bw}{Initial bandwidth value.  If NULL, this will be chosen automatically
  by the algorithm.}
  \item{constbw}{Logical:  If TRUE, the bandwidth is held constant throughout the
  algorithm; if FALSE, it adapts at each iteration according to the rules
  given in Hunter and Young (2012).}
  \item{bwmult}{
  Whenever it is updated automatically,
  the bandwidth is equal to \code{bwmult} divided by 
  the fifth root of \eqn{n} times the smaller of s and IQR/1.34, 
  where s and IQR are estimates of the standard deviation and interquartile 
  range of the residuals, as explained in Hunter and Young (2012).
  The value of 0.9 gives the rule of Silverman (1986) and
  the value of 1.06 gives the rule of Scott (1992).
  Larger values lead to greater smoothing, whereas smaller values lead to less 
  smoothing.  }
  \item{z.hat}{Initial nxm matrix of posterior probabilities.  If NULL, this
  is initialized randomly.  As long as a parametric estimation method like least
  squares is used to estimate \code{beta} in each M-step, the \code{z.hat}
  values are the only values necessary to begin the EM iterations.}
  \item{symm}{Logical:  If TRUE, the error density is assumed symmetric
  about zero.  If FALSE, it is not.  WARNING:  If FALSE, the intercept parameter
  is not uniquely identifiable if it is included in the linear model.}
  \item{betamethod}{Method of calculating beta coefficients in the
  M-step.  Current possible values are "LS" for least-squares; 
  "L1" for least absolute deviation; "NP" for fully nonparametric; 
  and "transition" for a transition
  from least squares to fully nonparametric. If something other than
  these four possibilities is used, then "NP" is assumed. For details
  of these methods, see Hunter and Young (2012).}
  \item{m}{Number of components in the mixture.}
  \item{epsilon}{Convergence is declared if the largest change in any lambda or
  beta coordinate is smaller than \code{epsilon}.}
  \item{maxit}{The maximum number of iterations; if convergence is never declared
  based on comparison with \code{epsilon}, then the algorithm stops after
  \code{maxit} iterations.}
  \item{verbose}{Logical: If TRUE, then various updates are printed during 
  each iteration of the algorithm.} 
  \item{\dots}{Additional parameters passed to the 
  \code{\link{model.frame}} and \code{\link{model.matrix}} functions,
  which are used to obtain the response and predictor of the regression.}
}
\value{
  \code{regmixEM} returns a list of class \code{npEM} with items:
  \item{x}{The set of predictors (which includes a column of 1's if \code{addintercept} = TRUE).}
  \item{y}{The response values.}
  \item{lambda}{The mixing proportions for every iteration in the form of a 
              matrix with m columns and (#iterations) rows}
  \item{beta}{The final regression coefficients.}
  \item{posterior}{An nxm matrix of posterior probabilities for
   observations.}
  \item{np.stdev}{Nonparametric estimate of the standard deviation, as given
  in Hunter and Young (2012)}
  \item{bandwidth}{Final value of the bandwidth}
  \item{density.x}{Points at which the error density is estimated}
  \item{density.y}{Values of the error density at the points \code{density.x}}
  \item{symmetric}{Logical:  Was the error density assumed symmetric?}
  \item{loglik}{A quantity similar to a log-likelihood, computed just like 
  a standard loglikelihood would be, conditional on the component density 
  functions being equal to the final density estimates.}
  \item{ft}{A character vector giving the name of the function.}
}
\seealso{
\code{\link{regmixEM}}, \code{\link{spEMsymloc}}, \code{\link{lm}}
}
\references{
  Hunter, D. R. and Young, D. S. (2012) Semi-parametric Mixtures of Regressions,
  Journal of Nonparametric Statistics 24(1): 19-38.

  Scott, D. W. (1992) \emph{Multivariate Density Estimation},
  John Wiley & Sons Inc., New York.

  Silverman, B. W. (1986). \emph{Density Estimation for Statistics and Data 
  Analysis}, Chapman & Hall, London.
}
\examples{
data(tonedata)
## By default, the bandwidth will adapt and the error density is assumed symmetric
set.seed(100)
a=spregmix(tuned~stretchratio, bw=.2, data=tonedata, verb=TRUE)

## Look at the sp mixreg solution:
plot(tonedata)
abline(a=a$beta[1,1],b=a$beta[2,1], col=2)
abline(a=a$beta[1,2],b=a$beta[2,2], col=3)

## Look at the nonparametric KD-based estimate of the error density, 
## constrained to be zero-symmetric:
plot(xx<-a$density.x, yy<-a$density.y, type="l")
## Compare to a normal density with mean 0 and NP-estimated stdev:
z <- seq(min(xx), max(xx), len=200)
lines(z, dnorm(z, sd=sqrt((a$np.stdev)^2+a$bandwidth^2)), col=2, lty=2)
# Add bandwidth^2 to variance estimate to get estimated var of KDE

## Now add the sp mixreg estimate without assuming symmetric errors:
b=spregmix(tuned~stretchratio, bw=.2, , symm=FALSE, data=tonedata, verb=TRUE)
lines(b$density.x, b$density.y, col=3)
}


\keyword{file}
