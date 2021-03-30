\name{regmixmodel.sel}
\title{Model Selection in Mixtures of Regressions}
\alias{regmixmodel.sel}
\usage{
regmixmodel.sel(x, y, w = NULL, k = 2, type = c("fixed", 
                "random", "mixed"), ...)
}
\description{
  Assess the number of components in a mixture of regressions model using the Akaike's information
  criterion (AIC), Schwartz's Bayesian information criterion (BIC), Bozdogan's consistent AIC (CAIC),
  and Integrated Completed Likelihood (ICL).
}
\arguments{
  \item{x}{An nxp matrix (or list) of predictors. If an intercept is required, then \code{x} must NOT include a column of 1's!
  Requiring an intercept may be controlled through arguments specified in \code{...}.}
  \item{y}{An n-vector (or list) of response values.}
  \item{w}{An optional list of fixed effects predictors for type "mixed" or "random".}
  \item{k}{The maximum number of components to assess.}
  \item{type}{The type of regression mixture to use.  If "fixed", then a mixture of regressions with fixed effects
  will be used.  If "random", then a mixture of regressions where the random effects regression coefficients are assumed
  to come from a mixture will be used.  If "mixed", the mixture structure used is the same as "random", except a coefficient
  of fixed effects is also assumed.}
  \item{...}{Additional arguments passed to the EM algorithm used for calculating the type of regression mixture specified
  in \code{type}.}
}
\value{
  \code{regmixmodel.sel} returns a matrix of the AIC, BIC, CAIC, and ICL values along with the winner (i.e., the highest
  value given by the model selection criterion) for various types of regression mixtures.
}
\seealso{
  \code{\link{regmixEM}}, \code{\link{regmixEM.mixed}}
}
\references{
  Biernacki, C., Celeux, G. and Govaert, G. (2000) Assessing a Mixture Model for Clustering with the
  Integrated Completed Likelihood, \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence} \bold{22(7)}, 719--725.

  Bozdogan, H. (1987) Model Selection and Akaike's Information Criterion (AIC): The General Theory and its
  Analytical Extensions, \emph{Psychometrika} \bold{52}, 345--370.
}
\examples{
## Assessing the number of components for NOdata.

data(NOdata)
attach(NOdata)
set.seed(100)
regmixmodel.sel(x = NO, y = Equivalence, k = 3, type = "fixed")
}

\keyword{file}
