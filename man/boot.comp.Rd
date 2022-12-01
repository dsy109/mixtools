\name{boot.comp}
\title{Performs Parametric Bootstrap for Sequentially Testing
the Number of Components in Various Mixture Models}
\alias{boot.comp}
\usage{
boot.comp(y, x = NULL, N = NULL, max.comp = 2, B = 100,
          sig = 0.05, arbmean = TRUE, arbvar = TRUE,
          mix.type = c("logisregmix", "multmix", "mvnormalmix",
          "normalmix", "poisregmix", "regmix", "regmix.mixed", 
          "repnormmix"), hist = TRUE, ...)
}

\description{
  Performs a parametric bootstrap by producing B bootstrap realizations of the likelihood ratio
  statistic for testing the null hypothesis of a k-component fit versus the alternative
  hypothesis of a (k+1)-component fit to various mixture models. This is performed for up to a
  specified number of maximum components, k. A p-value is calculated for each test and once the p-value
  is above a specified significance level, the testing terminates.  An optional histogram showing
  the distribution of the likelihood ratio statistic along with the observed statistic can also be
  produced.
}
\arguments{
  \item{y}{The raw data for \code{multmix}, \code{mvnormalmix}, \code{normalmix}, and \code{repnormmix} and
    the response values for \code{logisregmix}, \code{poisregmix}, and \code{regmix}.  See the documentation
    concerning their respective EM algorithms for specific structure of the raw data.}
  \item{x}{The predictor values required only for the regression mixtures \code{logisregmix}, 
  \code{poisregmix}, and \code{regmix}. A column of 1s for the intercept term must not be included! See the
    documentation concerning their respective EM algorithms for specific structure of the predictor values.}
  \item{N}{An n-vector of number of trials for the logistic regression type \code{logisregmix}.  
  If NULL, then \code{N} is an n-vector of 1s for binary logistic regression.}
  \item{max.comp}{The maximum number of components to test for.  The default is 2.  This function will 
    perform a test of k-components versus (k+1)-components sequentially until we fail to reject the null hypothesis.
    This decision rule is governed by the calculated p-value and \code{sig}.}
  \item{B}{The number of bootstrap realizations of the likelihood ratio statistic to produce.  The default is 100,
    but ideally, values of 1000 or more would be more acceptable.}
  \item{sig}{The significance level for which to compare the p-value against when performing the test of k-components
    versus (k+1)-components.}
  \item{arbmean}{If FALSE, then a scale mixture analysis can be performed for \code{mvnormalmix}, \code{normalmix}, \code{regmix},
    or \code{repnormmix}. The default is TRUE.}
  \item{arbvar}{If FALSE, then a location mixture analysis can be performed for \code{mvnormalmix}, \code{normalmix}, \code{regmix},
    or \code{repnormmix}. The default is TRUE.}
  \item{mix.type}{The type of mixture analysis you wish to perform. The data inputted for \code{y} and \code{x} depend on
    which type of mixture is selected.  \code{logisregmix} corresponds to a mixture of logistic regressions.  \code{multmix} corresponds
    to a mixture of multinomials with data determined by the cut-point method.  \code{mvnormalmix} corresponds to a mixture of
    multivariate normals.  \code{normalmix} corresponds to a mixture of univariate normals. \code{poisregmix} corresponds to a mixture of
    Poisson regressions.  \code{regmix} corresponds to a mixture of regressions with normal components.  \code{regmix.mixed}
    corresponds to a mixture of regressions with random or mixed effects.  \code{repnormmix} corresponds to a mixture
    of normals with repeated measurements.}
  \item{hist}{An argument to provide a matrix plot of histograms for the boostrapped likelihood ratio statistic.}
  \item{...}{Additional arguments passed to the various EM algorithms for the mixture of interest.} 
}
\value{
  \code{boot.comp} returns a list with items:
  \item{p.values}{The p-values for each test of k-components versus (k+1)-components.}
  \item{log.lik}{The B bootstrap realizations of the likelihood ratio statistic.}
  \item{obs.log.lik}{The observed likelihood ratio statistic for each test which is used in determining
    the p-values.}
}
\seealso{
\code{\link{logisregmixEM}}, \code{\link{multmixEM}}, \code{\link{mvnormalmixEM}}, \code{\link{normalmixEM}},
\code{\link{poisregmixEM}}, \code{\link{regmixEM}}, \code{\link{regmixEM.mixed}}, \code{\link{repnormmixEM}}
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
}
\examples{
## Bootstrapping to test the number of components on the RTdata.

data(RTdata)
set.seed(100)
x <- as.matrix(RTdata[, 1:3])
y <- makemultdata(x, cuts = quantile(x, (1:9)/10))$y
a <- boot.comp(y = y, max.comp = 1, B = 5, mix.type = "multmix", 
               epsilon = 1e-3)
a$p.values

}


\keyword{file}
