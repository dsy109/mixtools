\name{summary.spRMM}
\alias{summary.spRMM}
\title{Summarizing fits from Stochastic EM algorithm for semiparametric scaled mixture of censored data}
\usage{
\method{summary}{spRMM}(object, digits = 6, ...)
}
\arguments{
  \item{object}{an object of class \code{spRMM} such as a result of a call
  to \code{\link{spRMM_SEM}}}
  \item{digits}{Significant digits for printing values}
  \item{...}{Additional parameters passed to \code{print}.}
}
\description{
\code{\link[base]{summary}} method for class \code{spRMM}.
}
\details{
 \code{\link{summary.spRMM}} prints scalar parameter estimates for
 a fitted mixture model: each component weight and the scaling factor, see reference below.
 The functional (nonparametric) estimates of survival and hazard rate funcions can be obtained 
 using \code{\link{plotspRMM}}.
}
\value{
  The function \code{\link{summary.spRMM}} prints the final loglikelihood
  value at the solution as well as The estimated mixing weights and the scaling parameter.
}

\seealso{
Function for plotting functional (nonparametric) estimates:
  \code{\link{plotspRMM}}.
  
  Other models and algorithms for censored lifetime data
  (name convention is model_algorithm):
  \code{\link{expRMM_EM}},
  \code{\link{weibullRMM_SEM}}.
}

\references{
   \itemize{
   \item Bordes, L., and Chauveau, D. (2016),
   Stochastic EM algorithms for parametric and semiparametric mixture models 
   for right-censored lifetime data, 
   Computational Statistics, Volume 31, Issue 4, pages 1513-1538.
   \url{https://link.springer.com/article/10.1007/s00180-016-0661-7}
   }
}
\author{Didier Chauveau}

\examples{
# See example(spRMM_SEM)
}
\keyword{file}

