\name{plotspRMM}
\alias{plotspRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot output from Stochastic EM algorithm for semiparametric scaled mixture of censored data
}

\description{Function for plotting various results from an object returned by \code{\link{spRMM_SEM}}, a Stochastic EM algorithm for semiparametric scaled 
mixture of randomly right censored lifetime data.
Four plots of sequences of estimates along iterations, survival and density estimates
(see reference below).
}

\usage{
  plotspRMM(sem, tmax = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sem}{An object returned by \code{\link{spRMM_SEM}}.}
  \item{tmax}{The max time for \eqn{x} axis, set to some default value if \code{NULL}.}
  }


\value{The four plots returned}


\seealso{
Related functions: \code{\link{spRMM_SEM}}.

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


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
# See example(spRMM_SEM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
