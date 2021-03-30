\name{summary.mixEM}
\alias{summary.mixEM}
\title{Summarizing EM mixture model fits}
\usage{
\method{summary}{mixEM}(object, digits=6, \dots)
}
\arguments{
  \item{object}{an object of class \code{mixEM} such as a result of a call
  to \code{\link{normalmixEM}}}
  \item{digits}{Significant digits for printing values}
  \item{\dots}{further arguments passed to \code{print} method.}
}
\description{
\code{\link[base]{summary}} method for class \code{mixEM}.
}
\details{
 \code{\link{summary.mixEM}} prints parameter estimates for
 each component of a fitted mixture model.  
 The estimates printed vary with the type of model.
}
\value{
  The function \code{\link{summary.mixEM}} prints the final loglikelihood
  value at the solution as well as a matrix of values for each component
  that could include:
  \item{lambda}{The estimated mixing weights}
  \item{mu}{The estimated mean parameters}
  \item{sigma}{The estimated standard deviations}
  \item{theta}{The estimated multinomial parameters}
  \item{beta}{The estimated regression parameters}
}
\seealso{ 
\code{\link{normalmixEM}},
\code{\link{logisregmixEM}},
\code{\link{multmixEM}},
\code{\link{mvnormalmixEM}},
\code{\link{poisregmixEM}},
\code{\link{regmixEM}},
\code{\link{regmixEM.lambda}},
\code{\link{regmixEM.loc}},
\code{\link{regmixEM.mixed}},
\code{\link{regmixEM.chgpt}},
\code{\link{repnormmixEM}},
\code{\link{expRMM_EM}},
\code{\link{weibullRMM_SEM}}
}
\examples{
data(faithful)
attach(faithful)
set.seed(100)
out <- normalmixEM(waiting, mu=c(50,80), sigma=c(5,5), lambda=c(.5,.5))
summary(out)
}
\keyword{file}

