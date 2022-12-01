\name{boot.se}
\title{Performs Parametric Bootstrap for Standard Error 
Approximation}
\alias{boot.se}
\usage{
boot.se(em.fit, B = 100, arbmean = TRUE, arbvar = TRUE, 
        N = NULL, ...)
}

\description{
  Performs a parametric bootstrap by producing B bootstrap samples for the parameters in the specified
  mixture model.
}
\arguments{
  \item{em.fit}{An object of class \code{mixEM}.  The estimates produced in \code{em.fit} will be used as the 
    parameters for the distribution from which we generate the bootstrap data.}
  \item{B}{The number of bootstrap samples to produce.  The default is 100,
    but ideally, values of 1000 or more would be more acceptable.}
  \item{arbmean}{If FALSE, then a scale mixture analysis can be performed for \code{mvnormalmix}, \code{normalmix}, \code{regmix},
    or \code{repnormmix}. The default is TRUE.}
  \item{arbvar}{If FALSE, then a location mixture analysis can be performed for \code{mvnormalmix}, \code{normalmix}, \code{regmix},
    or \code{repnormmix}. The default is TRUE.}
  \item{N}{An n-vector of number of trials for the logistic regression type \code{logisregmix}.  
  If NULL, then \code{N} is an n-vector of 1s for binary logistic regression.}
  \item{...}{Additional arguments passed to the various EM algorithms for the mixture of interest.} 
}
\value{
  \code{boot.se} returns a list with the bootstrap samples and standard errors for the mixture of interest.
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
}
\examples{
## Bootstrapping standard errors for a regression mixture case.

data(NOdata)
attach(NOdata)
set.seed(100)
em.out <- regmixEM(Equivalence, NO, arbvar = FALSE)
out.bs <- boot.se(em.out, B = 10, arbvar = FALSE)
out.bs

}


\keyword{file}
