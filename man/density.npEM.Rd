\name{density.npEM}
\title{Normal kernel density estimate for nonparametric EM output}
\alias{density.npEM} 
\usage{
\method{density}{npEM}(x, u=NULL, component=1, block=1, scale=FALSE, \dots)
}

\description{
    Takes an object of class \code{npEM} and returns an object of class
    \code{\link{density}} giving the kernel density estimate for the selected
    component and, if applicable, the selected block.
} 
\arguments{
  \item{x}{An object of class \code{npEM} such as the output
  of the \code{\link{npEM}} or \code{\link{spEMsymloc}} functions.}
  \item{u}{Vector of points at which the density is to be evaluated}
  \item{component}{Mixture component number; should be an integer from 1 to the
  number of columns of \code{x$posteriors}.}
  \item{block}{Block of repeated measures.  Only applicable in repeated measures
  case, for which \code{x$blockid} exists; should be an integer from 1 to
  \code{max(x$blockid)}.}
  \item{scale}{Logical:  If TRUE, multiply the density values by the 
  corresponding mixing proportions found in \code{x$lambdahat}}
  \item{\dots}{Additional arguments; not used by this method.}
}
\details{
  The bandwidth is taken to be the same as that used to produce the \code{npEM}
  object, which is given by \code{x$bandwidth}.
}
\value{
  \code{density.npEM} returns a list of type \code{"density"}.  See 
  \code{\link{density}} for details.  In particular, the output of 
  \code{density.npEM} may be used directly by functions such as 
  \code{\link{plot}} or \code{\link{lines}}.
}

\seealso{ 
\code{\link{npEM}}, \code{\link{spEMsymloc}}, \code{\link{plot.npEM}}
} 

\examples{ 
## Look at histogram of Old Faithful waiting times
data(faithful)
Minutes <- faithful$waiting
hist(Minutes, freq=FALSE)

## Superimpose equal-variance normal mixture fit:
set.seed(100)
nm <- normalmixEM(Minutes, mu=c(50,80), sigma=5, arbvar=FALSE, fast=TRUE)
x <- seq(min(Minutes), max(Minutes), len=200)
for (j in 1:2) 
  lines(x, nm$lambda[j]*dnorm(x, mean=nm$mu[j], sd=nm$sigma), lwd=3, lty=2)
  
## Superimpose several semiparametric fits with different bandwidths:
bw <- c(1, 3, 5)
for (i in 1:3) {
  sp <- spEMsymloc(Minutes, c(50,80), bw=bw[i], eps=1e-3)
  for (j in 1:2) 
    lines(density(sp, component=j, scale=TRUE), col=1+i, lwd=2)    
}
legend("topleft", legend=paste("Bandwidth =",bw), fill=2:4)
}
\keyword{file}
