\name{density.spEM}
\title{Normal kernel density estimate for semiparametric EM output}
\alias{density.spEM} 
\usage{
\method{density}{spEM}(x, u=NULL, component=1, block=1, scale=FALSE, \dots)
}

\description{
    Takes an object of class \code{spEM} and returns an object of class
    \code{\link{density}} giving the kernel density estimate.
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
  \code{density.spEM} returns a list of type \code{"density"}.  See 
  \code{\link{density}} for details.  In particular, the output of 
  \code{density.spEM} may be used directly by functions such as 
  \code{\link{plot}} or \code{\link{lines}}.
}

\seealso{ 
\code{\link{spEM}}, \code{\link{spEMsymloc}}, \code{\link{plot.spEM}}
} 

\examples{ 
set.seed(100)
mu <- matrix(c(0, 15), 2, 3)
sigma <- matrix(c(1, 5), 2, 3)
x <- rmvnormmix(300, lambda = c(.4,.6), mu = mu, sigma = sigma)

d <- spEM(x, mu0 = 2, blockid = rep(1,3), constbw = TRUE) 
plot(d, xlim=c(-10, 40), ylim = c(0, .16), xlab = "", breaks = 30, 
     cex.lab=1.5, cex.axis=1.5) # plot.spEM calls density.spEM here
}
\keyword{file}
