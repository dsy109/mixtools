\name{wquantile}
\alias{wIQR}
\alias{wquantile}
\title{Weighted quantiles}
\description{
  Functions to compute weighted quantiles and the weighted interquartile range.
  }
\usage{
wquantile(wt = rep(1,length(x)), x, probs, already.sorted = FALSE, 
          already.normalized = FALSE)
wIQR(wt = rep(1,length(x)), x, already.sorted = FALSE, 
     already.normalized = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{wt}{ Vector of weights }
  \item{x}{ Vector of data, same length as \code{wt} }
  \item{probs}{ Numeric vector of probabilities with values in [0,1].}
  \item{already.sorted}{ If FALSE, sort \code{wt} and \code{x} in increasing
  order of \code{x}.  If TRUE, it is assumed that \code{wt} and \code{x}
  are already sorted.}
  \item{already.normalized}{ If FALSE, normalize \code{wt} by diving each
  entry by the sum of all entries.  If TRUE, it is assumed that \code{sum(wt)==1}}
}
\details{
  \code{wquantile} uses the \code{\link{findInterval}} function.  \code{wIQR}
  calls the \code{wquantile} function.
  }
\value{
  Returns the sample quantiles or interquartile range of a discrete distribution with
  support points \code{x} and corresponding probability masses \code{wt}
}
\seealso{ \code{\link{npEM}}}
\examples{
IQR(1:10)
wIQR(x=1:10) # Note:  Different algorithm than IQR function
wIQR(1:10,1:10) # Weighted quartiles are now 4 and 8
}

\keyword{robust}
