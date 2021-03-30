\name{print.npEM}
\alias{print.npEM}
\title{Printing non- and semi-parametric multivariate mixture model fits}
\usage{
\method{print}{npEM}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{npEM} such as a result of a call
  to \code{\link{npEM}}}
  \item{\dots}{Additional arguments to \code{\link{print}}}
}
\description{
\code{\link[base]{print}} method for class \code{npEM}.
}
\details{
 \code{print.npEM} prints the elements of an \code{npEM} object
 without printing the data or the posterior probabilities.
 (These may still be accessed as \code{x$data} and \code{x$posteriors}.)
}
\value{
  \code{print.npEM} returns (invisibly) the full value of \code{x} itself,
  including the \code{data} and \code{posteriors} elements.
}
\seealso{ 
\code{\link{npEM}},
\code{\link{plot.npEM}}
\code{\link{summary.npEM}}
}
\examples{
data(Waterdata)
set.seed(100)
\dontrun{npEM(Waterdata[,3:10], 3, bw=4, verb=FALSE) # Assume indep but not iid}
}
\keyword{file}

