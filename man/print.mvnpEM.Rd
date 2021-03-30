\name{print.mvnpEM}
\alias{print.mvnpEM}
\title{Printing of Results from the mvnpEM Algorithm Output}
\usage{
\method{print}{mvnpEM}(x, \dots)
}
\arguments{
  \item{x}{an object of class \code{mvnpEM} such as a result of a call
  to \code{\link{mvnpEM}}}
  \item{\dots}{Additional arguments to \code{\link{print}}}
}
\description{
\code{\link[base]{print}} method for class \code{mvnpEM}.
}
\details{
 \code{print.mvnpEM} prints the elements of an \code{mvnpEM} object
 without printing the data or the posterior probabilities.
 (These may still be accessed as \code{x$data} and \code{x$posteriors}.)
}
\value{
  \code{print.mvnpEM} returns (invisibly) the full value of \code{x} itself,
  including the \code{data} and \code{posteriors} elements.
}
\seealso{ 
\code{\link{mvnpEM}},
\code{\link{plot.mvnpEM}}
\code{\link{summary.mvnpEM}}
}
\examples{
# Example as in Chauveau and Hoang (2015) with 6 coordinates
\dontrun{
m=2; r=6; blockid <-c(1,1,2,2,3,3) # 3 bivariate blocks 
# generate some data x ...
a <- mvnpEM(x, mu0=2, blockid, samebw=F) # adaptive bandwidth
print(a)}
}
\keyword{file}

