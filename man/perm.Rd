\name{perm}
\title{Permutation Function}
\alias{perm}
\usage{
perm(n, r, v = 1:n)
}
\description{
Enumerates the possible permutations of a specified size from the elements of a vector having the same size.
}
\arguments{
  \item{n}{Size of the source vector.}
  \item{r}{Size of the target vector.}
  \item{v}{Source vector.  Must be a vector of length \code{n}.}
}
\value{
  \code{perm} returns a matrix where each row contains one of the permutations of length \code{r}.
}
\details{
This function is called by \code{segregmixEM} and the associated internal functions.  This is also a simplified version of the function \code{permutations} found in the package \code{gtools}.
}
\seealso{
\code{\link{segregmixEM}}
}
\examples{
perm(3, 3, 2:4)
}


\keyword{internal}
