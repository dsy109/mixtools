\name{lambda.pert}
\title{Perturbation of Mixing Proportions}
\alias{lambda.pert}
\usage{
lambda.pert(lambda, pert)
}
\description{
Perturbs a set of mixing proportions by first scaling the
mixing proportions, then taking the logit of the scaled values,
perturbing them, and inverting back to produce a set of
new mixing proportions.
}
\arguments{
  \item{lambda}{A vector of length k giving the mixing proportions which
  are to be perturbed.}
  \item{pert}{A vector (possibly of length k-1) for which to perturb \code{lambda}.
  If the length is less than k, then values of the vector are recycled.  If length
  is greater than k, then only the first k values are used.}
}
\value{
  \code{lambda.pert} returns new \code{lambda} values perturbed by \code{pert}.  
}
\details{
This function is called by \code{regmixMH}.
}
\seealso{
\code{\link{regmixMH}}
}
\examples{
set.seed(100)
x <- c(0.5, 0.2, 0.3)
lambda.pert(x, rcauchy(3))

}

\keyword{internal}
