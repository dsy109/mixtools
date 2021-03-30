\name{depth}
\alias{depth}
\title{Elliptical and Spherical Depth}
\description{
  Computation of spherical or elliptical depth.
}
\usage{
 depth(pts, x, Cx = var(x))

}
\arguments{
  \item{pts}{A kxd matrix containing the k points that one wants to compute the depth. Each row is a point. }
  \item{x}{A nxd matrix containing the reference data. Each row is an observation.}
  \item{Cx}{A dxd scatter matrix for the data x where the default is var(x). When Cx = I(d), it returns the sphercial depth.}

}
\value{
  \code{depth} returns a k-vector where each entry is the elliptical depth of a point in \code{pts}.

}
\references{
  Elmore, R. T., Hettmansperger, T. P. and Xuan, F. (2000) Spherical Data Depth and a Multivariate Median,
  \emph{Proceedings of Data Depth: Robust Multivariate Statistical Analysis, Computational Geometry and Applications}.
}
\seealso{\code{\link{regcr}}
}
\examples{
  set.seed(100)
  x <- matrix(rnorm(200),nc = 2)
  depth(x[1:3, ], x)
}
\note{\code{depth} is used in \code{regcr}.}

\keyword{file}
