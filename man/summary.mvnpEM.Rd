\name{summary.mvnpEM}
\alias{summary.mvnpEM}
\alias{print.summary.mvnpEM}
\title{Summarizing Fits for Nonparametric Mixture Models 
with Conditionally Independent Multivariate Component Densities}
\usage{
\method{summary}{mvnpEM}(object, \dots)
\method{print}{summary.mvnpEM}(x, digits=3, \dots)
}
\arguments{
  \item{object,x}{an object of class \code{mvnpEM} such as a result of a call
  to \code{\link{mvnpEM}}}
  \item{digits}{Significant digits for printing values}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
\code{\link[base]{summary}} method for class \code{mvnpEM}.
}
\details{
 \code{\link{summary.mvnpEM}} prints means and variances of each block for
 each component.  These quantities might not be part of the model, but they
 are estimated nonparametrically based on the posterior probabilities and the
 data.
}
\value{
  The function \code{\link{summary.mvnpEM}} returns a list of type \code{summary.mvnpEM}
  with the following components:
  \item{n}{The number of observations}
  \item{m}{The number of mixture components}
  \item{B}{The number of blocks}
  \item{blockid}{The block ID (from 1 through B) for each of the coordinates
  of the multivariate observations.  The \code{blockid} component is of length
  \eqn{r}, the dimension of each observation.}
  \item{means}{A \eqn{B\times m}{B x m} matrix giving the estimated mean of each
  block in each component.}
  \item{variances}{Same as \code{means} but giving the estimated variances instead.}
}
\references{
 Benaglia, T., Chauveau, D., and Hunter, D. R. (2009), An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   \emph{Journal of Computational and Graphical Statistics}, \bold{18(2)}, 505--526.
   
Chauveau, D., and Hoang, V. T. L. (2015),
   Nonparametric mixture models with conditionally independent multivariate component densities, 
   Preprint under revision.
   \url{https://hal.archives-ouvertes.fr/hal-01094837}
}
\seealso{ 
\code{\link{mvnpEM}},
\code{\link{plot.mvnpEM}}
}
\examples{
# Example as in Chauveau and Hoang (2015) with 6 coordinates
\dontrun{
m=2; r=6; blockid <-c(1,1,2,2,3,3) # 3 bivariate blocks 
# generate some data x ...
a <- mvnpEM(x, mu0=2, blockid, samebw=F) # adaptive bandwidth
plot(a) # this S3 method produces 6 plots of univariate marginals
summary(a)}
}
\keyword{file}

