\name{summary.npEM}
\alias{summary.npEM}
\alias{print.summary.npEM}
\title{Summarizing non- and semi-parametric multivariate mixture model fits}
\usage{
\method{summary}{npEM}(object, \dots)
\method{print}{summary.npEM}(x, digits=3, \dots)
}
\arguments{
  \item{object,x}{an object of class \code{npEM} such as a result of a call
  to \code{\link{npEM}}}
  \item{digits}{Significant digits for printing values}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
\code{\link[base]{summary}} method for class \code{npEM}.
}
\details{
 \code{\link{summary.npEM}} prints means and variances of each block for
 each component.  These quantities might not be part of the model, but they
 are estimated nonparametrically based on the posterior probabilities and the
 data.
}
\value{
  The function \code{\link{summary.npEM}} returns a list of type \code{summary.npEM}
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
}
\seealso{ 
\code{\link{npEM}},
\code{\link{plot.npEM}}
}
\examples{
data(Waterdata)
set.seed(100)
\dontrun{
a <- npEM(Waterdata[,3:10], 3, bw=4) # Assume indep but not iid
summary(a) 

b <- npEM(Waterdata[,3:10], 3, bw=4, blockid=rep(1,8)) # Now assume iid
summary(b)
}
}
\keyword{file}

