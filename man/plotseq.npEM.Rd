\name{plotseq.npEM}
\alias{plotseq.npEM}
\title{Plotting sequences of estimates from non- or semiparametric EM-like Algorithm}
\usage{
  \method{plotseq}{npEM}(x, \dots)
}
\description{
  Returns plots of the sequences of scalar parameter
  estimates along iterations from an object of class \code{npEM}.
}
\arguments{
  \item{x}{an object of class \code{npEM}, as output by \code{\link{npEM}} 
  or \code{\link{spEMsymloc}}}
  \item{\dots}{further parameters that are passed to \code{\link{plot}}}
}
\details{
\code{plotseq.npEM} returns a figure with one plot for each component 
 proportion, and, in the case of \code{\link{spEMsymloc}}, one plot for each 
component mean.
}
\seealso{
\code{\link{plot.npEM}}, \code{\link{rnormmix}},
\code{\link{npEM}}, \code{\link{spEMsymloc}}
}
\references{
   \itemize{
   \item Benaglia, T., Chauveau, D., and Hunter, D. R. (2009), An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   Journal of Computational and Graphical Statistics (to appear).

   \item Bordes, L., Chauveau, D., and Vandekerkhove, P. (2007),
   An EM algorithm for a semiparametric mixture model, 
   Computational Statistics and Data Analysis, 51: 5429-5443.   
   }
}
\examples{
## Example from a normal location mixture
n <- 200
set.seed(100)
lambda <- c(1/3,2/3)
mu <- c(0, 4); sigma<-rep(1, 2)
x <- rnormmix(n, lambda, mu, sigma)
b <- spEMsymloc(x, mu0=c(-1, 2), stochastic=FALSE)
plotseq(b)
bst <- spEMsymloc(x, mu0=c(-1, 2), stochastic=TRUE)
plotseq(bst)
}

\keyword{file}
