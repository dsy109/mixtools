\name{spEMsymloc}
\title{Semiparametric EM-like Algorithm for univariate symmetric location mixture}
\alias{spEMsymloc}
\usage{
spEMsymloc(x, mu0, bw = bw.nrd0(x), h=bw, eps = 1e-8, maxiter = 100, 
           stochastic = FALSE, verbose = FALSE)
}

\description{
  Returns semiparametric EM algorithm output (Bordes et al, 2007,
  and Benaglia et al, 2009)  for location mixtures
  of univariate data and symmetric component density.
}
\arguments{
  \item{x}{A vector of length \eqn{n}{n} consisting of the data.}
  \item{mu0}{Either a vector specifying the initial
      centers for the \link{kmeans} function, and from which the number of 
   component is obtained,  or an integer \eqn{m} specifying the 
      number of initial centers, which are then choosen randomly in
   \link{kmeans}.}

  \item{bw}{Bandwidth for density estimation, equal to the standard deviation 
  of the kernel density.}
  \item{h}{Alternative way to specify the bandwidth, to provide backward 
  compatibility.}

  \item{eps}{Tolerance limit for declaring algorithm convergence.  Convergence
  is declared before \code{maxiter} iterations
  whenever the maximum change in any coordinate of the 
  \code{lambda} (mixing proportion estimates) and \code{mu} (means) 
  vector  does not exceed \code{eps}.}

  \item{maxiter}{The maximum number of iterations allowed, for both 
  stochastic and non-stochastic versions; 
  for non-stochastic algorithms (\code{stochastic = FALSE}), convergence
  may be declared before \code{maxiter} iterations (see \code{eps} above).}

  \item{stochastic}{Flag, if FALSE (the default), runs the non-stochastic version
      of the algorithm, as in Benaglia et al (2009). Set to TRUE to
      run a stochastic version which simulates the posteriors at each
      iteration (as in Bordes et al, 2007),
   and runs for \code{maxiter} iterations.}


  \item{verbose}{If TRUE, print updates for every iteration of the algorithm as
  it runs}
}
\value{
\code{spEMsymloc} returns a list of class \code{npEM} with the following items:
  \item{data}{The raw data (an \eqn{n\times r}{n x r} matrix).}
  \item{posteriors}{An \eqn{n\times m}{n x m} matrix of posterior probabilities for
   observations. If \code{stochastic = TRUE}, this matrix is computed 
   from an average over the \code{maxiter} iterations.}
  \item{bandwidth}{Same as the \code{bw} input argument, returned because this
  information is needed by any method that produces density estimates from the
  output.}
  \item{lambda}{The sequence of mixing proportions over iterations.}
  \item{lambdahat}{The final estimate for mixing proportions if \code{stochastic = 
  FALSE}, the average over the sequence if \code{stochastic = TRUE}.}
  \item{mu}{the sequence of component means over iterations.}
  \item{muhat}{the final estimate of component means if \code{stochastic = 
   FALSE}, the average over the sequence if \code{stochastic = TRUE}.}
  \item{symmetric}{Flag indicating that the kernel density estimate is
   using a symmetry assumption.}
}
\seealso{
\code{\link{plot.npEM}}, \code{\link{rnormmix}},
\code{\link{npEM}}, \code{\link{spEMsymlocN01}}, \code{\link{plotseq.npEM}}
}
\references{
   \itemize{
   \item Benaglia, T., Chauveau, D., and Hunter, D. R., An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   Journal of Computational and Graphical Statistics, 18, 505-526, 2009.
   
   \item Benaglia, T., Chauveau, D., Hunter, D. R., and Young, D.
	mixtools: An R package for analyzing finite mixture models.
	Journal of Statistical Software, 32(6):1-29, 2009.


   \item Bordes, L., Chauveau, D., and Vandekerkhove, P. (2007),
   An EM algorithm for a semiparametric mixture model, 
   Computational Statistics and Data Analysis, 51: 5429-5443.
   }
}
\examples{
## Example from a normal location mixture
set.seed(100)
n <- 200
lambda <- c(1/3,2/3)
mu <- c(0, 4); sigma<-rep(1, 2)
x <- rnormmix(n, lambda, mu, sigma)
out.stoc <- spEMsymloc(x, mu0=c(-1, 2), stochastic=TRUE)
out.nonstoc <- spEMsymloc(x, mu0=c(-1, 2))
}

\keyword{file}
