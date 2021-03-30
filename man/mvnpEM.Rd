\name{mvnpEM}
\title{EM-like Algorithm for Nonparametric Mixture Models 
with Conditionally Independent Multivariate Component Densities}
\alias{mvnpEM}
\usage{
mvnpEM(x, mu0, blockid = 1:ncol(x), samebw = TRUE, 
       bwdefault = apply(x,2,bw.nrd0), init = NULL,
       eps = 1e-8, maxiter = 500, verb = TRUE)
}
                    
                    
                    
                    

\description{
  An extension of the original \code{\link{npEM}} algorithm, for mixtures
  of multivariate data where the coordinates of a row (case)
  in the data matrix are assumed to be made of independent but multivariate blocks (instead of just coordinates), 
  conditional on the mixture
  component (subpopulation) from which they are drawn (Chauveau and Hoang 2015).
}
\arguments{
  \item{x}{An \eqn{n\times r}{n x r} matrix of data.  Each of the \eqn{n} rows is a case,
  and each case has \eqn{r} repeated measurements.  These measurements are assumed
  to be conditionally independent, conditional on the mixture component (subpopulation) 
  from which the case is drawn.}
  \item{mu0}{Either an \eqn{m\times r}{m x r} matrix specifying the initial
      centers for the \link{kmeans} function, or an integer \eqn{m} specifying the 
      number of initial centers, which are then chosen randomly in \link{kmeans}}
  \item{blockid}{A vector of length \eqn{r} identifying coordinates 
  (columns of \code{x}) that are in the same block.  The default has all distinct elements, 
  indicating that the model has \eqn{r} blocks of dimension 1, in which case the model is handled 
  directly by the \code{\link{npEM}} algorithm.
  See example below for actual multivariate blocks example.}
  
  \item{samebw}{Logical:  If \code{TRUE}, use the same bandwidth per coordinate for
  all iteration and all components.  If \code{FALSE}, 
  use a separate bandwidth for each component and coordinate, and update
  this bandwidth at each iteration of the algorithm using a suitably
  modified \code{\link{bw.nrd0}} method as described in 
  Benaglia et al (2011) and Chauveau and Hoang (2015).}


  \item{bwdefault}{Bandwidth default for density estimation,a simplistic application of the 
  default \code{\link{bw.nrd0}} for each coordinate (column) of the data.} 
  
  \item{init}{Initialization method, based on an initial \eqn{n\times m}{n x m}
  matrix for the posterior probabilities. If \code{NULL},
  a  \code{\link{kmeans}} clustering with \code{mu0} initial centers is applied to the data and the 
  initial matrix of posteriors is built from the result.}
  
  \item{eps}{Tolerance limit for declaring algorithm convergence.  Convergence
  is declared whenever the maximum change in any coordinate of the 
  \code{lambda} vector (of mixing proportion estimates) does not exceed 
  \code{eps}.}

  \item{maxiter}{The maximum number of iterations allowed; convergence
  may be declared before \code{maxiter} iterations (see \code{eps} above).}
  
  \item{verb}{Verbose mode; if TRUE, print updates for every iteration of the algorithm as
  it runs}
}
\value{
\code{mvnpEM} returns a list of class \code{mvnpEM} with the following items:
  \item{data}{The raw data (an \eqn{n\times r}{n x r} matrix).}
  
  \item{posteriors}{An \eqn{n\times m}{n x m} matrix of posterior probabilities for
   each observation (row).}
   
  \item{lambda}{The sequence of mixing proportions over iterations.}
  
  \item{blockid}{The \code{blockid} input argument. Needed by any method that produces density 
  estimates from the output, like \code{\link{plot.mvnpEM}}.}


  \item{samebw}{The \code{samebw} input argument. 
  Needed by any method that produces density estimates from the output, like \code{\link{plot.mvnpEM}}.}

  \item{bandwidth}{The final bandwidth matrix 
  after convergence of the algorithm.
  Its shape depends on the \code{samebw} input argument. If \code{samebw = TRUE}, a 
  vectors with the bandwidth value for each of the \code{r} coordinates (same for all components and iterations).
  If \code{samebw = FALSE}, a \eqn{m\times r}{m x r} matrix, where each row is associated to one component 
  and gives the \eqn{r} bandwidth values, one for each coordinate.
  Needed by any method that produces density estimates from the output, like \code{\link{plot.mvnpEM}}.}

  \item{lambdahat}{The final mixing proportions.}
  \item{loglik}{The sequence of pseudo log-likelihood values over iterations.}
}
\seealso{
\code{\link{plot.mvnpEM}}, \code{\link{npEM}}
}
\references{
   \itemize{
   \item Benaglia, T., Chauveau, D., and Hunter, D. R. (2009), An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   Journal of Computational and Graphical Statistics, 18, 505-526.
   
   \item Benaglia, T., Chauveau, D. and Hunter, D.R. (2011),
   Bandwidth Selection in an EM-like algorithm for nonparametric multivariate mixtures.
  Nonparametric Statistics and Mixture Models: A Festschrift in
  	Honor of Thomas P. Hettmansperger. World Scientific Publishing Co., 
  	pages 15-27.

   \item  Chauveau, D., and Hoang, V. T. L. (2015),
   Nonparametric mixture models with conditionally independent multivariate component densities, 
   Preprint under revision.
   \url{https://hal.archives-ouvertes.fr/hal-01094837}
   }
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
