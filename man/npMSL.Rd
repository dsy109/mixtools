\name{npMSL}
\title{Nonparametric EM-like Algorithm for Mixtures of Independent Repeated Measurements - Maximum Smoothed Likelihood version}
\alias{npMSL}
\usage{
npMSL(x, mu0, blockid = 1:ncol(x), 
      bw = bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE, 
      bwmethod = "S", h = bw, eps = 1e-8, 
      maxiter=500, bwiter = maxiter, nbfold = NULL,
      ngrid=200, post=NULL, verb = TRUE)
}

\description{
  Returns nonparametric Smoothed Likelihood algorithm output 
  (Levine et al, 2011) for mixtures
  of multivariate (repeated measures) data where the coordinates of a row (case)
  in the data matrix are assumed to be independent, conditional on the mixture
  component (subpopulation) from which they are drawn.
}
\arguments{
  \item{x}{An \eqn{n\times r}{n x r} matrix of data.  Each of the \eqn{n} rows is a case,
  and each case has \eqn{r} repeated measurements.  These measurements are assumed
  to be conditionally independent, conditional on the mixture component (subpopulation) 
  from which the case is drawn.}
  \item{mu0}{Either an \eqn{m\times r}{m x r} matrix specifying the initial
      centers for the \link{kmeans} function, or an integer \eqn{m} specifying the 
      number of initial centers, which are then choosen randomly in
  \link{kmeans}}
  \item{blockid}{A vector of length \eqn{r} identifying coordinates 
  (columns of \code{x}) that are
  assumed to be identically distributed (i.e., in the same block).  
  For instance,
  the default has all distinct elements, indicating that no two coordinates 
  are assumed identically distributed and thus a separate set of \eqn{m} 
  density estimates is produced for each column of \eqn{x}.  On the other hand,
  if \code{blockid=rep(1,ncol(x))}, then the coordinates in each row 
  are assumed conditionally i.i.d.}
  
  \item{bw}{Bandwidth for density estimation, equal to the standard deviation 
  of the kernel density.  By default, a simplistic application of the 
  default \code{\link{bw.nrd0}} 
  bandwidth used by \code{\link{density}} to the entire dataset.}
  
  \item{samebw}{Logical:  If \code{TRUE}, use the same bandwidth for
  each iteration and for each component and block.  If \code{FALSE}, 
  use a separate bandwidth for each component and block, and update
  this bandwidth at each iteration of the algorithm 
  until \code{bwiter} is reached (see below). Two adaptation methods are provided, 
  see \code{bwmethod} below.}
  
  \item{bwmethod}{Define the adaptive bandwidth strategy when \code{samebw = FALSE}, in which case
  the bandwidth depends on each component, block, and iteration of the algorithm.
  If set to "S" (the default), adaptation is done using a suitably
  modified \code{\link{bw.nrd0}} method as described in 
  Benaglia et al (2011).
  If set to "CV", an adaptive \eqn{k}-fold Cross Validation method is applied,
  as described in Chauveau et al (2014), where \code{nbfold} is the number of subsamples.
  This corresponds to a Leave-\eqn{[n/nbfold]}-Out CV.
  }
  \item{h}{Alternative way to specify the bandwidth, to provide backward 
  compatibility.}
  \item{eps}{Tolerance limit for declaring algorithm convergence.  Convergence
  is declared whenever the maximum change in any coordinate of the 
  \code{lambda} vector (of mixing proportion estimates) does not exceed 
  \code{eps}.}

  \item{maxiter}{The maximum number of iterations allowed, convergence
  may be declared before \code{maxiter} iterations (see \code{eps} above).}
  
  \item{bwiter}{The maximum number of iterations allowed for adaptive bandwidth stage,
  when \code{samebw = FALSE}. If set to \code{0}, then the initial bandwidth matrix is used without adaptation.}

  \item{nbfold}{A parameter passed to the internal function \code{wbs.kCV}, which controls the weighted bandwidth selection by k-fold cross-validation.}

  
  \item{ngrid}{Number of points in the discretization of the intervals over which are approximated the (univariate) integrals for non linear smoothing of the log-densities, as required in the E step of the npMSL algorithm,
  see Levine et al (2011).}

  \item{post}{If non-NULL, an \eqn{n\times m}{n x m} matrix specifying the initial 
  posterior probability vectors for each of the observations, i.e., the
  initial values to start the EM-like algorithm.}
  \item{verb}{If TRUE, print updates for every iteration of the algorithm as
  it runs}
}
\value{
\code{npMSL} returns a list of class \code{npEM} with the following items:
  \item{data}{The raw data (an \eqn{n\times r}{n x r} matrix).}
  \item{posteriors}{An \eqn{n\times m}{n x m} matrix 
   of posterior probabilities for
   observation.}
  \item{bandwidth}{If \code{samebw==TRUE}, 
  same as the \code{bw} input argument; otherwise, value of \code{bw} matrix
  at final iteration.  This
  information is needed by any method that produces density estimates from the
  output.}
  \item{blockid}{Same as the \code{blockid} input argument, but recoded to have
  positive integer values.  Also needed by any method that produces density 
  estimates from the output.}
  \item{lambda}{The sequence of mixing proportions over iterations.}
  \item{lambdahat}{The final mixing proportions.}
  \item{loglik}{The sequence of log-likelihoods over iterations.}
  \item{f}{An array of size \eqn{ngrid \times m \times l}{ngrid x m x l},
  returning last values of density for component \eqn{j}
  and block \eqn{k} over \code{grid} points.}
  \item{meanNaN}{Average number of \code{NaN} that occured over iterations (for internal testing and control purpose).}
  \item{meanUdfl}{Average number of "underflow" that occured over iterations (for internal testing and control purpose).}
}
\seealso{
\code{\link{npEM}}, \code{\link{plot.npEM}}, 
\code{\link{normmixrm.sim}}, \code{\link{spEMsymloc}},
\code{\link{spEM}}, \code{\link{plotseq.npEM}}
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

   
   \item Chauveau D., Hunter D. R. and Levine M. (2014),
  Semi-Parametric Estimation for Conditional Independence Multivariate Finite Mixture Models.
  Preprint (under revision).
    
  \item Levine, M., Hunter, D. and Chauveau, D. (2011),
   Maximum Smoothed Likelihood for Multivariate Mixtures,
   Biometrika 98(2): 403-416.
   }
}

\examples{
## Examine and plot water-level task data set.
## Block structure pairing clock angles that are directly opposite one
## another (1:00 with 7:00, 2:00 with 8:00, etc.)
set.seed(111) # Ensure that results are exactly reproducible
data(Waterdata)
blockid <- c(4,3,2,1,3,4,1,2) # see Benaglia et al (2009a)

\dontrun{
a <- npEM(Waterdata[,3:10], mu0=3, blockid=blockid, bw=4)  # npEM solution
b <- npMSL(Waterdata[,3:10], mu0=3, blockid=blockid, bw=4) # smoothed version

# Comparisons on the 4 default plots, one for each block
par(mfrow=c(2,2))
for (l in 1:4){
plot(a, blocks=l, breaks=5*(0:37)-92.5,
	xlim=c(-90,90), xaxt="n",ylim=c(0,.035), xlab="")
plot(b, blocks=l, hist=FALSE, newplot=FALSE, addlegend=FALSE, lty=2,
	dens.col=1)
axis(1, at=30*(1:7)-120, cex.axis=1)
legend("topleft",c("npMSL"),lty=2, lwd=2)}
}

}


\keyword{file}
