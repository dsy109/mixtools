\name{spRMM_SEM}
\title{Stochastic EM algorithm for semiparametric scaled mixture of censored data}
\alias{spRMM_SEM}
\usage{
  spRMM_SEM(t, d = NULL, lambda = NULL, scaling = NULL, 
          centers = 2, kernelft = triang_wkde, 
          bw = rep(bw.nrd0(t),length(t)), averaged = TRUE,
          epsilon = 1e-08, maxit = 100, batchsize = 1, verb = FALSE) }

\description{
  Stochastic EM algorithm for semiparametric scaled mixture for randomly right censored data.
}

\arguments{
  \item{t}{A vector of \eqn{n}{n} real positive lifetime (possibly censored) durations.
  If \code{d} is not \code{NULL} then a vector of random censoring times
  \code{c} occurred,  
  so that \eqn{x= min(x,c)} and \eqn{d = I(x <= c)}.}
  
  \item{d}{The vector of censoring indication, where 1 means observed lifetime data, 
  and 0 means censored lifetime data.}
  
  \item{lambda}{Initial value of mixing proportions.
  If \code{NULL}, then \code{lambda} is set to \code{rep(1/k,k)}.}
  
  \item{scaling}{Initial value of scaling between components, 
        set to 1 if \code{NULL}.}

  \item{centers}{initial centers for initial call to kmeans for initialization.}

  \item{kernelft}{.}
  
  \item{bw}{Bandwidth in the kernel hazard estimates.}
  
  \item{averaged}{averaged.}
  
  \item{epsilon}{Tolerance limit.}

  \item{maxit}{The number of iterations allowed.}
  
  \item{batchsize}{The batchsize (see reference below).}

  \item{verb}{If TRUE, print updates for every iteration of the algorithm as
  it runs}
}

\value{
\code{spRMM_SEM} returns a list of class \code{"spRMM"} with the following items:
  \item{t}{The input data.}
  \item{d}{The input censoring indicator.}
  \item{lambda}{The estimates for the mixing proportions.}
  \item{scaling}{The estimates for the components scaling.}
  \item{posterior}{An \eqn{n\times k}{n x 2} matrix of posterior probabilities for
   observation, after convergence of the algorithm.}
  
  \item{loglik}{The (pseudo) log-likelihood value at convergence of the algorithm.}

  \item{all.loglik}{The sequence of log-likelihood values over iterations.}

  \item{all.lambda}{The sequence of mixing proportions over iterations.}
  
  \item{all.scaling}{The sequence of scaling parameter over iterations.}
  
  \item{meanpost}{Posterior probabilities averaged over iterations.}
  
  \item{survival}{Kaplan-Meier last iteration estimate (a \code{stepfun} object).}
  
  \item{hazard}{Hazard rate last iteration estimate evaluated at \code{final.t}.}
  
  \item{final.t}{Last iteration unscaled sample (see reference).}
  
  \item{s.hat}{Kaplan-Meier average estimate.}
  
  \item{t.hat}{Ordered unscaled sample, for testing purpose.}
  
  \item{avg.od}{For testing purpose only.}
  
  \item{hazard.hat}{Hazard rate average estimate on \code{t.hat}.}
  
  \item{batch.t}{Batch sample (not ordered), see reference.}
  
  \item{batch.d}{Associated event indicators just \code{rep(d,batchsize)}, for testing purpose.}
  
  \item{sumNaNs}{Internal control of numerical stability.}

   \item{ft}{A character vector giving the name of the function.}
}
\seealso{
Related functions: 
  \code{\link{plotspRMM}},
  \code{\link{summary.spRMM}}.
  
  Other models and algorithms for censored lifetime data
  (name convention is model_algorithm):
  \code{\link{expRMM_EM}},
  \code{\link{weibullRMM_SEM}}.
}
\references{
   \itemize{
   \item Bordes, L., and Chauveau, D. (2016),
   Stochastic EM algorithms for parametric and semiparametric mixture models 
   for right-censored lifetime data, 
   Computational Statistics, Volume 31, Issue 4, pages 1513-1538.
   \url{https://link.springer.com/article/10.1007/s00180-016-0661-7}
   }
}
\author{Didier Chauveau}

\examples{
\dontrun{
n=500 # sample size
m=2 # nb components
lambda=c(0.4, 0.6) # parameters
meanlog=3; sdlog=0.5; scale=0.1
set.seed(12)
# simulate a scaled mixture of lognormals
x <- rlnormscalemix(n, lambda, meanlog, sdlog, scale)
cs=runif(n,20,max(x)+400) # Censoring (uniform) and incomplete data
t <- apply(cbind(x,cs),1,min)
d <- 1*(x <= cs)
tauxc <- 100*round( 1-mean(d),3)
cat(tauxc, "percents of data censored.\n")

c0 <- c(25, 180) # data-driven initial centers (visible modes)
sc0 <- 25/180    # and scaling
s <- spRMM_SEM(t, d, scaling = sc0, centers = c0, bw = 15, maxit = 100)

plotspRMM(s) # default
summary(s)   # S3 method for class "spRMM"
}
}


\keyword{file}
