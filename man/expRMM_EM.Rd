\name{expRMM_EM}
\title{EM algorithm for Reliability Mixture Models (RMM) with right Censoring}
\alias{expRMM_EM}
\usage{
  expRMM_EM(x, d=NULL, lambda = NULL, rate = NULL, k = 2, 
		    complete = "tdz", epsilon = 1e-08, maxit = 1000, verb = FALSE) }

\description{
  Parametric EM algorithm for univariate finite mixture of exponentials 
  distributions with randomly right censored data.
}
\arguments{
  \item{x}{A vector of \eqn{n}{n} real positive lifetime (possibly censored) durations.
  If \code{d} is not \code{NULL} then a vector of random censoring times
  \code{c} occurred,  
  so that \eqn{x= min(x,c)} and \eqn{d = I(x <= c)}.}
  
  \item{d}{The vector of censoring indication, where 1 means observed lifetime data, 
  and 0 means censored lifetime data.}
  
  \item{lambda}{Initial value of mixing proportions.
  If \code{NULL}, then \code{lambda} is set to \code{rep(1/k,k)}.}
  
  \item{rate}{Initial value of component exponential rates, 
    all set to 1 if  \code{NULL}.}

  \item{k}{Number of components of the mixture.}

  \item{complete}{Nature of complete data involved within the EM machinery,
  can be "tdz" for \code{(t,d,z)} (the default), or "xz" for \code{(x,z)}
  (see  Bordes L. and Chauveau D. (2016) reference below).}
  
  \item{epsilon}{Tolerance limit for declaring algorithm convergence based on
  the change between two consecutive iterations.}

  \item{maxit}{The maximum number of iterations allowed, convergence
  may be declared before \code{maxit} iterations (see \code{epsilon} above).}
  
  \item{verb}{If TRUE, print updates for every iteration of the algorithm as
  it runs}
}

\value{
\code{expRMM_EM} returns a list of class "mixEM" with the following items:
  \item{x}{The input data.}
  \item{d}{The input censoring indicator.}
  \item{lambda}{The estimates for the mixing proportions.}
  \item{rate}{The estimates for the component rates.}
  \item{loglik}{The log-likelihood value at convergence of the algorithm.}
  \item{posterior}{An \eqn{n\times k}{n x k} matrix of posterior probabilities for
   observation, after convergence of the algorithm.}
  
  \item{all.loglik}{The sequence of log-likelihoods over iterations.}
  \item{all.lambda}{The sequence of mixing proportions over iterations.}
  
  \item{all.rate}{The sequence of component rates over iterations.}

   \item{ft}{A character vector giving the name of the function.}
}
\seealso{
Related functions: 
  \code{\link{plotexpRMM}},
  \code{\link{summary.mixEM}}.
  
  Other models and algorithms for censored lifetime data:
  \code{\link{weibullRMM_SEM}}, 
  \code{\link{spRMM_SEM}}.
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
n <- 300 # sample size
m <- 2   # number of mixture components
lambda <- c(1/3,1-1/3); rate <- c(1,1/10) # mixture parameters
set.seed(1234)
x <- rexpmix(n, lambda, rate) # iid ~ exponential mixture
cs <- runif(n,0,max(x)) # Censoring (uniform) and incomplete data
t <- apply(cbind(x,cs),1,min) # observed or censored data
d <- 1*(x <= cs)              # censoring indicator

###### EM for RMM, exponential lifetimes
l0 <- rep(1/m,m); r0 <- c(1, 0.5) # "arbitrary" initial values
a <- expRMM_EM(t, d, lambda = l0, rate = r0, k = m)
summary(a)                 # EM estimates etc
plotexpRMM(a, lwd=2) # default plot of EM sequences
plot(a, which=2) # or equivalently, S3 method for "mixEM" object
%%\dontrun{
%%}
}


\keyword{file}
