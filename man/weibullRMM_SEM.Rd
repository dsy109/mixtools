\name{weibullRMM_SEM}
\title{St-EM algorithm for Reliability Mixture Models (RMM) of Weibull with right Censoring}
\alias{weibullRMM_SEM}

\usage{
weibullRMM_SEM(x, d = NULL, lambda = NULL, shape = NULL, scale = NULL,
        k = 2, maxit = 200, maxit.survreg = 200, epsilon = 1e-03,
        averaged = TRUE, verb = FALSE) }
\description{
  Parametric Stochastic EM (St-EM) algorithm for univariate finite mixture of Weibull 
  distributions with randomly right censored data.
}
\arguments{
  \item{x}{A vector of \eqn{n}{n} real positive lifetime (possibly censored) durations.
  If \code{d} is not \code{NULL} then a vector of random censoring times
  \code{c} occured,  
  so that \eqn{x= min(x,c)} and \eqn{d = I(x <= c)}.}
  
  \item{d}{The vector of censoring indication, where 1 means observed lifetime data, 
  and 0 means censored lifetime data.}
  
  \item{lambda}{Initial value of mixing proportions.
  If \code{NULL}, then \code{lambda} is set to \code{rep(1/k,k)}.}
  
  \item{shape}{Initial value of Weibull component shapes, 
    all set to 1 if  \code{NULL}.}
    
  \item{scale}{Initial value of Weibull component scales, 
    all set to 1 if  \code{NULL}.}

  \item{k}{Number of components of the mixture.}

  \item{maxit}{The number of iterations allowed, since for St-EM algorithms convergence
  is not based on stabilization, exactly \code{maxit} iterations are performed
  (see  Bordes L. and Chauveau D. (2016) reference below).}
  
    \item{maxit.survreg}{The number of iterations allowed in the computations of the
    MLE for censored weibull data from the  \code{survival} package
  (see  Bordes L. and Chauveau D. (2016) reference below).}
  
    \item{epsilon}{Tolerance parameter used in the numerical computations of the
    MLE for censored weibull data by \code{survreg} from the  \code{survival} package
  (see  Bordes L. and Chauveau D. (2016) reference below).}


  \item{averaged}{The way of updating parameters at each iteration: if \code{TRUE}, current values of the parameters 
  are obtained by averaging the sequence
  (see Bordes L. and Chauveau D. (2016) reference below).}

  \item{verb}{If TRUE, print updates for every iteration of the algorithm as
  it runs}
}

\details{This St-EM algorithm calls functions from  the \code{survival} package to compute 
        parametric MLE for censored weibull data.}
\value{
\code{weibullRMM_SEM} returns a list of class "mixEM" with the following items:
  \item{x}{The input data.}
  \item{d}{The input censoring indicator.}
  \item{lambda}{The estimates for the  mixing proportions.}
  \item{scale}{The estimates for the Weibull component scales.}
  \item{shape}{The estimates for the Weibull component shapes.}
  \item{loglik}{The log-likelihood value at convergence of the algorithm.}
  \item{posterior}{An \eqn{n\times k}{n x k} matrix of posterior probabilities for
   observation, after convergence of the algorithm.}
  \item{all.loglik}{The sequence of log-likelihoods over iterations.}
  \item{all.lambda}{The sequence of mixing proportions over iterations.}
  \item{all.scale}{The sequence of component scales over iterations.}
  \item{all.shape}{The sequence of component shapes over iterations.}
   \item{ft}{A character vector giving the name of the function called.}
}
\seealso{
Related functions:
\code{\link{plotweibullRMM}},   \code{\link{summary.mixEM}}.

Other models and algorithms for censored lifetime data 
(name convention is model_algorithm):
\code{\link{expRMM_EM}},
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
n = 500 # sample size
m = 2 # nb components
lambda=c(0.4, 0.6)
shape <- c(0.5,5); scale <- c(1,20) # model parameters
set.seed(321)
x <- rweibullmix(n, lambda, shape, scale) # iid ~ weibull mixture
cs=runif(n,0,max(x)+10) # iid censoring times
t <- apply(cbind(x,cs),1,min) # censored observations
d <- 1*(x <= cs)              # censoring indicator

## set arbitrary or "reasonable" (e.g., data-driven) initial values
l0 <- rep(1/m,m); sh0 <- c(1, 2); sc0 <- c(2,10)
# Stochastic EM algorithm 
a <- weibullRMM_SEM(t, d, lambda = l0, shape = sh0, scale = sc0, maxit = 200)

summary(a) # Parameters estimates etc
plotweibullRMM(a) # plot of St-EM sequences
plot(a, which=2) # or equivalently, S3 method for "mixEM" object
%%\dontrun{
%%}
}

\keyword{file}
