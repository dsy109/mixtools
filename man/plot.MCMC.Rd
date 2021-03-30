\name{plot.mixMCMC}
\title{Various Plots Pertaining to Mixture Model Output Using MCMC Methods}
\alias{plot.mixMCMC} 
\usage{ 
\method{plot}{mixMCMC}(x, trace.plots = TRUE, 
     summary.plots = FALSE, burnin = 2000, \dots) 
}

\description{
    Takes an object of class \code{mixMCMC} and returns various graphical output for select mixture models.
} 
\arguments{
  \item{x}{An object of class \code{mixMCMC}.}
  \item{trace.plots}{If TRUE, trace plots of the various parameters estimated by the MCMC methods is given.}
  \item{summary.plots}{Graphics pertaining to certain mixture models.  The details are given below.}
  \item{burnin}{The values 1 to \code{burnin} are dropped when producing the plots in \code{summary.plots}.}
  \item{...}{Graphical parameters passed to \code{regcr} function.}
}
\value{
  \code{plot.mixMCMC} returns trace plots of the various parameters estimated by the MCMC methods for all objects of class
  \code{mixMCMC}.  In addition, other plots may be produced for the following k-component mixture model functions:
  \item{regmixMH}{Credible bands for the regression lines in a mixture of linear regressions.  See \code{regcr} for more details.}
} 

\seealso{ 
\code{\link{regcr}} 
} 

\examples{ 
## M-H algorithm for NOdata with acceptance rate about 40\%.

data(NOdata)
attach(NOdata)
set.seed(100)
beta <- matrix(c(1.3, -0.1, 0.6, 0.1), 2, 2)
sigma <- c(.02, .05)
MH.out <- regmixMH(Equivalence, NO, beta = beta, s = sigma, 
                   sampsize = 2500, omega = .0013)
plot(MH.out, summary.plots = TRUE, burnin = 2450, 
     alpha = 0.01)
}

\keyword{file}
