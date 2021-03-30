\name{plot.mixEM}
\title{Various Plots Pertaining to Mixture Models}
\alias{plot.mixEM} 
\usage{
\method{plot}{mixEM}(x, whichplots = 1, 
     loglik = 1 \%in\% whichplots,
     density = 2 \%in\% whichplots,
     xlab1="Iteration", ylab1="Log-Likelihood",
     main1="Observed Data Log-Likelihood", col1=1, lwd1=2,
     xlab2=NULL, ylab2=NULL, main2=NULL, col2=NULL, 
     lwd2=2, alpha = 0.05, marginal = FALSE, ...)
}

\description{
    Takes an object of class \code{mixEM} and returns various graphical output for select mixture models.
} 
\arguments{
  \item{x}{An object of class \code{mixEM}.}
  \item{whichplots}{vector telling which plots to produce:  1 = loglikelihood
        plot, 2 = density plot.  Irrelevant if \code{loglik} and \code{density}
        are specified.}
  \item{loglik}{If TRUE, a plot of the log-likelihood versus the EM iterations is given.}
  \item{density}{Graphics pertaining to certain mixture models.  The details are given below.}
  \item{xlab1, ylab1, main1, col1, lwd1}{Graphical parameters \code{xlab}, ..., \code{lwd}
  to be passed to the loglikelihood plot.  Trying to change these parameters using
  \code{xlab}, ..., \code{lwd} will result in an error, but all other graphical parameters
  are passed directly to the plotting functions via ...}
  \item{xlab2, ylab2, main2, col2, lwd2}{Same as \code{xlab1} etc. but for the
  density plot}
  \item{alpha}{A vector of significance levels when constructing confidence ellipses and confidence bands for the mixture
  of multivariate normals and mixture of regressions cases, respectively.  The default is 0.05.}
  \item{marginal}{For the mixture of bivariate normals, should optional marginal histograms be included?}
  \item{...}{Graphical parameters passed to \code{plot} command.}
}
\value{
  \code{plot.mixEM} returns a plot of the log-likelihood versus the EM iterations by default for all objects of class
  \code{mixEM}.  In addition, other plots may be produced for the following k-component mixture model functions:
  \item{normalmixEM}{A histogram of the raw data is produced along with k density curves determined by \code{normalmixEM}.}
  \item{repnormmixEM}{A histogram of the raw data produced in a similar manner as for \code{normalmixEM}.}
  \item{mvnormalmixEM}{A 2-dimensional plot with each point color-coded to denote its most probable component membership. In
  addition, the estimated component means are plotted along with (1 - \code{alpha})\% bivariate normal density contours.  These ellipses are
  constructed by assigning each value to their component of most probable membership and then using normal theory. Optional marginal histograms
  may also be produced.}
  \item{regmixEM}{A plot of the response versus the predictor with each point color-coded to denote its most probable component
  membership.  In addition, the estimated component regression lines are plotted along with (1 - \code{alpha})\% Working-Hotelling 
  confidence bands. These bands are constructed by assigning each value to their component of most probable membership and then
  performing least squares estimation.}
  \item{logisregmixEM}{A plot of the binary response versus the predictor with each point color-coded to denote its most probable
  compopnent membership.  In addition, the estimate component logistic regression lines are plotted.}
  \item{regmixEM.mixed}{Provides a 2x2 matrix of plots summarizing the posterior slope and posterior intercept terms from a
  mixture of random effects regression.  See \code{post.beta} for a more detailed description.}
} 

\seealso{ 
\code{\link{post.beta}} 
} 

\examples{ 
##Analyzing the Old Faithful geyser data with a 2-component mixture of normals.

data(faithful)
attach(faithful)
set.seed(100)
out <- normalmixEM(waiting, arbvar = FALSE, verb = TRUE,
                   epsilon = 1e-04)
plot(out, density = TRUE, w = 1.1)

##Fitting randomly generated data with a 2-component location mixture of bivariate normals.

x.1 <- rmvnorm(40, c(0, 0))
x.2 <- rmvnorm(60, c(3, 4))
X.1 <- rbind(x.1, x.2)

out.1 <- mvnormalmixEM(X.1, arbvar = FALSE, verb = TRUE,
                       epsilon = 1e-03)
plot(out.1, density = TRUE, alpha = c(0.01, 0.05, 0.10), 
     marginal = TRUE)

}

\keyword{file}
