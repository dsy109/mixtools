\name{plotly_mixMCMC}
\title{Various Plots Pertaining to Mixture Model Output Using MCMC Methods using \code{plotly}}
\alias{plotly_mixMCMC} 
\usage{ 
plotly_mixMCMC(x, trace.plot = TRUE, summary.plot = FALSE, burnin = 2000, 
               credit.region = 0.95, col.cr = NULL,
               cex.trace = 3, width.trace = 3, 
               cex.summary = 3, width.summary = 1,
               title.trace = "", title.trace.x = 0.5, 
               title.trace.y = 0.95, title.trace.size = 15,
               xlab.trace = "Index", xlab.trace.size = 15, xtick.trace.size = 15,
               ylab.trace = NULL, ylab.trace.size = 15, ytick.trace.size = 15,
               title.summary = "Credible Regions", title.summary.x = 0.5, 
               title.summary.y = 0.95, title.summary.size = 15,
               xlab.summary = "Predictor", xlab.summary.size = 15, 
               xtick.summary.size = 15,
               ylab.summary = "Response", ylab.summary.size = 15, 
               ytick.summary.size = 15
) 
}

\description{
    This is an updated version of \code{plot.mixMCMC}. For technical details, please refer to \code{plot.mixMCMC}.
} 
\arguments{
  \item{x}{An object of class \code{mixMCMC}.}
  \item{trace.plot}{If TRUE, trace plots of the various parameters estimated by the MCMC methods is given.}
  \item{summary.plot}{Graphics pertaining to certain mixture models.  The details are given below.}
  \item{burnin}{The values 1 to \code{burnin} are dropped when producing the plots in \code{summary.plots}.}
  \item{credit.region}{Confidence level of credit region.}
  \item{col.cr}{Color of credit region. Number of color specified needs to be consistent with number of components.}
  \item{cex.trace}{Dot size of trace plots.}
  \item{width.trace}{Line width of trace plots.}
  \item{cex.summary}{Dot size of summary plots.}
  \item{width.summary}{Line width of summary plots.}
  \item{title.trace}{Text of the main title of trace plots.}
  \item{title.trace.x}{Horizontal position of main title of trace plots.}
  \item{title.trace.y}{Vertical position of main title of trace plots.}
  \item{title.trace.size}{Text sise of main title of trace plots.}
  \item{xlab.trace}{Label of X-axis of trace plots.}
  \item{xlab.trace.size}{Size of the lable of X-axis of trace plots.}
  \item{xtick.trace.size}{Size of tick lables of X-axis of trace plots.}
  \item{ylab.trace}{Label of Y-axis of trace plots.}
  \item{ylab.trace.size}{Size of the lable of Y-axis of trace plots.}
  \item{ytick.trace.size}{Size of tick lables of Y-axis of trace plots.}
  \item{title.summary}{Text of the main title of summar plot.}
  \item{title.summary.x}{Horizontal position of main title of summary plot.}
  \item{title.summary.y}{Vertical position of main title of summary plot.}
  \item{title.summary.size}{Text sise of main title of summary plot.}
  \item{xlab.summary}{Label of X-axis of summary plot.}
  \item{xlab.summary.size}{Size of the lable of X-axis of summary plot.}
  \item{xtick.summary.size}{Size of tick lables of X-axis of summary plot.}
  \item{ylab.summary}{Label of Y-axis of summary plot.}
  \item{ylab.summary.size}{Size of the lable of Y-axis of summary plot.}
  \item{ytick.summary.size}{Size of tick lables of Y-axis of summary plot.}

}
\value{
  \code{plotly_mixMCMC} returns trace plots of the various parameters estimated by the MCMC methods for all objects of class
  \code{mixMCMC}.  In addition, other plots may be produced for the following k-component mixture model functions:
  \item{regmixMH}{Credible bands for the regression lines in a mixture of linear regressions.  See \code{regcr} for more details.}
} 

\seealso{ 
\code{\link{regcr}}, \code{plot.mixMCMC} 
} 

\examples{ 
\dontrun{
data(NOdata)
attach(NOdata)
set.seed(100)
beta <- matrix(c(1.3, -0.1, 0.6, 0.1), 2, 2)
sigma <- c(.02, .05)
MH.out <- regmixMH(Equivalence, NO, beta = beta, s = sigma,
                   sampsize = 2500, omega = .0013)
plotly_mixMCMC(x = MH.out, summary.plot = TRUE, col.cr = c("red", "green"))
}
}

\keyword{file}
