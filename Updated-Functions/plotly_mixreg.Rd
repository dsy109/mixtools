\name{plotly_mixreg}
\alias{plotly_mixreg}

\title{Visualization of output of regressional \code{mixEM} function using \code{plotly}}

\description{This is an updated version of \code{plot.mixEM}. For more technical details, please refer to \code{plot.mixEM}.}

\usage{
plotly_mixreg <-function (x,
                          xlab="X", xlab.size=15 , xtick.size=15,
                          ylab="Y", ylab.size=15 , ytick.size=15,
                          title="Estimated Mixed Regression", title.size=15,
                          title.x = 0.5,title.y=0.95,
                          col="#1f77b4", lwd=3, cex=6)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{mixEM}.}
  \item{xlab}{Label of x-axis to be passed to the loglikelihood plot. Trying to change these parameters using \code{xlab} result in an error.} 
  \item{xlab.size}{Font of \code{xlab1}.}
  \item{xtick.size}{Font of tick labels of x-axis to be passed to the loglikelihood plot.}
  \item{ylab}{Label of y-axis to be passed to the loglikelihood plot. Trying to change these parameters using \code{ylab} result in an error.} 
  \item{ylab.size}{Font of \code{ylab1}.}
  \item{ytick.size}{Font of tick labels of y-axis to be passed to the loglikelihood plot.}
  \item{title}{Title to be passed to the loglikelihood plot.} 
  \item{title.size}{Tile size of the loglikelihood plot.}
  \item{title.x}{Horizontal position of the loglikelihood plot.}
  \item{title.y}{Verticle position of the loglikelihood plot.}
  \item{col}{Color of different components.}
  \item{lwd}{Width of the density curve of the loglikelihood plot.}
  \item{cex}{Dot size of the loglikelihood plot.}
  }

\value{A plot of the output of \code{mixEM} function is presented depends on output type.}

\author{Kedai Cheng}

\examples{
\dontrun{
data(NOdata)
attach(NOdata)
set.seed(100)
em.out <- regmixEM(Equivalence, NO, verb = TRUE, epsilon = 1e-04)
plotly_mixreg(em.out)  
}
}
%
\keyword{file}