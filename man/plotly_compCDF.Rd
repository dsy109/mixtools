\name{plotly_compCDF}
\title{Plot the Component CDF using \code{plotly}}
\alias{plotly_compCDF}
\usage{
plotly_compCDF(data, weights, x=seq(min(data, na.rm=TRUE), max(data, na.rm=TRUE), 
               len=250), comp=1:NCOL(weights), makeplot=TRUE,
               cex = 3, width = 3,
               legend.text = "Composition", legend.text.size = 15, legend.size = 15,
               title = "Empirical CDF", title.x = 0.5, title.y = 0.95, title.size = 15,
               xlab = "Data", xlab.size = 15, xtick.size = 15,
               ylab = "Probability", ylab.size = 15, ytick.size = 15,
               col.comp = NULL) 
}
\description{
Plot the components' CDF via the posterior probabilities using \code{plotly}.
}
\arguments{
  \item{data}{A matrix containing the raw data. Rows are subjects and columns
   are repeated measurements.}
  \item{weights}{The weights to compute the empirical CDF; however, most of
   time they are the posterior probabilities.}
  \item{x}{The points at which the CDFs are to be evaluated.}
  \item{comp}{The mixture components for which CDFs are desired.}
  \item{makeplot}{Logical:  Should a plot be produced as a side effect?}
  \item{cex}{Size of markers.}
  \item{width}{Line width.}
  \item{title}{Text of the main title.}
  \item{title.size}{Size of the main title.}
  \item{title.x}{Horsizontal position of the main title.}
  \item{title.y}{Vertical posotion of the main title.}
  \item{xlab}{Label of X-axis.}
  \item{xlab.size}{Size of the lable of X-axis.}
  \item{xtick.size}{Size of tick lables of X-axis.}
  \item{ylab}{Label of Y-axis.}
  \item{ylab.size}{Size of the lable of Y-axis.}
  \item{ytick.size}{Size of tick lables of Y-axis.}
  \item{legend.text}{Title of legend.}
  \item{legend.text.size}{Size of the legend title.}
  \item{legend.size}{Size of legend.}
  \item{col.comp}{Color of compositions. Number of color specified needs to be consistent with number of compositions.}
}
\value{
  A matrix with \code{length(comp)} rows and \code{length(x)} columns 
  in which each row gives the CDF evaluated at each point of \code{x}.
}
\details{
  When \code{makeplot} is \code{TRUE}, a line plot is produced of the
  CDFs evaluated at \code{x}.  The plot is not a step function plot; 
  the points \eqn{(x, CDF(x))} are simply joined by line segments.
}
\references{
  McLachlan, G. J. and Peel, D. (2000) \emph{Finite Mixture Models}, John Wiley and Sons, Inc.
  
  Elmore, R. T., Hettmansperger, T. P. and Xuan, F. (2004) The Sign Statistic, One-Way Layouts
  and Mixture Models, \emph{Statistical Science} \bold{19(4)}, 579--587.
}
\seealso{
\code{\link{makemultdata}}, \code{\link{multmixmodel.sel}}, \code{\link{multmixEM}}, \code{\link{compCDF}}.
}
\examples{
## The sulfur content of the coal seams in Texas
set.seed(100)
A <- c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B <- c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59)
C <- c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49)
D <- c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32)
E <- c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)
dis.coal <- makemultdata(A, B, C, D, E,
                         cuts = median(c(A, B, C, D, E)))
temp <- multmixEM(dis.coal)
## Now plot the components' CDF via the posterior probabilities
plotly_compCDF(dis.coal$x, temp$posterior, xlab="Sulfur")
}

\keyword{file}

