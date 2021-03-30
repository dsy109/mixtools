\name{plot.npEM}
\title{Plot Nonparametric or Semiparametric EM Output}
\alias{plot.npEM} 
\alias{plot.spEM}
\usage{
\method{plot}{npEM}(x, blocks = NULL, hist=TRUE, addlegend = TRUE,
      scale=TRUE, title=NULL, breaks="Sturges", ylim=NULL, dens.col,
      newplot = TRUE, pos.legend = "topright", cex.legend = 1, \dots)         
\method{plot}{spEM}(x, blocks = NULL, hist=TRUE, addlegend = TRUE,
      scale=TRUE, title=NULL, breaks="Sturges", ylim=NULL, dens.col,
      newplot = TRUE, pos.legend = "topright", cex.legend = 1, \dots)
}
\description{
    Takes an object of class \code{npEM} and returns a set of plots of the 
    density estimates for each block and each component.  There is one plot
    per block, with all the components displayed on each block so it is possible
    to see the mixture structure for each block.      
} 
\arguments{
  \item{x}{An object of class \code{npEM} such as the output
  of the \code{\link{npEM}} function}
  \item{blocks}{Blocks (of repeated measures coordinates) to plot; not relevant
  for univariate case.  Default is to plot all blocks.}
  \item{hist}{If TRUE, superimpose density estimate plots on a histogram
    of the data}
  \item{addlegend}{If TRUE, adds legend to the plot.}
  \item{scale}{If TRUE, scale each density estimate by its corresponding estimated
  mixing proportion, so that the total area under all densities equals 1 and the 
  densities plotted may be added to produce an estimate of the mixture density.
  When FALSE, each density curve has area 1 in the plot.}
  \item{title}{Alternative vector of main titles for plots (recycled as many times
  as needed)}
  \item{breaks}{Passed directly to the \code{\link{hist}} function}
  \item{ylim}{\code{ylim} parameter to use for all plots, if desired.  If not
  given, each plot uses its own ylim that ensures that no part of the plot will
  go past the top of the plotting area.}
  \item{dens.col}{Color values to use for the individual component density
    functions, repeated as necessary.  Default value is \code{2:(m+1)}.}
  \item{newplot}{If TRUE, creates a new plot.}
  \item{pos.legend}{Single argument specifying the 
    position of the legend.  See `Details' section of
  \code{\link{legend}}.}
  \item{cex.legend}{Character expansion factor for \code{\link{legend}}.}
  \item{\dots}{Any remaining arguments are passed to the \code{\link{hist}} and 
  \code{\link{lines}} functions.}
}
\value{
  \code{plot.npEM} returns a list with two elements:
%  \item{means}{A \eqn{B\times m}{B x m} matrix of estimated population means,
%  where \eqn{B} is the number of blocks and \eqn{m} is the number of mixture components 
%  (subpopulations)}
%  \item{variances}{A \eqn{B\times m}{B x m} matrix of estimated population variances}
  \item{x}{List of matrices.  The \eqn{j}th column of the \eqn{i}th matrix is the vector of 
  \eqn{x}-values for the \eqn{j}th density in the \eqn{i}th plot.}
  \item{y}{\eqn{y}-values, given in the same form as the \eqn{x}-values.}
} 

\seealso{ 
\code{\link{npEM}}, \code{\link{density.npEM}}, \code{\link{spEMsymloc}},
\code{\link{plotseq.npEM}}
} 

\examples{ 
## Examine and plot water-level task data set.

## First, try a 3-component solution where no two coordinates are
## assumed i.d.
data(Waterdata)
set.seed(100)
\dontrun{
a <- npEM(Waterdata[,3:10], 3, bw=4)
par(mfrow=c(2,4))
plot(a) # This produces 8 plots, one for each coordinate
}

\dontrun{
## Next, same thing but pairing clock angles that are directly opposite one
## another (1:00 with 7:00, 2:00 with 8:00, etc.)
b <- npEM(Waterdata[,3:10], 3, blockid=c(4,3,2,1,3,4,1,2), bw=4)
par(mfrow=c(2,2))
plot(b) # Now only 4 plots, one for each block
}

}

\keyword{file}
