\name{plotly_npEM}
\title{Plot Nonparametric or Semiparametric EM Output}
\alias{plotly_npEM} 
\alias{plotly_spEM}
\usage{
plotly_npEM(x, blocks = NULL, hist=TRUE, addlegend=TRUE,
            scale = TRUE, title=NULL, breaks="Sturges", 
            dens.col = NULL, newplot=TRUE, ylim = NULL ,
            col.hist = "#1f77b4",
            width = 3, title.x = 0.5 , title.y = 0.95, title.size = 15,
            xlab = "X" , xlab.size = 15 , xtick.size = 15,
            ylab = "Density" , ylab.size = 15 , ytick.size = 15,
            legend.text = "Posteriors",
            legend.text.size = 15,
            legend.size = 15)         
plotly_spEM(x, blocks = NULL, hist=TRUE, addlegend=TRUE,
            scale = TRUE, title=NULL, breaks="Sturges", 
            dens.col = NULL, newplot=TRUE, ylim = NULL ,
            col.hist = "#1f77b4",
            width = 3, title.x = 0.5 , title.y = 0.95, title.size = 15,
            xlab = "X" , xlab.size = 15 , xtick.size = 15,
            ylab = "Density" , ylab.size = 15 , ytick.size = 15,
            legend.text = "Posteriors",
            legend.text.size = 15,
            legend.size = 15)
}
\description{
  This is an updater version of \code{plot.npEM} function by using \code{plotly}. For technical details, please refer to \code{plot.npEM}.     
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
  \item{col.hist}{Color of the histogram to plot.}
  \item{width}{Line width.}
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
}
\value{
  \code{plotly_npEM} returns a list with two elements:
  \item{x}{List of matrices.  The \eqn{j}th column of the \eqn{i}th matrix is the vector of 
  \eqn{x}-values for the \eqn{j}th density in the \eqn{i}th plot.}
  \item{y}{\eqn{y}-values, given in the same form as the \eqn{x}-values.}
} 

\seealso{ 
\code{\link{npEM}}, \code{\link{density.npEM}}, \code{\link{spEMsymloc}},
\code{\link{plotseq.npEM}}, \code{plot.npEM}
} 

\examples{ 
\dontrun{
## Examine and plot water-level task data set.

## First, try a 3-component solution where no two coordinates are
## assumed i.d.
data(Waterdata)
set.seed(100)
a <- npEM(Waterdata[,3:10], 3, bw=4)
plotly_npEM(a , newplot = FALSE)

## Next, same thing but pairing clock angles that are directly opposite one
## another (1:00 with 7:00, 2:00 with 8:00, etc.)
b <- npEM(Waterdata[,3:10], 3, blockid=c(4,3,2,1,3,4,1,2), bw=4)
plotly_npEM(b , newplot = FALSE)
}

}

\keyword{file}
