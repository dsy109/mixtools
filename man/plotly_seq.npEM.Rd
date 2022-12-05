\name{plotly_seq.npEM}
\alias{plotly_seq.npEM}
\title{Plotting sequences of estimates from non- or semiparametric EM-like Algorithm using \code{plotly}}

\description{This is an updated version of \code{\link{plotseq.npEM}}. For technical details, please refer to \code{\link{plotseq.npEM}}.
}

\usage{
  plotly_seq.npEM (x, col = '#1f77b4' , width = 6,
                   xlab = "Iteration" , xlab.size = 15 , xtick.size = 15,
                   ylab.size = 15 , ytick.size = 15,
                   title.size = 15 , title.x = 0.5 , title.y = 0.95)
}
\arguments{
  \item{x}{an object of class \code{npEM}, as output by \code{\link{npEM}} 
  or \code{\link{spEMsymloc}}}
  \item{col}{Line color.}
  \item{width}{Line width.}
  \item{title}{Text of the main title.}
  \item{title.size}{Size of the main title.}
  \item{title.x}{Horsizontal position of the main title.}
  \item{title.y}{Vertical posotion of the main title.}
  \item{xlab}{Label of X-axis.}
  \item{xlab.size}{Size of the lable of X-axis.}
  \item{xtick.size}{Size of tick lables of X-axis.}
  \item{ylab.size}{Size of the lable of Y-axis.}
  \item{ytick.size}{Size of tick lables of Y-axis.}
  }


\value{\code{plotly_seq.npEM} returns a figure with one plot for each component 
 proportion, and, in the case of \code{\link{spEMsymloc}}, one plot for each 
component mean.}

\seealso{
\code{\link{plot.npEM}}, \code{\link{rnormmix}},
\code{\link{npEM}}, \code{\link{spEMsymloc}}, \code{\link{plotly_seq.npEM}}
}


\references{
   \itemize{
   \item Benaglia, T., Chauveau, D., and Hunter, D. R. (2009), An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   Journal of Computational and Graphical Statistics (to appear).

   \item Bordes, L., Chauveau, D., and Vandekerkhove, P. (2007),
   An EM algorithm for a semiparametric mixture model, 
   Computational Statistics and Data Analysis, 51: 5429-5443.   
   }
}
\author{Didier Chauveau}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
\dontrun{
## Examine and plot water-level task data set.
## First, try a 3-component solution where no two coordinates are
## assumed i.d.
data(Waterdata)
set.seed(100)
## Not run:
a <- npEM(Waterdata[,3:10], mu0=3, bw=4) # Assume indep but not iid
plotly_seq.npEM(a)
}
}
\keyword{file}
