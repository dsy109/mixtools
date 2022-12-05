\name{plotly_ise.npEM}
\title{Visualization of Integrated Squared Error for a selected density from npEM output using \code{plotly}}
\alias{plotly_ise.npEM}
\usage{
plotly_ise.npEM(npEMout, component=1, block=1, truepdf=dnorm, lower=-Inf,
                upper=Inf, plots = TRUE ,
                col = NULL , width = 3,
                title = NULL , title.size = 15 , title.x = 0.5 , title.y = 0.95,
                xlab = "t" , xlab.size = 15 , xtick.size = 15,
                ylab = "" , ylab.size = 15 , ytick.size = 15,
                legend.text = "" , legend.text.size = 15 , legend.size = 15, ...)
}
\description{
This is an updated visualization function for \code{ise.npEM}. For more technical details, please refer to \code{ise.npEM}.
}
\arguments{
  \item{npEMout}{An object of class \code{npEM} such as the output
  of the \code{\link{npEM}} function}
  \item{component, block}{Component and block of particular density to analyze
  from \code{npEMout}.}
  \item{truepdf}{an \R function taking a numeric first argument and 
  returning a numeric vector of the same length. Returning a 
  non-finite element will generate an error.}
  \item{lower, upper}{the limits of integration.  Can be infinite.}
  \item{plots}{logical:  Should plots be produced?}
  \item{...}{additional arguments to be passed to \code{truepdf}
  	(and that may be mandatory like, e.g., the \code{df = } argument of \code{dt}). 
  Remember to use argument names not matching those of \code{ise.npRM}.}
  \item{col}{Color of traces.}
  \item{width}{Line width of traces.}
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
}
\value{
  Just as for the \code{\link{integrate}} function,
  a list of class \code{"integrate"} with components
  \item{value}{the final estimate of the integral.}
  \item{abs.error}{estimate of the modulus of the absolute error.}
  \item{subdivisions}{the number of subintervals produced in the
    subdivision process.}
  \item{message}{\code{"OK"} or a character string giving the error message.}
  \item{call}{the matched call.}
}
\details{This function calls the \code{\link{wkde}} (weighted kernel
density estimate) function.
}
\seealso{
\code{\link{npEM}}, \code{\link{wkde}}, \code{\link{integrate}}, \code{\link{ise.npEM}}
}
\references{
   \itemize{
   
  \item Benaglia, T., Chauveau, D., and Hunter, D. R. (2009), An EM-like algorithm
   for semi- and non-parametric estimation in multivariate mixtures, 
   Journal of Computational and Graphical Statistics, 18, 505-526.

   \item Benaglia, T., Chauveau, D., Hunter, D. R., and Young, D. (2009),
  mixtools: An R package for analyzing finite mixture models.
	Journal of Statistical Software, 32(6):1-29.


   }
}
\examples{
\dontrun{
data(Waterdata)
set.seed(100)
a <- npEM(Waterdata[,3:10], mu0=3, bw=4) # Assume indep but not iid
plotly_ise.npEM(a , plots = TRUE)
}
}


\keyword{file}
