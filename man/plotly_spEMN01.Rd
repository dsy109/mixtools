\name{plotly_spEMN01}
\alias{plotly_spEMN01}

\title{Plot mixture pdf for the semiparametric mixture model output by \code{spEMsymlocN01} using \code{plotly}.}

\description{This is an updated version of \code{plotlspEMN01} function by using \code{plotly}. For technical details, please refer to \code{\link{plot.spEMN01}}.}

\usage{
plotly_spEMN01(x, bw=x$bandwidth, knownpdf=dnorm, add.plot=FALSE,
               width = 3 , col.dens = NULL, col.hist =  '#1f77b4',
               title = NULL , title.size = 15 , 
               title.x = 0.5 , title.y = 0.95,
               xlab = "t" , xlab.size = 15 , xtick.size = 15,
               ylab = "Density" , ylab.size = 15 , ytick.size = 15,
               legend.text = "Densities" , legend.text.size = 15 , 
               legend.size = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "spEMN01" as returned by spEMsymlocN01}
  \item{bw}{Bandwidth for weighted kernel density estimation.}
  \item{knownpdf}{The known density of component 1, default to \code{dnorm}.}
  \item{add.plot}{Set to TRUE to add to an existing plot.}
  \item{width}{Line width.}
  \item{col.dens}{Color of density lines. Number of colors specified needs to be consistent with number of components.}
  \item{col.hist}{Color of histogram.}
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

\value{A plot of the density of the mixture}

\references{
	   \itemize{
	\item Chauveau, D., Saby, N., Orton, T. G., Lemercier B., Walter, C. and Arrouys, D.
	Large-scale simultaneous hypothesis testing in soil monitoring:
	A semi-parametric mixture approach, preprint (2013).
	}
}

\author{Didier Chauveau}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{spEMsymlocN01}}, \code{\link{plot.spEMN01}}}
%
\examples{
## Probit transform of p-values
## from a Beta-Uniform mixture model
## comparion of parametric and semiparametric EM fit
## Note: in actual situations n=thousands
set.seed(50)
n=300 # nb of multiple tests
m=2 # 2 mixture components
a=c(1,0.1); b=c(1,1); lambda=c(0.6,0.4) # parameters
z=sample(1:m, n, rep=TRUE, prob = lambda)
p <- rbeta(n, shape1 = a[z], shape2 = b[z]) # p-values
o <- order(p)
cpd <- cbind(z,p)[o,] # sorted complete data, z=1 if H0, 2 if H1
p <- cpd[,2] # sorted p-values
y <- qnorm(p) # probit transform of the pvalues
# gaussian EM fit with component 1 constrained to N(0,1)
s1 <- normalmixEM(y, mu=c(0,-4),
                  mean.constr = c(0,NA), sd.constr = c(1,NA))
s2 <- spEMsymlocN01(y, mu0 = c(0,-3)) # spEM with N(0,1) fit
plotly_spEMN01(s2 , add.plot = FALSE)
%%\dontrun{
%%}
}
\keyword{file}