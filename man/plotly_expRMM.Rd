\name{plotly_expRMM}
\alias{plotly_expRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot sequences from the EM algorithm for censored mixture of exponentials using \code{plotly}
}
\description{This is an updated function of \code{plotexpRMM}. For more technical details, please refer to \code{plotexpRMM}.
}
\usage{
  plotly_expRMM(a , title = NULL , rowstyle = TRUE , subtitle=NULL,
  width = 2 , cex = 2 , col.comp = NULL,
  legend.text = NULL, legend.text.size = 15, legend.size = 15,
  title.x = 0.5, title.y = 0.95, title.size = 15,
  xlab.size = 15, xtick.size = 15, 
  ylab.size = 15, ytick.size = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{An object returned by \code{\link{expRMM_EM}}.}
  \item{title}{The title of the plot, set to some default value if \code{NULL}.}
  \item{rowstyle}{Window organization, for plots in rows (the default) or columns.}
  \item{subtitle}{A subtitle for the plot, set to some default value if \code{NULL}.}
  \item{width}{Line width.}
  \item{cex}{Size of dots.}
  \item{col.comp}{Color of different components. Number of color specified needs to be consistent with number of components.}
  \item{legend.text}{Title of legend.}
  \item{legend.text.size}{Size of the legend title.}
  \item{legend.size}{Size of legend.}
  \item{title.size}{Size of the main title.}
  \item{title.x}{Horsizontal position of the main title.}
  \item{title.y}{Vertical posotion of the main title.}
  \item{xlab.size}{Size of the lable of X-axis.}
  \item{xtick.size}{Size of tick lables of X-axis.}
  \item{ylab.size}{Size of the lable of Y-axis.}
  \item{ytick.size}{Size of tick lables of Y-axis.}
  }


\value{The plot returned}


\seealso{
  Related functions: 
  \code{\link{expRMM_EM}},   \code{\link{summary.mixEM}}, \code{\link{plot.mixEM}}, \code{\link{plotexpRMM}}.

  Other models and algorithms for censored lifetime data 
(name convention is model_algorithm):
  \code{\link{weibullRMM_SEM}}, \code{\link{spRMM_SEM}}.
}

\references{
   \itemize{
   \item Bordes, L., and Chauveau, D. (2016),
   Stochastic EM algorithms for parametric and semiparametric mixture models 
   for right-censored lifetime data, 
   Computational Statistics, Volume 31, Issue 4, pages 1513-1538.
   \url{https://link.springer.com/article/10.1007/s00180-016-0661-7}
   }
}
\author{Didier Chauveau}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n=300 # sample size
m=2 # number of mixture components
lambda <- c(1/3,1-1/3); rate <- c(1,1/10) # mixture parameters
set.seed(1234)
x <- rexpmix(n, lambda, rate) # iid ~ exponential mixture
cs=runif(n,0,max(x)) # Censoring (uniform) and incomplete data
t <- apply(cbind(x,cs),1,min) # observed or censored data
d <- 1*(x <= cs) # censoring indicator
###### EM for RMM, exponential lifetimes
l0 <- rep(1/m,m); r0 <- c(1, 0.5) # "arbitrary" initial values
a <- expRMM_EM(t, d, lambda=l0, rate=r0, k = m)
summary(a) # EM estimates etc
plotly_expRMM(a , rowstyle = TRUE) # plot of EM sequences
%%\dontrun{
%%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
