\name{plotly_spRMM}
\alias{plotly_spRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot output from Stochastic EM algorithm for semiparametric scaled mixture of censored data using \code{plotly}.
}

\description{This is an updated version of \code{plotspRMM} function. For technical details, please refer to \code{\link{plotspRMM}.
}}

\usage{
  plotly_spRMM(sem, tmax = NULL,
               width = 3 , col = '#1f77b4', cex = 3,
               title.size = 15 , 
               title.x = 0.5 , title.y = 0.95,
               xlab.size = 15 , xtick.size=15 ,
               ylab.size = 15 , ytick.size=15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sem}{An object returned by \code{\link{spRMM_SEM}}.}
  \item{tmax}{The max time for \eqn{x} axis, set to some default value if \code{NULL}.}
  \item{width}{Width of lines.}
  \item{col}{Color of lines.}
  \item{cex}{Size of dots.}
  \item{title.size}{Size of the main title.}
  \item{title.x}{Horizontal position of the main title.}
  \item{title.y}{Vertical position of the main title.}
  \item{xlab.size}{Size of the label of X-axis.}
  \item{xtick.size}{Size of the tick of X-axis.}
  \item{ylab.size}{Size of the label of Y-axis.}
  \item{ytick.size}{Size of the tick of Y-axis.}
  }
  


\value{The four plots returned.}


\seealso{
Related functions: \code{\link{spRMM_SEM}} , \code{\link{plotspRMM}}.

Other models and algorithms for censored lifetime data 
(name convention is model_algorithm):
\code{\link{expRMM_EM}},
\code{\link{weibullRMM_SEM}}.
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
\dontrun{
n=500 # sample size
m=2 # nb components
lambda=c(0.4, 0.6) # parameters
meanlog=3; sdlog=0.5; scale=0.1
set.seed(12)
# simulate a scaled mixture of lognormals
x <- rlnormscalemix(n, lambda, meanlog, sdlog, scale)
cs=runif(n,20,max(x)+400) # Censoring (uniform) and incomplete data
t <- apply(cbind(x,cs),1,min)
d <- 1*(x <= cs)
tauxc <- 100*round( 1-mean(d),3)
cat(tauxc, "percents of data censored.\n")

c0 <- c(25, 180) # data-driven initial centers (visible modes)
sc0 <- 25/180    # and scaling
s <- spRMM_SEM(t, d, scaling = sc0, centers = c0, bw = 15, maxit = 100)

plotly_spRMM(s) # default
summary(s)   # S3 method for class "spRMM"
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
