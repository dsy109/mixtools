\name{plotly_weibullRMM}
\alias{plotly_weibullRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot sequences from the Stochastic EM algorithm for mixture of Weibull using \code{plotly}
}

\description{This is an updated version of \code{plotweibullRMM} function by using \code{plotly} function. For technical details, please refer to \code{\link{plotweibullRMM}}.
}

\usage{
  plotly_weibullRMM(a, title=NULL, rowstyle=TRUE, subtitle=NULL,
                    width = 3 , col = NULL , 
                    title.size = 15 , title.x = 0.5 , title.y = 0.95,
                    xlab = "Iterations" , xlab.size = 15 , xtick.size = 15,
                    ylab = "Estimates" , ylab.size = 15 , ytick.size = 15,
                    legend.size = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{An object returned by \code{\link{weibullRMM_SEM}}.}
  \item{title}{The title of the plot, set to some default value if \code{NULL}.}
  \item{rowstyle}{Window organization, for plots in rows (the default) or columns.}
  \item{subtitle}{A subtitle for the plot, set to some default value if \code{NULL}.}
  \item{width}{Line width.}
  \item{col}{Color of lines. Number of colors specified needs to be consistent with number of components.}
  \item{title.size}{Size of the main title.}
  \item{title.x}{Horsizontal position of the main title.}
  \item{title.y}{Vertical posotion of the main title.}
  \item{xlab}{Label of X-axis.}
  \item{xlab.size}{Size of the lable of X-axis.}
  \item{xtick.size}{Size of tick lables of X-axis.}
  \item{ylab}{Label of Y-axis.}
  \item{ylab.size}{Size of the lable of Y-axis.}
  \item{ytick.size}{Size of tick lables of Y-axis.}
  \item{legend.size}{Size of legend.}
  }


\value{The plot returned.}

\seealso{
Related functions:
\code{\link{weibullRMM_SEM}},   \code{\link{summary.mixEM}}, \code{\link{plotweibullRMM}}.

Other models and algorithms for censored lifetime data 
(name convention is model_algorithm):
\code{\link{expRMM_EM}},
\code{\link{spRMM_SEM}}
.
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
n = 500 # sample size
m = 2 # nb components
lambda=c(0.4, 0.6)
shape <- c(0.5,5); scale <- c(1,20) # model parameters
set.seed(321)
x <- rweibullmix(n, lambda, shape, scale) # iid ~ weibull mixture
cs=runif(n,0,max(x)+10) # iid censoring times
t <- apply(cbind(x,cs),1,min) # censored observations
d <- 1*(x <= cs) # censoring indicator
## set arbitrary or "reasonable" (e.g., data-driven) initial values
l0 <- rep(1/m,m); sh0 <- c(1, 2); sc0 <- c(2,10)
# Stochastic EM algorithm
a <- weibullRMM_SEM(t, d, lambda = l0, shape = sh0, scale = sc0, maxit = 200)
summary(a) # Parameters estimates etc
plotly_weibullRMM(a , legend.size = 20) # plot of St-EM sequences
%%\dontrun{
%%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
