\name{plotweibullRMM}
\alias{plotweibullRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot sequences from the Stochastic EM algorithm for mixture of Weibull
}

\description{Function for plotting sequences of estimates along iterations, from an object returned by \code{\link{weibullRMM_SEM}}, a Stochastic EM algorithm for mixture of Weibull 
distributions with randomly right censored data (see reference below).
}

\usage{
  plotweibullRMM(a, title = NULL, rowstyle = TRUE, subtitle = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{An object returned by \code{\link{weibullRMM_SEM}}.}
  \item{title}{The title of the plot, set to some default value if \code{NULL}.}
  \item{rowstyle}{Window organization, for plots in rows (the default) or columns.}
  \item{subtitle}{A subtitle for the plot, set to some default value if \code{NULL}.}
  \item{...}{Other parameters (such as \code{lwd}) passed to \code{plot}, \code{lines}, and  
          \code{legend} commands.}
  }


\value{The plot returned}

\seealso{
Related functions:
\code{\link{weibullRMM_SEM}},   \code{\link{summary.mixEM}}.

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
   \url{http://link.springer.com/article/10.1007/s00180-016-0661-7}
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
d <- 1*(x <= cs)              # censoring indicator

## set arbitrary or "reasonable" (e.g., data-driven) initial values
l0 <- rep(1/m,m); sh0 <- c(1, 2); sc0 <- c(2,10)
# Stochastic EM algorithm 
a <- weibullRMM_SEM(t, d, lambda = l0, shape = sh0, scale = sc0, maxit = 200)

summary(a) # Parameters estimates etc
plotweibullRMM(a) # default plot of St-EM sequences
%%\dontrun{
%%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
