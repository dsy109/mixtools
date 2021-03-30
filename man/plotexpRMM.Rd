\name{plotexpRMM}
\alias{plotexpRMM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plot sequences from the EM algorithm for censored mixture of exponentials
}
\description{Function for plotting sequences of estimates along iterations, from an object returned by the \code{\link{expRMM_EM}}, an EM algorithm for mixture of exponential 
distributions with randomly right censored data (see reference below).
}
\usage{
  plotexpRMM(a, title=NULL, rowstyle=TRUE, subtitle=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{a}{An object returned by \code{\link{expRMM_EM}}.}
  \item{title}{The title of the plot, set to some default value if \code{NULL}.}
  \item{rowstyle}{Window organization, for plots in rows (the default) or columns.}
  \item{subtitle}{A subtitle for the plot, set to some default value if \code{NULL}.}
  \item{...}{Other parameters (such as \code{lwd}) passed to \code{plot}, \code{lines}, and  
          \code{legend} commands.}
  }


\value{The plot returned}


\seealso{
  Related functions: 
  \code{\link{expRMM_EM}},   \code{\link{summary.mixEM}}, \code{\link{plot.mixEM}}.

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
   \url{http://link.springer.com/article/10.1007/s00180-016-0661-7}
   }
}
\author{Didier Chauveau}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
n=300 # sample size
m=2   # number of mixture components
lambda <- c(1/3,1-1/3); rate <- c(1,1/10) # mixture parameters
set.seed(1234)
x <- rexpmix(n, lambda, rate) # iid ~ exponential mixture
cs=runif(n,0,max(x)) # Censoring (uniform) and incomplete data
t <- apply(cbind(x,cs),1,min) # observed or censored data
d <- 1*(x <= cs)              # censoring indicator

###### EM for RMM, exponential lifetimes
l0 <- rep(1/m,m); r0 <- c(1, 0.5) # "arbitrary" initial values
a <- expRMM_EM(t, d, lambda=l0, rate=r0, k = m)
summary(a)             # EM estimates etc
plotexpRMM(a, lwd=2) # plot of EM sequences
%%\dontrun{
%%}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{file}
