\name{mixtools-internal}
\alias{inv.logit}
\alias{dexpmixt}
\alias{HRkde}
\alias{kern.B}
\alias{kern.C}
\alias{kern.G}
\alias{kern.O}
\alias{kern.T}
\alias{kfoldCV}
\alias{KMintegrate}
\alias{KMod}
\alias{ldc}
\alias{logit}
\alias{npMSL_old}
\alias{plotseq}
\alias{rlnormscalemix}
\alias{splitsample}
\alias{triang_wkde}
\alias{wbw.kCV}

\title{Internal 'mixtools' Functions}
\description{
  Internal kernel, semiparametric-related, and miscellaneous functions for the package \code{mixtools}.
}
\usage{
dexpmixt(t, lam, rate)
HRkde(cpd, u = cpd[,1], kernelft = triang_wkde, 
      bw = rep(bw.nrd0(as.vector(cpd[,1])), length(cpd[,1])))
inv.logit(eta)
kern.B(x, xi, h, g = 0)
kern.C(x, xi, h)
kern.G(x, xi, h)
kern.O(x, xi, h)
kern.T(x, xi, h)
kfoldCV(h, x, nbsets = 2, w = rep(1, length(x)), 
        lower = mean(x) - 5*sd(x), upper = mean(x) + 5*sd(x))
KMintegrate(s) 
KMod(cpd, already.ordered = TRUE)
ldc(data, class, score)
logit(mu)
npMSL_old(x, mu0, blockid = 1:ncol(x),
          bw=bw.nrd0(as.vector(as.matrix(x))), samebw = TRUE,
          h=bw, eps=1e-8, maxiter=500, bwiter = maxiter,
          ngrid = 200, post = NULL, verb = TRUE)
plotseq(x, ...)
rlnormscalemix(n, lambda=1, meanlog=1, sdlog=1, scale=0.1)
splitsample(n, nbsets = 2)
triang_wkde(t, u=t, w=rep(1/length(t),length(t)), bw=rep(bw.nrd0(t), length(t)))
wbw.kCV(x, nbfold = 5, w = rep(1, length(x)), 
        hmin = 0.1*hmax, hmax = NULL)
}

\arguments{
  \item{x}{A vector of values to which local modeling techniques are applied.}
  \item{xi}{An n-vector of data values.}
  \item{h}{The bandwidth controlling the size of the window used for the
  local estimation around \code{x}.  This pertains to its usage in the kernel functionns
  \code{kern.B}, \code{kern.C}, \code{kern.G}, \code{kern.O}, and \code{kern.T}.  For its
  usage in the \code{kfoldCV} function, see updated arguments in the \code{npMSL} function.}
  \item{g}{A shape parameter required for the symmetric beta kernel.  The default
  is \code{g} = 0 which yields the uniform kernel.  Some common values are \code{g} = 1 for the
  Epanechnikov kernel, \code{g} = 2 for the biweight kernel, and \code{g} = 3 for the triweight kernel.}
  \item{mu0}{See updated arguments in the \code{npMSL} function.}
  \item{blockid}{See updated arguments in the \code{npMSL} function.}
  \item{bw}{See updated arguments in the \code{npMSL} function.}
  \item{samebw}{See updated arguments in the \code{npMSL} function.}
  \item{eps}{See updated arguments in the \code{npMSL} function.}
  \item{maxiter}{See updated arguments in the \code{npMSL} function.}
  \item{bwiter}{See updated arguments in the \code{npMSL} function.}
  \item{ngrid}{See updated arguments in the \code{npMSL} function.}
  \item{post}{See updated arguments in the \code{npMSL} function.}
  \item{verb}{See updated arguments in the \code{npMSL} function.}
  \item{n}{See updated arguments in the \code{npMSL} function.}
  \item{nbsets}{See updated arguments in the \code{npMSL} function.}
  \item{w}{See updated arguments in the \code{npMSL} function.}
  \item{lower}{See updated arguments in the \code{npMSL} function.}
  \item{upper}{See updated arguments in the \code{npMSL} function.}
  \item{nbfold}{See updated arguments in the \code{npMSL} function.}
  \item{hmin}{See updated arguments in the \code{npMSL} function.}
  \item{hmax}{See updated arguments in the \code{npMSL} function.}
  \item{data}{Data, possibly multivariate, fed to the \code{mixturegram} function.}
  \item{class}{The number of classes, inputted based on number of components in the \code{mixturegram} function.}
  \item{score}{The score vector from LDA used in constructing a mixturegram.}
  \item{lam}{A vector of mixture proportions, should sum to one.}
  \item{rate}{A vector of mixture component rates.}
  \item{t}{Argument for \code{dexpmixt}.}
  \item{mu}{A proportion for which to calculate the logit function; i.e., \code{log(mu / (1 - mu))}.}
  \item{eta}{Any real value for which to calculate the inverse logit function; 
        i.e., \code{1 / (1 + exp(eta))}.}  
  \item{cpd}{Argument for \code{HRkde}.}
  \item{kernelft}{Argument for \code{HRkde}.} 
  \item{s}{Argument for \code{KMintegrate}.}
  \item{meanlog}{Argument for \code{rlnormscalemix}.}
  \item{sdlog}{Argument for \code{rlnormscalemix}.}
}

\details{
  These are usually not to be called by the user.
}

\seealso{
  \code{\link{npMSL}}
}


\keyword{internal}
