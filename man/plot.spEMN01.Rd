\name{plot.spEMN01}
\alias{plot.spEMN01}

\title{Plot mixture pdf for the semiparametric mixture model output by spEMsymlocN01}

\description{Plot mixture density for the semiparametric mixture model output by spEMsymlocN01, with one component known and set to normal(0,1), and a symmetric nonparametric density with location parameter.}

\usage{
\method{plot}{spEMN01}(x, bw = x$bandwidth, knownpdf = dnorm, add.plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class "spEMN01" as returned by spEMsymlocN01}
  \item{bw}{Bandwidth for weighted kernel density estimation.}
  \item{knownpdf}{The known density of component 1, default to \code{dnorm}.}
  \item{add.plot}{Set to TRUE to add to an existing plot.}
  \item{...}{further arguments passed to \code{plot} if \code{add.plot = FALSE},
  	and to \code{lines} if \code{add.plot = TRUE}.}
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
\code{\link{spEMsymlocN01}}}
%
\keyword{file}