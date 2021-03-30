\name{plotFDR}
\alias{plotFDR}

\title{Plot False Discovery Rate (FDR) estimates from output by EM-like strategies}

\description{Plot FDR\eqn{(p_i)} estimates against index of sorted p-values
	from, e.g., normalmixEM or the semiparametric mixture model posterior probabilities
	 output by \code{\link{spEMsymlocN01}}, or any EM-algorithm like
\code{\link{normalmixEM}} which returns posterior probabilities. The function
can simultaneously plot FDR estimates from two strategies for comparison.
Plot of the true FDR can be added if complete data are available 
(typically in simulation studies).}

\usage{
plotFDR(post1, post2 = NULL, lg1 = "FDR 1", lg2 = NULL, title = NULL, 
        compH0 = 1, alpha = 0.1, complete.data = NULL, pctfdr = 0.3)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{post1}{The matrix of posterior probabilities from objects such as the output
  	from \code{\link{spEMsymlocN01}}. The rows need to be sorted by increasing pvalues.}
  \item{post2}{A second object like \code{post1} if comparison is desired, also sorted by increasing pvalues.}
  \item{lg1}{Text describing the FDR estimate in \code{post1}.}
  \item{lg2}{Text describing the FDR estimate in \code{post2} if provided.}
  \item{title}{Plot title, a default is provided if \code{NULL}.}
  \item{compH0}{The component indicator associated to the null hypothesis H0,
  	normally 1 since it is defined in this way in \code{\link{spEMsymlocN01}}, but in case of label switching in other algorithms it can be set to \code{2}.}
  \item{alpha}{The target FDR level; the index at which the FDR estimate crosses the horizontal line for level \code{alpha} gives the maximum number of cases to reject.}
  \item{complete.data}{An array with \eqn{n} lines and 2 columns, with the component indicator in column 1 and the p-values in column 2, sorted by p-values.}
  \item{pctfdr}{The level up to which the FDR is plotted, i.e. 
  	the scale of the vertical axis.}
 }

\value{A plot of one or two FDR estimates, with the true FDR if available}

\references{
	   \itemize{
	\item Chauveau, D., Saby, N., Orton, T. G., Lemercier B., Walter, C. and Arrouys, D.
	Large-scale simultaneous hypothesis testing in monitoring carbon content 
  from French soil database -- A semi-parametric mixture approach, Geoderma 219-220 (2014), 117-124.
	}
}

\author{Didier Chauveau}

\seealso{
\code{\link{spEMsymlocN01}}}
%
\keyword{file}