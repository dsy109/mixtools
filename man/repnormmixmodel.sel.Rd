\name{repnormmixmodel.sel}
\title{Model Selection in Mixtures of Normals with Repeated Measures}
\alias{repnormmixmodel.sel}
\usage{
repnormmixmodel.sel(x, k = 2, ...)
}
\description{
  Assess the number of components in a mixture model with normal components and repeated measures 
  using the Akaike's information criterion (AIC), Schwartz's Bayesian information criterion (BIC), 
  Bozdogan's consistent AIC (CAIC), and Integrated Completed Likelihood (ICL).
}
\arguments{
  \item{x}{An mxn matrix of observations. The rows correspond to the repeated measures and the columns correspond
  to the subject.}
  \item{k}{The maximum number of components to assess.}
  \item{...}{Additional arguments passed to \code{repnormmixEM}.}
}
\value{
  \code{repnormmixmodel.sel} returns a matrix of the AIC, BIC, CAIC, and ICL values along with the winner (i.e., the highest
  value given by the model selection criterion) for a mixture of normals with repeated measures.
}
\seealso{
  \code{\link{repnormmixEM}}
}
\references{
  Biernacki, C., Celeux, G., and Govaert, G. (2000). Assessing a Mixture Model for Clustering with the
  Integrated Completed Likelihood. \emph{IEEE Transactions on Pattern Analysis and Machine Intelligence}, 22(7):719-725.
  
  Bozdogan, H. (1987). Model Selection and Akaike's Information Criterion (AIC): The General Theory and its
  Analytical Extensions. \emph{Psychometrika}, 52:345-370.
}
\examples{
## Assessing the number of components for the water-level task data set.

data(Waterdata)
water<-t(as.matrix(Waterdata[,3:10]))
set.seed(100)
out <- repnormmixmodel.sel(water, k = 3, epsilon = 5e-01)
out
}

\keyword{file}
