\name{plotly_FDR}
\alias{plotly_FDR}

\title{Plot False Discovery Rate (FDR) estimates from output by EM-like strategies using \code{plotly}}

\description{This is an updated version of \code{plotFDR}. For more technical details, please refer to \code{plotFDR}.}

\usage{
plotly_FDR(post1, post2=NULL, lg1="FDR 1", lg2=NULL, 
          compH0=1, alpha=0.1, complete.data =NULL, pctfdr=0.3,
          col = NULL, width = 3 ,
          title = NULL , title.size = 15 , title.x = 0.5 , title.y = 0.95,
          xlab = "Index" , xlab.size = 15 , xtick.size = 15,
          ylab = "Probability" , ylab.size = 15 , ytick.size = 15,
          legend.text = "" , legend.text.size = 15 , legend.size = 15)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{post1}{The matrix of posterior probabilities from objects such as the output
  	from \code{\link{spEMsymlocN01}}. The rows need to be sorted by increasing pvalues.}
  \item{post2}{A second object like \code{post1} if comparison is desired, also sorted by increasing pvalues.}
  \item{lg1}{Text describing the FDR estimate in \code{post1}.}
  \item{lg2}{Text describing the FDR estimate in \code{post2} if provided.}
  \item{compH0}{The component indicator associated to the null hypothesis H0,
  	normally 1 since it is defined in this way in \code{\link{spEMsymlocN01}}, but in case of label switching in other algorithms it can be set to \code{2}.}
  \item{alpha}{The target FDR level; the index at which the FDR estimate crosses the horizontal line for level \code{alpha} gives the maximum number of cases to reject.}
  \item{complete.data}{An array with \eqn{n} lines and 2 columns, with the component indicator in column 1 and the p-values in column 2, sorted by p-values.}
  \item{pctfdr}{The level up to which the FDR is plotted, i.e. 
  	the scale of the vertical axis.}
  \item{col}{Color of traces.}
  \item{width}{Width of traces.}
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
\code{\link{spEMsymlocN01}}, \code{\link{plotFDR}}}

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
plotly_FDR(s1$post, s2$post, lg1 = "normalmixEM", lg2 = "spEMsymlocN01",
           complete.data = cpd) # with true FDR computed from z
}
%
\keyword{file}