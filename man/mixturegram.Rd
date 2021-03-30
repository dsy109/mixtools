\name{mixturegram}
\title{Mixturegrams}
\alias{mixturegram}
\usage{
mixturegram(data, pmbs, method = c("pca", "kpca", "lda"), all.n = FALSE,
			id.con = NULL, score = 1, iter.max = 50, nstart = 25, ...)
}

\description{
  Construct a mixturegram for determining an apporpriate number of components.
}
\arguments{
  \item{data}{The data, which must either be a vector or a matrix.  If a matrix, then the rows correspond to the observations.}
  \item{pmbs}{A list of length (K-1) such that each element is an nxk matrix of the posterior membership probabilities.  These are obtained from each of the "best" estimated k-component mixture models, k = 2,...,K.
}
  \item{method}{The dimension reduction method used.  \code{method = "pca"} implements principal components analysis.  \code{method = "kpca"} implements kernel principal components analysis. \code{method = "lda"} implements reduced rank linear discriminant analysis.
}
  \item{all.n}{A logical specifying whether the mixturegram should plot the profiles of all observations (\code{TRUE}) or just the K-profile summaries (\code{FALSE}).  The default is \code{FALSE}.}
  \item{id.con}{An argument that allows one to impose some sort of (meaningful) identifiability constraint so that the mixture components are in some sort of comparable order between mixture models with different numbers of components.  If \code{NULL}, then the components are ordered by the component means for univariate data or ordered by the first dimension of the component means for multivariate data.}
  \item{score}{The value for the specified dimension reduction technique's score, which is used for constructing the mixturegram.  By default, this value is \code{1}, which is the value that will typically be used.  Larger values will result in more variability displayed on the mixturegram.  Note that the largest value that can be calculated at each value of k>1 on the mixturegram is p+k-1, where p is the number of columns of \code{data}. 
}
  \item{iter.max}{The maximum number of iterations allowed for the k-means clustering algorithm, which is passed to the \code{\link{kmeans}} function.  The default is \code{50}.
}
  \item{nstart}{The number of random sets chosen based on k centers, which is passed to the \code{\link{kmeans}} function.  The default is \code{25}.}
  \item{...}{Additional arguments that can be passed to the underlying \code{\link{plot}} function.}
 }
\value{
  \code{mixturegram} returns a mixturegram where the profiles are plotted over component values of k = 1,...,K.
}
 
\seealso{
\code{\link{boot.comp}}
}
\references{
  Young, D. S., Ke, C., and Zeng, X. (2018) The Mixturegram: A Visualization Tool for Assessing the
  Number of Components in Finite Mixture Models, \emph{Journal 
  of Computational and Graphical Statistics}, \bold{27(3)}, 564--575.
}
\examples{
##Data generated from a 2-component mixture of normals.


set.seed(100)
n <- 100
w <- rmultinom(n,1,c(.3,.7))
y <- sapply(1:n,function(i) w[1,i]*rnorm(1,-6,1) +
            w[2,i]*rnorm(1,0,1))

selection <- function(i,data,rep=30){
  out <- replicate(rep,normalmixEM(data,epsilon=1e-06,
  				   k=i,maxit=5000),simplify=FALSE)
  counts <- lapply(1:rep,function(j) 
                   table(apply(out[[j]]$posterior,1,
                   which.max)))
  counts.length <- sapply(counts, length)
  counts.min <- sapply(counts, min)
  counts.test <- (counts.length != i)|(counts.min < 5)
  if(sum(counts.test) > 0 & sum(counts.test) < rep) 
  	out <- out[!counts.test]
  l <- unlist(lapply(out, function(x) x$loglik))
  tmp <- out[[which.max(l)]]
}

all.out <- lapply(2:5, selection, data = y, rep = 2)
pmbs <- lapply(1:length(all.out), function(i) 
			   all.out[[i]]$post)
mixturegram(y, pmbs, method = "pca", all.n = FALSE,
		    id.con = NULL, score = 1, 
		    main = "Mixturegram (Well-Separated Data)")
}

\keyword{file}
