\name{makemultdata}
\title{Produce Cutpoint Multinomial Data}
\alias{makemultdata}
\usage{
makemultdata(..., cuts)
}

\description{
  Change data into a matrix of multinomial counts using the
  cutpoint method and generate EM algorithm starting values for 
  a k-component mixture of multinomials.
}
\arguments{
  \item{...}{Either vectors (possibly of different lengths) of raw data 
    or an nxm matrix (or data frame) of data. If \code{...} are vectors of varying length, 
    then \code{makemultdata} will create a matrix of size nxm where n is the
    sample size and m is the length of the vector with maximum length.  Those 
    vectors with length less than m will have \code{NA}s to make the 
    corresponding row in the matrix of length m.  If \code{...} is a matrix (or data frame), then
    the rows must correspond to the sample and the columns the repeated measures.}
  \item{cuts}{A vector of cutpoints.  This vector is sorted by the algorithm.}
}
\value{
  \code{makemultdata} returns an object which is a list with components:
  \item{x}{An nxm matrix of the raw data.}
  \item{y}{An nxp matrix of the discretized data where p is one more than the
    number of cutpoints. Each row is a multinomial vector of counts.  In particular,
    each row should sum to the number of repeated measures for that sample.}
} 
\details{
  The (i, j)th entry of the matrix \code{y} (for j < p)
  is equal to the number of entries
  in the ith column of \code{x} that are less than or equal to \code{cuts}[j].
  The (i, p)th entry is equal to the number of entries greater than
  \code{cuts}[j].
}
\seealso{
\code{\link{compCDF}}, \code{\link{multmixmodel.sel}}, \code{\link{multmixEM}}
}
\references{
  Elmore, R. T., Hettmansperger, T. P. and Xuan, F. (2004) The Sign Statistic, One-Way Layouts
  and Mixture Models, \emph{Statistical Science} \bold{19(4)}, 579--587.
}
\examples{
## Randomly generated data.

set.seed(100)
y <- matrix(rpois(70, 6), 10, 7)
cuts <- c(2, 5, 7)
out1 <- makemultdata(y, cuts = cuts)
out1

## The sulfur content of the coal seams in Texas.

A <- c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B <- c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59)
C <- c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49)
D <- c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32)
E <- c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)

out2 <- makemultdata(A, B, C, D, E, 
                     cuts = median(c(A, B, C, D, E)))
out2

## The reaction time data.

data(RTdata)
out3 <- makemultdata(RTdata, cuts = 
                     100*c(5, 10, 12, 14, 16, 20, 25, 30, 40, 50))
dim(out3$y)
out3$y[1:10,]
}

\keyword{file}
