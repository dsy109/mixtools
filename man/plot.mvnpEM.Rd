\name{plot.mvnpEM}
\title{Plots of Marginal Density Estimates from the mvnpEM Algorithm Output}
\alias{plot.mvnpEM} 
\usage{
\method{plot}{mvnpEM}(x, truenorm = FALSE, lambda = NULL, mu = NULL, v = NULL, 
            lgdcex = 1, \dots)
}
\description{
    Takes an object of class \code{mvnpEM}, as the one returned by the  \code{\link{mvnpEM}}
    algorithm, and returns a set of plots of the 
    density estimates for each coordinate within each multivariate block.
    All the components are displayed on each plot so it is possible
    to see the mixture structure for each coordinate and block. The final bandwidth values are also displayed,
    in a format depending on the bandwidth strategy .     
} 
\arguments{
  \item{x}{An object of class \code{mvnpEM} such as the output
  of the \code{\link{mvnpEM}} function}
  
  \item{truenorm}{Mostly for checking purpose, if the nonparametric model is to be compared
  with a multivariate Gaussian mixture as the true model.}
  
  \item{lambda}{true weight parameters, for Gaussian models only (see above)}
  \item{mu}{true mean parameters, for Gaussian models only (see above)}
  \item{v}{true covariance matrices, for Gaussian models only (see above)}

  \item{lgdcex}{Character expansion factor for \code{\link{legend}}.}
  \item{\dots}{Any remaining arguments are passed to \code{\link{hist}}.}
}
\value{
  \code{plot.mvnpEM} currently just plots the figure.
} 

\seealso{ 
\code{\link{mvnpEM}}, \code{\link{npEM}}, \code{\link{density.npEM}}
}
 

\examples{ 
# example as in Chauveau and Hoang (2015) with 6 coordinates
\dontrun{
m=2; r=6; blockid <-c(1,1,2,2,3,3) # 3 bivariate blocks 
# generate some data x ...
a <- mvnpEM(x, mu0=2, blockid, samebw=F) # adaptive bandwidth
plot(a) # this S3 method produces 6 plots of univariate marginals
summary(a)}
}

\keyword{file}
