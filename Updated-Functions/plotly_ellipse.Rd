\name{plotly_ellipse}
\title{Draw Two-Dimensional Ellipse Based on Mean and Covariance using \code{plotly}}
\alias{plotly_ellipse}
\usage{
plotly_ellipse(mu, sigma, alpha=.05, npoints=250,
               draw=TRUE, cex = 3, col = "#1f77b4", lwd = 3,
               title = "", title.x = 0.5, title.y = 0.95, title.size = 15,
               xlab = "X", xlab.size = 15, xtick.size = 15,
               ylab = "Y", ylab.size = 15, ytick.size = 15)
}

\description{
This is an updated version of \code{ellipse}. For more technical details, please refer to \code{ellipse}.
}
\arguments{
  \item{mu}{A 2-vector giving the mean.}
  \item{sigma}{A 2x2 matrix giving the covariance matrix.}
  \item{alpha}{Probability to be excluded from the ellipse. The
default value is alpha = .05, which results in a 95\% ellipse.}
  \item{npoints}{Number of points comprising the border of the ellipse.}
  \item{draw}{If TRUE, draw the ellipse.}
  \item{cex}{Size of markers.}
  \item{lwd}{Line width of the ellipse.}
  \item{col}{Color of both markers and lines.}
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
}
\value{
  \code{plotly_ellipse} returns an \code{npoints}x2 matrix of the points forming the
border of the ellipse.
}
\references{
  Johnson, R. A. and Wichern, D. W. (2002) \emph{Applied Multivariate Statistical Analysis, Fifth Edition},
  Prentice Hall.
}
\seealso{
\code{\link{regcr}}, \code{\link{ellipse}}
}
\examples{
## Produce a 95% ellipse with the specified mean and covariance structure.
mu <- c(1, 3)
sigma <- matrix(c(1, .3, .3, 1.5), 2, 2)
plotly_ellipse(mu, sigma, npoints = 200)
}

\keyword{file}
