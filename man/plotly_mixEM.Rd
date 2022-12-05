\name{plotly_mixEM}
\alias{plotly_mixEM}

\title{Visualization of output of \code{mixEM} function using \code{plotly}}

\description{This is an updated version of \code{plot.mixEM}. For more technical details, please refer to \code{plot.mixEM}.}

\usage{
plotly_mixEM(x, 
             loglik = TRUE,
             density = FALSE,
             xlab1="Iteration", xlab1.size=15 , xtick1.size=15,
             ylab1="Log-Likelihood", ylab1.size=15 , ytick1.size=15,
             title1="Observed Data Log-Likelihood", title1.size=15,
             title1.x = 0.5,title1.y=0.95,
             col1="#1f77b4", lwd1=3, cex1=6,
             xlab2=NULL, xlab2.size=15 , xtick2.size=15,
             ylab2=NULL, ylab2.size=15 , ytick2.size=15,
             title2=NULL, title2.size=15,
             title2.x = 0.5,title2.y=0.95, col.hist = "#1f77b4",
             col2=NULL, lwd2=3, cex2=6,
             alpha = 0.05, marginal = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{An object of class \code{mixEM}.}
  \item{loglik}{If TRUE, a plot of the log-likelihood versus the EM iterations is given.}
  \item{density}{Graphics pertaining to certain mixture models. The details are given below.}
  \item{xlab1}{Label of x-axis to be passed to the loglikelihood plot. Trying to change these parameters using \code{xlab} result in an error.} 
  \item{xlab1.size}{Font of \code{xlab1}.}
  \item{xtick1.size}{Font of tick labels of x-axis to be passed to the loglikelihood plot.}
  \item{ylab1}{Label of y-axis to be passed to the loglikelihood plot. Trying to change these parameters using \code{ylab} result in an error.} 
  \item{ylab1.size}{Font of \code{ylab1}.}
  \item{ytick1.size}{Font of tick labels of y-axis to be passed to the loglikelihood plot.}
  \item{title1}{Title to be passed to the loglikelihood plot.} 
  \item{title1.size}{Tile size of the loglikelihood plot.}
  \item{title1.x}{Horizontal position of the loglikelihood plot.}
  \item{title1.y}{Verticle position of the loglikelihood plot.}
  \item{col1}{Color of the loglikelihood plot.}
  \item{lwd1}{Width of the density curve of the loglikelihood plot.}
  \item{cex1}{Dot size of the loglikelihood plot.}
  
  \item{xlab2}{Label of x-axis to be passed to the density plot. Trying to change these parameters using \code{xlab} result in an error.} 
  \item{xlab2.size}{Font of \code{xlab2}.}
  \item{xtick2.size}{Font of tick labels of x-axis to be passed to the density plot.}
  \item{ylab2}{Label of y-axis to be passed to the density plot. Trying to change these parameters using \code{ylab} result in an error.} 
  \item{ylab2.size}{Font of \code{ylab2}.}
  \item{ytick2.size}{Font of tick labels of y-axis to be passed to the density plot.}
  \item{title2}{Title to be passed to the density plot.} 
  \item{title2.size}{Tile size of the density plot.}
  \item{title2.x}{Horizontal position of the density plot.}
  \item{title2.y}{Verticle position of the density plot.}
  \item{col2}{Color of the density plot.}
  \item{lwd2}{Width of the density curve of the density plot.}
  \item{cex2}{Dot size of the density plot.}
  \item{col.hist}{Color of the histogram of the density plot}
  \item{alpha}{A vector of significance levels when constructing confidence ellipses and confidence bands for the mixture of multivariate normals and mixture of regressions
cases, respectively. The default is 0.05} 
  \item{marginal}{If \code{TRUE}, marginal density is presented on the side of the corresponding variable.}
 }

\value{A plot of the output of \code{mixEM} function is presented depends on output type.}

\seealso{
\code{\link{post.beta}}}

\examples{
\dontrun{
## EM output for data generated from a 2-component binary logistic regression model.
beta <- matrix(c(-10, .1, 20, -.1), 2, 2)
x <- runif(500, 50, 250)
x1 <- cbind(1, x)
xbeta <- x1%*%beta
w <- rbinom(500, 1, .3)
y <- w*rbinom(500, size = 1, prob = (1/(1+exp(-xbeta[, 1]))))+
  (1-w)*rbinom(500, size = 1, prob =
                 (1/(1+exp(-xbeta[, 2]))))
out.2 <- logisregmixEM(y, x, beta = beta, lambda = c(.3, .7),
                       verb = TRUE, epsilon = 1e-01)
plotly_mixEM(out.2 , col2 = c("red" , "green") , density = TRUE)

## Fitting randomly generated data with a 2-component location mixture of bivariate normals.
set.seed(100)
x.1 <- rmvnorm(40, c(0, 0))
x.2 <- rmvnorm(60, c(3, 4))
X.1 <- rbind(x.1, x.2)
mu <- list(c(0, 0), c(3, 4))
out.1 <- mvnormalmixEM(X.1, arbvar = FALSE, mu = mu,
                       epsilon = 1e-02)
plotly_mixEM(out.1 , col2 = c("brown" , "blue") ,
             alpha = c(0.01 , 0.05 , 0.1),
             density = TRUE , marginal = FALSE)

## Fitting randomly generated data with a 2-component scale mixture of bivariate normals.
x.3 <- rmvnorm(40, c(0, 0), sigma =
                 matrix(c(200, 1, 1, 150), 2, 2))
x.4 <- rmvnorm(60, c(0, 0))
X.2 <- rbind(x.3, x.4)
lambda <- c(0.40, 0.60)
sigma <- list(diag(1, 2), matrix(c(200, 1, 1, 150), 2, 2))
out.2 <- mvnormalmixEM(X.2, arbmean = FALSE,
                       sigma = sigma, lambda = lambda,
                       epsilon = 1e-02)
plotly_mixEM(out.1 , col2 = c("brown" , "blue") ,
             alpha = c(0.01 , 0.05 , 0.1),
             density = TRUE , marginal = TRUE)

## EM output for simulated data from 2-component mixture of random effects.
data(RanEffdata)
set.seed(100)
x <- lapply(1:length(RanEffdata), function(i)
  matrix(RanEffdata[[i]][, 2:3], ncol = 2))
x <- x[1:20]
y <- lapply(1:length(RanEffdata), function(i)
  matrix(RanEffdata[[i]][, 1], ncol = 1))
y <- y[1:20]
lambda <- c(0.45, 0.55)
mu <- matrix(c(0, 4, 100, 12), 2, 2)
sigma <- 2
R <- list(diag(1, 2), diag(1, 2))
em.out <- regmixEM.mixed(y, x, sigma = sigma, arb.sigma = FALSE,
                         lambda = lambda, mu = mu, R = R,
                         addintercept.random = FALSE,
                         epsilon = 1e-02, verb = TRUE)
plotly_mixEM(em.out , col2 = c("gold" , "purple") , 
             density = TRUE , lwd2 = 1 , cex2 =9)

## Analyzing the Old Faithful geyser data with a 2-component mixture of normals.
data(faithful)
attach(faithful)
set.seed(100)
out <- normalmixEM(waiting, arbvar = FALSE, verb = TRUE,
                   epsilon = 1e-04)
plotly_mixEM(out, density = TRUE , col2 = c("gold" , "purple"))

## EM output for the water-level task data set.
data(Waterdata)
set.seed(100)
water <- t(as.matrix(Waterdata[,3:10]))
em.out <- repnormmixEM(water, k = 2, verb = TRUE, epsilon = 1e-03)
plotly_mixEM(em.out, density = TRUE , col2 = c("gold" , "purple"))
}
}
%
\keyword{file}