library(plotly)
library(mixtools)
library(scales)

plotly.mixEM <-function (x, 
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
                         alpha = 0.05, marginal = FALSE){
  def.par <- par(ask=(loglik + density > 1), "mar") # only ask and mar are changed
  mix.object <- x
  if (!inherits(mix.object, "mixEM"))
    stop("Use only with \"mixEM\" objects!")
  ### iteration plot ###
  if (loglik) {
    plot.loglik <- plot_ly()%>%
      add_trace(x= seq(from=1 , to=length(mix.object$all.loglik) , by=1), 
                y= mix.object$all.loglik , type = 'scatter' , mode = 'lines+markers',
                marker = list(color = col1 , size = cex1),
                line = list(color = col1 , width = lwd1),
                name = "Log-Likelihood" , showlegend = FALSE)%>%
      layout(
        title = list(text = title1,
                     x = title1.x,
                     y = title1.y,
                     font = list(size=title1.size)
        ),
        xaxis = list(title = xlab1,
                     tickfont = list(size = xtick1.size),
                     titlefont = list(size = xlab1.size)),
        yaxis = list(title = ylab1,
                     tickfont = list(size = ytick1.size),
                     titlefont = list(size = ylab1.size))
      )
    print(plot.loglik)
  }
  ### density plot ###
  if (density){
    if (mix.object$ft == "logisregmixEM") {
      if (ncol(mix.object$x) != 2) {
        stop("The predictors must have 2 columns!")
      }
      if (sum((mix.object$y == 1) + (mix.object$y == 0)) != length(mix.object$y)) {
        stop("The response must be binary!")
      }
      k = ncol(mix.object$beta)
      x = mix.object$x[, 2]
      if(is.null(title2)) { title2 <- "Most Probable Component Membership" }
      if(is.null(xlab2)) { xlab2 <- "Predictor" }
      if(is.null(ylab2)) { ylab2 <- "Response" }
      if (is.null(col2)){
        col2 <- hue_pal()(k)
      }
      if (length(col2) != k){
        print(paste("Please specify" , k , "colors in 'col2'."))
      }
      plot.density <- plot_ly()%>%
        add_trace(x=x , 
                  y=mix.object$y, 
                  type = 'scatter' , mode = 'markers',
                  marker = list(color = col2[apply(mix.object$posterior,1, which.max)] , size = cex2),
                  name = "Data" , showlegend = FALSE)%>%
        layout(
          title = list(text = title2,
                       x = title2.x,
                       y = title2.y,
                       font = list(size=title2.size)
          ),
          xaxis = list(title = xlab2,
                       tickfont = list(size = xtick2.size),
                       titlefont = list(size = xlab2.size)),
          yaxis = list(title = ylab2,
                       tickfont = list(size = ytick2.size),
                       titlefont = list(size = ylab2.size))
        )
        
      a = cbind(x, mix.object$y)
      a = a[order(a[, 1]), ]
      
      for (i in 1:k) {
        plot.density <- add_trace(plot.density,
                                  x=a[,1] , 
                                  y=plogis(mix.object$beta[1, i]+mix.object$beta[2,i] * a[,1]), 
                                  type = 'scatter' , mode = 'lines',
                                  line = list(width = lwd2 , color = col2[i]),
                                  name = paste("Component" , i) , showlegend = FALSE)
      }
    }
    if (mix.object$ft == "normalmixEM") {
      k <- ncol(mix.object$posterior)
      x <- sort(mix.object$x)
      a <- hist(x, plot = FALSE)
      maxy <- max(max(a$density), 0.3989*mix.object$lambda/mix.object$sigma)
      if(is.null(title2)) { title2 <- "Density Curves" }
      if(is.null(xlab2)) { xlab2 <- "Data" }
      if (is.null(col2)){
        col2 <- hue_pal()(k)
      }
      if (length(col2) != k){
        print(paste("Please specify" , k , "colors in 'col2'."))
      }
      
      plot.density <- plot_ly()%>%
        add_trace(x=x ,
                  type = 'histogram', histnorm = "probability density",
                  name = 'Data' , showlegend = FALSE,
                  marker = list(color = col.hist,
                                line = list(color = col.hist))
                  )%>%
        layout(
          title = list(text = title2,
                       x = title2.x,
                       y = title2.y,
                       font = list(size=title2.size)),
          xaxis = list(title = xlab2,
                       tickfont = list(size = xtick2.size),
                       titlefont = list(size = xlab2.size)),
          yaxis = list(title = ylab2,
                       tickfont = list(size = ytick2.size),
                       titlefont = list(size = ylab2.size),
                       range = c(0 , maxy)),
          bargap = 0.01
        )
      if (length(mix.object$mu) == 1) {
        arbvar <- TRUE
        mix.object$sigma <- mix.object$scale * mix.object$sigma
        arbmean <- FALSE
      }
      if (length(mix.object$mu) == k && length(mix.object$sigma) == 1) {
        arbmean <- TRUE
        arbvar <- FALSE
      }
      if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      for (i in 1:k) {
        plot.density <- add_trace(plot.density,
                                  x=x , 
                                  y=mix.object$lambda[i] *
                                    dnorm(x, mean = mix.object$mu[i * arbmean + (1 - arbmean)],
                                          sd = mix.object$sigma[i * arbvar + (1 - arbvar)]), 
                                  type = 'scatter' , mode = 'lines',
                                  line = list(width = lwd2 , color = col2[i]),
                                  name = paste("Component" , i) , showlegend = FALSE)
      }
    }
    if (mix.object$ft == "repnormmixEM") {
      x <- as.vector(as.matrix(mix.object$x))
      k <- ncol(mix.object$posterior)
      x.sort <- sort(x)
      a <- hist(x.sort, plot = FALSE)
      maxy <- max(max(a$density), .3989*mix.object$lambda/mix.object$sigma)
      if (is.null(title2)) { title2 <- "Density Curves" }
      if(is.null(xlab2)) { xlab2 <- "Data" }
      if (is.null(col2)){
        col2 <- hue_pal()(k)
      }
      if (length(col2) != k){
        print(paste("Please specify" , k , "colors in 'col2'."))
      }
      
      plot.density <- plot_ly()%>%
        add_trace(x=x ,
                  type = 'histogram', histnorm = "probability density",
                  name = 'Data' , showlegend = FALSE,
                  marker = list(color = col.hist,
                                line = list(color = col.hist))
        )%>%
        layout(
          title = list(text = title2,
                       x = title2.x,
                       y = title2.y,
                       font = list(size=title2.size)),
          xaxis = list(title = xlab2,
                       tickfont = list(size = xtick2.size),
                       titlefont = list(size = xlab2.size)),
          yaxis = list(title = ylab2,
                       tickfont = list(size = ytick2.size),
                       titlefont = list(size = ylab2.size),
                       range = c(0 , maxy)),
          bargap = 0.01
        )
      if (length(mix.object$mu) == 1) {
        arbvar <- TRUE
        mix.object$sigma = mix.object$scale * mix.object$sigma
        arbmean <- FALSE
      }
      if (length(mix.object$mu) == k && length(mix.object$sigma) == 1) {
        arbmean <- TRUE
        arbvar <- FALSE
      }
      if (length(mix.object$sigma) == k && length(mix.object$mu) == k) {
        arbmean <- TRUE
        arbvar <- TRUE
      }
      for (i in 1:k) {
        plot.density <- add_trace(plot.density,
                                  x=x.sort , 
                                  y=mix.object$lambda[i] * 
                                    dnorm(x.sort, mean = mix.object$mu[i * arbmean + (1 - arbmean)], 
                                          sd = mix.object$sigma[i * arbvar + (1 - arbvar)]), 
                                  type = 'scatter' , mode = 'lines',
                                  line = list(width = lwd2 , color = col2[i]),
                                  name = paste("Component" , i) , showlegend = FALSE)
      }
    }
    if (mix.object$ft == "regmixEM.mixed") {
      if (is.null(col2)){
        col2 <- hue_pal()(ncol(x$posterior.z))
      }
      if (length(col2) != ncol(x$posterior.z)){
        print(paste("Please specify", ncol(x$posterior.z) ,"color in 'col2'."))
      }
      x.1 = mix.object$x
      n = sum(sapply(x.1, nrow))
      x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][,1])))
      if (x.1.sum == n) {
        x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][,-1], ncol = 1))
      }else {
        x = x.1
      }
      plot.density <- plotly.post.beta(x = x, y = mix.object$y, p.beta = mix.object$posterior.beta, 
                                       p.z = mix.object$posterior.z ,
                                       cex = cex2,lwd=lwd2,
                                       title.size = title2.size,
                                       xlab.size = xlab2.size , xtick.size = xtick2.size,
                                       ylab.size = ylab2.size , ytick.size = ytick2.size,
                                       col.comp = col2) 
    }
    if (mix.object$ft == "mvnormalmixEM") {
      x = mix.object$x
      if (ncol(x) != 2) {
        stop("The data must have 2 columns!")
      }
      post = apply(mix.object$posterior, 1, which.max)
      k <- ncol(mix.object$posterior)
      if (is.null(col2)){
        col2 <- hue_pal()(k)
      }
      if (length(col2) != k){
        print(paste("Please specify" ,k," colors in 'col2'."))
      }
      if (is.list(mix.object$sigma)) {
        sigma = mix.object$sigma
      } else {
        sigma = lapply(1:k, function(i) mix.object$sigma)
      }
      if (is.list(mix.object$mu)) {
        mu = mix.object$mu
      } else {
        mu = lapply(1:k, function(i) mix.object$mu)
      }
      if(is.null(xlab2)) { xlab2 <- "X.1" }
      if(is.null(ylab2)) { ylab2 <- "X.2" }
      if (!marginal) {
        if (is.null(title2)) { title2 <- "Density Curves" }
      }
      if (marginal) {
        title2 <- ""
      }
      
      plot.main <- plot_ly()%>%
        add_trace(x=mix.object$x[,1] , 
                  y=mix.object$x[,2], 
                  type = 'scatter' , mode = 'markers',
                  marker = list(color = col2[post] , size = cex2),
                  name = "Data" , showlegend = FALSE)%>%
        layout(
          title = list(text = title2,
                       x = title2.x,
                       y = title2.y,
                       font = list(size=title2.size)
          ),
          xaxis = list(title = xlab2,
                       tickfont = list(size = xtick2.size),
                       titlefont = list(size = xlab2.size)),
          yaxis = list(title = ylab2,
                       tickfont = list(size = ytick2.size),
                       titlefont = list(size = ylab2.size))
        )
      for (i in 1:k){
        plot.main <- add_markers(plot.main,
                                 x = mu[[i]][1],
                                 y = mu[[i]][2],
                                 marker = list(color = "black" , size = cex2+3),
                                 name = paste("Center" , i) , showlegend = FALSE)
      }
      es.multi <- lapply(sigma, eigen)
      e1.multi <- lapply(es.multi, function(x){x$vectors%*%diag(sqrt(x$values))})
      r1.multi <- sapply(alpha, function(x){sqrt(qchisq(1-x,2))})
      theta <- seq(0,2*pi,len=300)
      v1.multi <- lapply(r1.multi , function(x){cbind(x*cos(theta),x*sin(theta))})
      pts.multi <- rep(list(NA),length(sigma))
      for (i in 1:length(sigma)){
        pts.multi[[i]] <- rep(list(NA) , length(alpha))
        for (j in 1:length(alpha)){
          pts.multi[[i]][[j]] <- t(mu[[i]]-e1.multi[[i]]%*%t(v1.multi[[j]]))
        }
      }
      for (i in 1:k) {
        for (j in 1:length(alpha)) {
          plot.main <- add_trace(
            plot.main,
            x=pts.multi[[i]][[j]][,1] , 
            y=pts.multi[[i]][[j]][,2] , type = 'scatter' , mode = 'lines',
            line = list(color = col2[i] , width = lwd2),
            name = paste((1-alpha[j])*100,'% Ellipse'),showlegend = FALSE)
        }
      }  
      if (!marginal){
        plot.density <- plot.main
      }
      if (marginal){
        x.marginal <- plot_ly()%>%
          add_trace(x=mix.object$x[, 1],
                    type = 'histogram',
                    name = "Dist X",
                    showlegend = FALSE,
                    marker = list(color = col.hist,
                                  line = list(color = col.hist))
          )%>%
          layout(
            bargap = 0.01
          )
        
        y.marginal <- plot_ly()%>%
          add_trace(y=mix.object$x[, 2],
                    type = 'histogram',
                    name = "Dist X",
                    showlegend = FALSE,
                    marker = list(color = col.hist,
                                  line = list(color = col.hist))
          )%>%
          layout(
            bargap = 0.01
          )
        
        plot.density <- subplot(
          x.marginal,
          plotly_empty(type = 'scatter' , mode = "markers"),
          plot.main,
          y.marginal,
          nrows = 2, heights = c(.2, .8), widths = c(.8,.2), margin = 0,
          shareX = TRUE, shareY = TRUE) %>%
          layout(showlegend = F)
      }
    }
    print(plot.density)
    
    if (mix.object$ft == "expRMM_EM") {plotexpRMM(mix.object, ...)} # all default
    if (mix.object$ft == "weibullRMM_SEM") {plotweibullRMM(mix.object, ...)} # all default
  }
  par(def.par) # reset ask and mar to original values
}

#####################
###### Example ######
#####################
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
plotly.mixEM(out.2 , col2 = c("red" , "green") , density = TRUE)

##Fitting randomly generated data with a 2-component location mixture of bivariate normals.
set.seed(100)
x.1 <- rmvnorm(40, c(0, 0))
x.2 <- rmvnorm(60, c(3, 4))
X.1 <- rbind(x.1, x.2)
mu <- list(c(0, 0), c(3, 4))
out.1 <- mvnormalmixEM(X.1, arbvar = FALSE, mu = mu,
                       epsilon = 1e-02)
plotly.mixEM(out.1 , col2 = c("brown" , "blue") ,
             alpha = c(0.01 , 0.05 , 0.1),
             density = TRUE , marginal = FALSE)

##Fitting randomly generated data with a 2-component scale mixture of bivariate normals.
x.3 <- rmvnorm(40, c(0, 0), sigma =
                 matrix(c(200, 1, 1, 150), 2, 2))
x.4 <- rmvnorm(60, c(0, 0))
X.2 <- rbind(x.3, x.4)
lambda <- c(0.40, 0.60)
sigma <- list(diag(1, 2), matrix(c(200, 1, 1, 150), 2, 2))
out.2 <- mvnormalmixEM(X.2, arbmean = FALSE,
                       sigma = sigma, lambda = lambda,
                       epsilon = 1e-02)
plotly.mixEM(out.1 , col2 = c("brown" , "blue") ,
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
plotly.mixEM(em.out , col2 = c("gold" , "purple") , 
             density = TRUE , lwd2 = 1 , cex2 =9)

## Analyzing the Old Faithful geyser data with a 2-component mixture of normals.
data(faithful)
attach(faithful)
set.seed(100)
out <- normalmixEM(waiting, arbvar = FALSE, verb = TRUE,
                   epsilon = 1e-04)
plotly.mixEM(out, density = TRUE , col2 = c("gold" , "purple"))

## EM output for the water-level task data set.
data(Waterdata)
set.seed(100)
water <- t(as.matrix(Waterdata[,3:10]))
em.out <- repnormmixEM(water, k = 2, verb = TRUE, epsilon = 1e-03)
plotly.mixEM(em.out, density = TRUE , col2 = c("gold" , "purple"))
