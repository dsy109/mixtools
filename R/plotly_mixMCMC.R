plotly_mixMCMC <- function(
  x , trace.plot = TRUE , summary.plot = FALSE , burnin = 2000, credit.region = 0.95, col.cr = NULL,
  cex.trace = 3 , width.trace = 3 , 
  cex.summary = 3 , width.summary = 1,
  title.trace = "", title.trace.x = 0.5 , title.trace.y = 0.95, title.trace.size = 15,
  xlab.trace = "Index" , xlab.trace.size = 15, xtick.trace.size = 15,
  ylab.trace = NULL , ylab.trace.size = 15, ytick.trace.size = 15,
  title.summary = "Credible Regions", title.summary.x = 0.5 , title.summary.y = 0.95, title.summary.size = 15,
  xlab.summary = "Predictor" , xlab.summary.size = 15, xtick.summary.size = 15,
  ylab.summary = "Response" , ylab.summary.size = 15, ytick.summary.size = 15
){
  mix.object <- x
  
  if (!inherits(mix.object, "mixMCMC")){
    stop("Use only with \"mixMCMC\" objects!")
  }
  if (trace.plot){
    k<-mix.object$components
    theta<-mix.object$theta
    p.k=ncol(theta)
    p=p.k/k
    name.theta<-colnames(theta)
    
    for (i in 1:dim(theta)[2]){
      if (is.null(ylab.trace)){
        ylab.trace <- name.theta[i]
      }
      plot.trace <- plot_ly()%>%
        add_trace(x=seq(from = 1 , to = dim(theta)[1] , by=1) , 
                  y=theta[,i] , type = 'scatter' , mode = 'lines+markers',
                  marker = list(size = cex.trace),
                  line = list(width = width.trace),
                  showlegend = FALSE) %>%
        plotly::layout(
          title = list(text = title.trace,
                       x = title.trace.x,
                       y = title.trace.y,
                       font = list(size=title.trace.size)),
          xaxis = list(title = list(text = xlab.trace,
                                    font = list(size = xlab.trace.size)),
                       tickfont = list(size = xtick.trace.size)),
          yaxis = list(title = list(text = ylab.trace,
                                    font = list(size = ylab.trace.size)),
                       tickfont = list(size = ytick.trace.size))
          )
      print(plot.trace)
    }
  }
  if(is.matrix(mix.object$x) == TRUE && is.null(mix.object$y) == FALSE && summary.plot == TRUE){
    y<-mix.object$y
    n<-length(y)
    x<-mix.object$x
    p<-ncol(x)
    k<-mix.object$components
    theta<-mix.object$theta
    if(p!=2 || sum(x[,1])!=n){                                                                  
      stop(paste("This only works for simple linear regression!","\n"))
    }
    plot.summary <- plot_ly()%>%
      add_trace(x=x[,2] , 
                y=y , type = 'scatter' , mode = 'markers',
                marker = list(size = cex.summary),
                showlegend = FALSE) %>%
      plotly::layout(
        title = list(text = paste(credit.region*100, "% " , title.summary ,sep = ""),
                     x = title.summary.x,
                     y = title.summary.y,
                     font = list(size=title.summary.size)),
        xaxis = list(title = list(text = xlab.summary,
                                  font = list(size = xlab.summary.size)),
                     tickfont = list(size = xtick.summary.size)
        ),
        yaxis = list(title = list(text = ylab.summary,
                                  font = list(size = ylab.summary.size)),
                     tickfont = list(size = ytick.summary.size)
        )
      )
    if(is.null(col.cr)){
      col.cr <- hue_pal()(k)
    }
    if (length(col.cr) != k){
      print(paste("Please specify",k,"colors in 'col.cr'."))
    }
    for (i in 1:k){
      beta.summary <- cbind(theta[-c(1:burnin),2*i-1],theta[-c(1:burnin),2*i])
      xbar = apply(beta.summary, 2, mean)
      n = nrow(beta.summary)
      cbeta = t(t(beta.summary) - xbar)
      S = t(cbeta) %*% cbeta/(n - 1)
      eS = eigen(S)
      B = eS$vec %*% diag(sqrt(eS$val))
      theta.summary = seq(0, 2 * pi, len = 250)
      alpha = 1-credit.region
      v = cbind(cos(theta.summary), sin(theta.summary)) * sqrt(qchisq(1 - alpha,2))
      h = t(B %*% t(v) + xbar)
      nh = nrow(h)
      m = which.max(h[, 2])
      h = rbind(h[m:nh, ], h[((1:nh) < m), ], h[m, ])
      bound = h
      x.summary <- x[,2]
      z <- length(x.summary) * 5
      u <- seq(min(x.summary), max(x.summary), length = z)
      lower <- c()
      upper <- c()
      v <- c()
      for (q in 1:z) {
        for (l in 1:nrow(beta.summary)) {
          v[l] <- as.matrix(beta.summary[, 1][l] + beta.summary[, 2][l] * u[q])
        }
        uv <- cbind(u[q], v)
        lower <- rbind(lower, uv[order(v), ][1, ])
        upper <- rbind(upper, uv[order(v), ][nrow(beta.summary), ])
      }
      plot.summary <- plot.summary%>%
        add_trace(x=lower[,1] , 
                  y=lower[,2] , type = 'scatter' , mode = 'lines',
                  line = list(width = width.summary , color = col.cr[i]),
                  name = paste("Lower Bound for Component",i), showlegend = FALSE)%>%
        add_trace(x=upper[,1] , 
                  y=upper[,2] , type = 'scatter' , mode = 'lines',
                  line = list(width = width.summary , color = col.cr[i]),
                  name = paste("Upper Bound for Component",i), showlegend = FALSE)
    }
    print(plot.summary)
  }
}