library(mixtools)
library(plotly)
library(scales)

plotlyweibullRMM <- function(a, title=NULL, rowstyle=TRUE, subtitle=NULL,
                             width = 3 , col = NULL , 
                             title.size = 15 , title.x = 0.5 , title.y = 0.95,
                             xlab = "Iterations" , xlab.size = 15 , xtick.size = 15,
                             ylab = "Estimates" , ylab.size = 15 , ytick.size = 15,
                             legend.size = 15){
  n <- length(a$x)
  m <- dim(a$all.lambda)[2]
  if (is.null(col)){
    col <- hue_pal()(m)
  }
  if (length(col) != m){
    print("Please specify",m,"colors in 'col'.")
  }
  pcc <- round(100*(1-mean(a$d)),2)
  if (is.null(subtitle)) {
    subtitle <- paste("n=",n,", ", pcc, "% censored", sep="")}
  if (is.null(title)) {
    tt1 <- "Shape parameters"; tt2 <- "Scale parameters" 
    tt3 <- "Weight parameters"
  } else tt1 <- tt2 <- tt3 <- title
  
  plot1 <- plot_ly()%>%
    add_trace(x = seq(from = 1 , to = length(a$all.shape[,1]) , by = 1),
              y = a$all.shape[,1] , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col[1]),
              name = "Shape 1",showlegend = TRUE)%>%
    layout(
      legend = list(font = list(size=legend.size)),
      title = list(text = paste(tt1 , "\n(", subtitle,")"),
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = xlab,
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = ylab,
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (j in 2:m){
    plot1 <- plot1%>%
      add_trace(x = seq(from = 1 , to = length(a$all.shape[,j]) , by = 1),
                y = a$all.shape[,j] , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = col[j]),
                name = paste("Shape",j),showlegend = TRUE)
  }
  
  plot2 <- plot_ly()%>%
    add_trace(x = seq(from = 1 , to = length(a$all.scale[,1]) , by = 1),
              y = a$all.scale[,1] , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col[1]),
              name = "Scale 1",showlegend = TRUE)%>%
    layout(
      legend = list(font = list(size=legend.size)),
      title = list(text = paste(tt2 , "\n(", subtitle,")"),
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = xlab,
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = ylab,
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (j in 2:m){
    plot2 <- plot2%>%
      add_trace(x = seq(from = 1 , to = length(a$all.scale[,j]) , by = 1),
                y = a$all.scale[,j] , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = col[j]),
                name = paste("Scale",j),showlegend = TRUE)
  } 
  
  plot3 <- plot_ly()%>%
    add_trace(x = seq(from = 1 , to = length(a$all.lambda[,1]) , by = 1),
              y = a$all.lambda[,1] , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col[1]),
              name = paste('&#955;',"<sub>",1,"</sub>" , sep=""),showlegend = TRUE)%>%
    layout(
      legend = list(font = list(size=legend.size)),
      title = list(text = paste(tt3 , "\n(", subtitle,")"),
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = xlab,
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = ylab,
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (j in 2:m){
    plot3 <- plot3%>%
      add_trace(x = seq(from = 1 , to = length(a$all.lambda[,j]) , by = 1),
                y = a$all.lambda[,j] , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = col[j]),
                name = paste('&#955;',"<sub>",j,"</sub>" , sep=""),showlegend = TRUE)
  } 
  
  print(plot1)
  print(plot2)
  print(plot3)
}

###### Example ######
n = 500 # sample size
m = 2 # nb components
lambda=c(0.4, 0.6)
shape <- c(0.5,5); scale <- c(1,20) # model parameters
set.seed(321)
x <- rweibullmix(n, lambda, shape, scale) # iid ~ weibull mixture
cs=runif(n,0,max(x)+10) # iid censoring times
t <- apply(cbind(x,cs),1,min) # censored observations
d <- 1*(x <= cs) # censoring indicator
## set arbitrary or "reasonable" (e.g., data-driven) initial values
l0 <- rep(1/m,m); sh0 <- c(1, 2); sc0 <- c(2,10)
# Stochastic EM algorithm
a <- weibullRMM_SEM(t, d, lambda = l0, shape = sh0, scale = sc0, maxit = 200)
summary(a) # Parameters estimates etc
plotlyweibullRMM(a , legend.size = 20) # plot of St-EM sequences
