library(mixtools)
library(plotly)
library(scales)

plotly_expRMM <- function(a , title = NULL , 
                          rowstyle = TRUE , subtitle=NULL,
                          width = 2 , cex = 2 , col.comp = NULL,
                          legend.text = NULL,
                          legend.text.size = 15,
                          legend.size = 15,
                          title.x = 0.5,
                          title.y = 0.95,
                          title.size = 15,
                          xlab.size = 15,
                          xtick.size = 15,
                          ylab.size = 15,
                          ytick.size = 15
){
  n <- length(a$x) 
  m <- dim(a$all.lambda)[2]
  if (is.null(col.comp)){
    col.comp <- hue_pal()(m)
  }
  if (length(col.comp) != m){
    print(paste("Please specify",m,"colors in 'col.comp'."))
  }
  pcc <- round(100*(1-mean(a$d)),2)
  sizes <- paste("n=",n,", ", pcc, "% censored", sep="")
  if (is.null(subtitle)) {
    subtitle <- paste("n=",n,", ", pcc, "% censored", sep="")}
  if (is.null(title)) {
    tt1 <- "Rate parameters" 
    tt2 <- "Weight parameters"
  } else {tt1 <- tt2 <- title}
  # lgd <- expression(xi[1]) ## Cannot be converted to plotly.
  
  plot1 <- plot_ly() %>%
    layout(
      legend = list(title=legend.text,
                    font = list(size=legend.size)),
      xaxis = list(title = "Iterations",
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = "Estimates",
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (j in 1:m){
    plot1 <- plot1 %>%
      add_trace(x=seq(from = 1 , to = length(a$all.rate[,j]) , by = 1) , 
                y=a$all.rate[,j] , type = 'scatter' , mode = 'lines+markers',
                marker = list(size = cex , color = col.comp[j]),
                line = list(width = width , color = col.comp[j]),
                name = paste("Comp" , j),showlegend = TRUE)
  }
  
  plot2 <- plot_ly() %>%
    layout(
      legend = list(title=legend.text,
                    font = list(size=legend.size)),
      xaxis = list(title = "Iterations",
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = "Estimates",
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (j in 1:m){
    plot2 <- plot2 %>%
      add_trace(x=seq(from = 1 , to = length(a$all.lambda[,j]) , by = 1) , 
                y=a$all.lambda[,j] , type = 'scatter' , mode = 'lines+markers',
                marker = list(size = cex , color = col.comp[j]),
                line = list(width = width , color = col.comp[j]),
                name = paste("Comp" , j),showlegend = FALSE)
  }
  
  if (rowstyle) {
    nrow <- 1
    x.1 <- 0.25
    x.2 <- 0.75
    y.1 <- y.2 <- 0.95
    share.X <- FALSE
  } else if (!rowstyle){
    nrow <- 2
    x.1 <- x.2 <- 0.5
    y.1 <- 0.95
    y.2 <- 0.45
    share.X <- TRUE
  }
  plot.all <- subplot(
    plot1 , plot2 , nrows = nrow , shareX = share.X , shareY = FALSE,
    titleX = TRUE , titleY = TRUE
  )%>%
    layout(annotations = list(
      list(
        x = x.1, 
        y = y.1, 
        font = list(size = title.size), 
        text = paste(tt1 , "\n(", subtitle,")"), 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      ), 
      list(
        x = x.2, 
        y = y.2, 
        font = list(size = title.size), 
        text = paste(tt2 , "\n(", subtitle,")"), 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      )
    )
  )
  print(plot.all)
}

###### Example ######
n=300 # sample size
m=2 # number of mixture components
lambda <- c(1/3,1-1/3); rate <- c(1,1/10) # mixture parameters
set.seed(1234)
x <- rexpmix(n, lambda, rate) # iid ~ exponential mixture
cs=runif(n,0,max(x)) # Censoring (uniform) and incomplete data
t <- apply(cbind(x,cs),1,min) # observed or censored data
d <- 1*(x <= cs) # censoring indicator
###### EM for RMM, exponential lifetimes
l0 <- rep(1/m,m); r0 <- c(1, 0.5) # "arbitrary" initial values
a <- expRMM_EM(t, d, lambda=l0, rate=r0, k = m)
summary(a) # EM estimates etc
plotly_expRMM(a , rowstyle = TRUE) # plot of EM sequences
