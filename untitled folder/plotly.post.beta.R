plotly.post.beta <- function(y, x, p.beta, p.z,
                             cex = 6,lwd=1,
                             title.size = 15,
                             xlab.size = 15 , xtick.size = 15,
                             ylab.size = 15 , ytick.size = 15,
                             col.data = "#1f77b4",
                             col.comp = c("#f54242" , "#1ae865")){
  N = length(y)
  k = ncol(p.z)
  g = apply(p.z, 1, function(i) (1:length(i))[i == max(i)])
  
  if (length(col.comp) != k){
    print(paste("Please specify" , k , "colors in 'col.comp'."))
  }
  
  abline.fun <- function(x , intercept , slope){
    y <- slope * x + intercept
    return(y)
  }
  
  
  min.int=min(sapply(1:N, function(i) min(p.beta[[i]][1,])))
  max.int=max(sapply(1:N, function(i) max(p.beta[[i]][1,])))
  min.slp=min(sapply(1:N, function(i) min(p.beta[[i]][2,])))
  max.slp=max(sapply(1:N, function(i) max(p.beta[[i]][2,])))
  
  # All posterior regression lines.
  fig1 <- plot_ly()%>%
    layout(
      xaxis = list(title = "x-values",
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = "y-values",
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (i in 1:length(x)){
    fig1 <- add_trace(fig1,
      x=as.vector(x[[i]]) , 
      y=as.vector(y[[i]]), 
      type = 'scatter' , mode = 'markers',
      marker = list(size = cex , color = col.data),
      name = paste("Comp" , i) , showlegend = FALSE)
  }
  for (j in 1:k){
    for (l in 1:length(p.beta)){
      fig1 <- add_trace(fig1,
                        x=c(min(unlist(x)) , max(unlist(x))) , 
                        y=c(abline.fun(x=min(unlist(x)) , 
                                       intercept = p.beta[[l]][1,j] , 
                                       slope = p.beta[[l]][2,j]),
                            abline.fun(x=max(unlist(x)) , 
                                       intercept = p.beta[[l]][1,j] , 
                                       slope = p.beta[[l]][2,j])), 
                        type = 'scatter' , mode = 'lines',
                        line = list(width = lwd , color = col.comp[j]),
                        name = paste("Fit" , j) , showlegend = FALSE)
    }
  }
  
  #Posterior regression lines chosen according to the membership probabilities. 
  fig2 <- plot_ly()%>%
    layout(
      xaxis = list(title = "x-values",
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = "y-values",
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (i in 1:length(x)){
    fig2 <- add_trace(fig2,
                      x=as.vector(x[[i]]) , 
                      y=as.vector(y[[i]]), 
                      type = 'scatter' , mode = 'markers',
                      marker = list(size = cex , color = col.data),
                      name = paste("Comp" , i) , showlegend = FALSE)
  }
  for (l in 1:length(p.beta)){
    fig2 <- add_trace(fig2,
                      x=c(min(unlist(x)) , max(unlist(x))) , 
                      y=c(abline.fun(x=min(unlist(x)) , 
                                     intercept = p.beta[[l]][1,g[l]] , 
                                     slope = p.beta[[l]][2,g[l]]),
                          abline.fun(x=max(unlist(x)) , 
                                     intercept = p.beta[[l]][1,g[l]] , 
                                     slope = p.beta[[l]][2,g[l]])), 
                      type = 'scatter' , mode = 'lines',
                      line = list(width = lwd , color = col.comp[g[l]]),
                      name = paste("Fit" , j) , showlegend = FALSE)
  }
  
  #All posterior beta values.
  fig3 <- plot_ly()%>%
    layout(
      xaxis = list(title = "Posterior Intercepts",
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = "Posterior Slopes",
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (j in 1:k){
    for (l in 1:length(p.beta)){
      fig3 <- add_trace(fig3,
                        x=p.beta[[l]][1,j], 
                        y=p.beta[[l]][2,j], 
                        type = 'scatter' , mode = 'markers',
                        marker = list(size = cex , color = col.comp[j]),
                        name = paste("Comp" , j) , showlegend = FALSE)
    }
  }
  #Posterior beta values chosen according to the membership probabilities.
  fig4 <- plot_ly()%>%
    layout(
      xaxis = list(title = "Posterior Intercepts",
                   tickfont = list(size = xtick.size),
                   titlefont = list(size = xlab.size)),
      yaxis = list(title = "Posterior Slopes",
                   tickfont = list(size = ytick.size),
                   titlefont = list(size = ylab.size))
    )
  for (l in 1:length(p.beta)){
    fig4 <- add_trace(fig4,
                      x=p.beta[[l]][1,g[l]] , 
                      y=p.beta[[l]][2,g[l]], 
                      type = 'scatter' , mode = 'markers',
                      marker = list(size = cex , color = col.comp[g[l]]),
                      name = paste("Comp" , l) , showlegend = FALSE)
  }
  
  fig <- subplot(fig1, fig2, fig3 , fig4 , nrows = 2,
                 titleX = TRUE , titleY = TRUE,
                 margin = c(0.03,0.03,0.15,0.15)) %>%
    layout(annotations = list(
      list(
        x = 0.225, 
        y = 1.0, 
        font = list(size = title.size), 
        text = "Data and Posterior Regression Lines", 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      ), 
      list(
        x = 0.775, 
        y = 1.0, 
        font = list(size = title.size), 
        text = "Data and Most Probable Posterior Regression Lines", 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      ), 
      list(
        x = 0.225, 
        y = 0.375, 
        font = list(size = title.size), 
        text = "All Posterior Regression Coefficients", 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      ), 
      list(
        x = 0.775, 
        y = 0.375, 
        font = list(size = title.size), 
        text = "Most Probable Posterior Regression Coefficients", 
        xref = "paper", 
        yref = "paper", 
        xanchor = "center", 
        yanchor = "bottom", 
        showarrow = FALSE
      )
    )
  )
  print(fig)
}

### Example ###
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

x.1 = em.out$x
n = sum(sapply(x.1, nrow))
x.1.sum = sum(sapply(1:length(x.1), function(i) length(x.1[[i]][,1])))
if (x.1.sum == n) {
  x = lapply(1:length(x.1), function(i) matrix(x.1[[i]][,-1], ncol = 1))
} else {
  x = x.1
}

plotly.post.beta(x = x, y = em.out$y, p.beta = em.out$posterior.beta, 
                 p.z = em.out$posterior.z , 
                 col.comp = c("brown" , "gold") , col.data = "black")
