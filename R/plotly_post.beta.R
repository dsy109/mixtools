plotly_post.beta <- function(y, x, p.beta, p.z,
                             cex = 6,lwd=1,
                             title.size = 15,
                             xlab.size = 15 , xtick.size = 15,
                             ylab.size = 15 , ytick.size = 15,
                             col.data = "#1f77b4",
                             col.comp = NULL){
  N = length(y)
  k = ncol(p.z)
  g = apply(p.z, 1, function(i) (1:length(i))[i == max(i)])
  
  if(is.null(col.comp)){
    col.comp <- hue_pal()(k)
  }
  
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
    plotly::layout(
      xaxis = list(title = list(text = "x-values",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = "y-values",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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
    plotly::layout(
      xaxis = list(title = list(text = "x-values",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = "y-values",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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
    plotly::layout(
      xaxis = list(title = list(text = "Posterior Intercepts",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = "Posterior Slopes",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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
    plotly::layout(
      xaxis = list(title = list(text = "Posterior Intercepts",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = "Posterior Slopes",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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
    plotly::layout(annotations = list(
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