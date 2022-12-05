plotly_spRMM <- function(sem, tmax = NULL,
                        width = 3 , col = '#1f77b4', cex = 3,
                        title.size = 15 , title.x = 0.5 , title.y = 0.95,
                        xlab.size = 15 , xtick.size=15 ,
                        ylab.size = 15 , ytick.size=15){
  t <- sem$t   
  ym <- max(sem$all.scaling)
  
  plot1 <- plot_ly()%>%
    add_trace(x = seq(from = 1 , to = length(sem$all.scaling) , by = 1),
              y = sem$all.scaling , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col),
              showlegend = FALSE)%>%
    plotly::layout(
      title = list(text = "Scaling",
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = "Iteration",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = "",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size),
                   range = c(0,ym)
      )
    )
  
  plot2 <- plot_ly()%>%
    add_trace(x = seq(from = 1 , to = length(sem$all.lambda[,1]) , by = 1),
              y = sem$all.lambda[,1] , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col),
              showlegend = FALSE)%>%
    plotly::layout(
      title = list(text = "Weight of Component 1",
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = "Iteration",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = "",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size),
                   range = c(0,1)
      )
    )
  
  if (is.null(tmax)){tmax <- max(sem$scaling*t) + 2}
  u <- seq(0, tmax, len=200)
  fhat <- sem$s.hat(sem$t.hat)*sem$hazard.hat
  ffinal <- sem$survival(sem$final.t)*sem$hazard
  
  plot3 <- plot_ly()%>%
    add_trace(x = seq(from = 1 , to = length(sem$survival) , by = 1),
              y = sem$survival , type = 'scatter' , mode = 'markers',
              marker = list(size = cex , color = col),
              showlegend = FALSE)%>%
    plotly::layout(
      title = list(text = "Survival Function Estimate",
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = "Time",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size),
                   range = c(0,tmax)
      ),
      yaxis = list(title = list(text = "",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size),
                   range = c(0,1)
      )
    )
  
  plot4 <- plot_ly()%>%
    add_trace(x = sem$final.t,
              y = ffinal , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col),
              showlegend = FALSE)%>%
    plotly::layout(
      title = list(text = "Density Estimate",
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = "Time",
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size),
                   range = c(0,tmax)
      ),
      yaxis = list(title = list(text = "Density",
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
    )
  
  print(plot1)
  print(plot2)
  print(plot3)
  print(plot4)
}