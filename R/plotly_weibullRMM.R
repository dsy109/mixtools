plotly_weibullRMM <- function(a, title=NULL, rowstyle=TRUE, subtitle=NULL,
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
    plotly::layout(
      legend = list(font = list(size=legend.size)),
      title = list(text = paste(tt1 , "\n(", subtitle,")"),
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = xlab,
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = ylab,
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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
    plotly::layout(
      legend = list(font = list(size=legend.size)),
      title = list(text = paste(tt2 , "\n(", subtitle,")"),
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = xlab,
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = ylab,
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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
    plotly::layout(
      legend = list(font = list(size=legend.size)),
      title = list(text = paste(tt3 , "\n(", subtitle,")"),
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = xlab,
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = ylab,
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
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