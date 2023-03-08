plotly_mixreg <-function (x,
                          xlab="X", xlab.size=15 , xtick.size=15,
                          ylab="Y", ylab.size=15 , ytick.size=15,
                          title="Estimated Mixed Regression", title.size=15,
                          title.x = 0.5,title.y=0.95,
                          col="#1f77b4", lwd=3, cex=6){
  mix.object <- x
  if (!inherits(mix.object, "mixEM"))
    stop("Use only with \"mixEM\" objects!")
  
  n.component <- dim(mix.object$beta)[2]
  line.x.min <- min(mix.object$x)-1
  line.x.max <- max(mix.object$x)+1
  
  if (length(col) < n.component){
    col <- hue_pal()(n.component)
  }
  
  plot.reg <- plot_ly()%>%
    add_trace(x=mix.object$x[,2] , 
              y=mix.object$y, 
              type = 'scatter' , mode = 'markers',
              marker = list(color = col[apply(mix.object$posterior,1, which.max)] , size = cex),
              name = "Data" , showlegend = FALSE)%>%
    layout(
      title = list(text = title,
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)
      ),
      xaxis = list(title = list(text = xlab,
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size),
                   zeroline = FALSE
      ),
      yaxis = list(title = list(text = ylab,
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size),
                   zeroline = FALSE
      )
    )
  
  for (i in 1:n.component) {
    b0 <- mix.object$beta[1,i]
    b1 <- mix.object$beta[2,i]
    line.y.min <- b0+b1*line.x.min
    line.y.max <- b0+b1*line.x.max
    
    plot.reg <- plot.reg %>%
      add_trace(x=c(line.x.min , line.x.max) ,
                y=c(line.y.min , line.y.max),
                type = 'scatter' , mode = 'lines', 
                line = list(width = lwd , color = col[i] , dash = "dash"),
                name = paste("Component" , i) , showlegend = FALSE)
  }
  
  print(plot.reg)
}

    
    
    
  