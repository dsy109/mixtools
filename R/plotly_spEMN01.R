plotly_spEMN01 <- function(x, bw=x$bandwidth, knownpdf=dnorm, add.plot=FALSE,
                           width = 3 , col.dens = NULL, col.hist =  '#1f77b4',
                           title = NULL , title.size = 15 , 
                           title.x = 0.5 , title.y = 0.95,
                           xlab = "t" , xlab.size = 15 , xtick.size = 15,
                           ylab = "Density" , ylab.size = 15 , ytick.size = 15,
                           legend.text = "Densities" , legend.text.size = 15 , legend.size = 15
){
  t <- seq(min(x$data), max(x$data), len=200)
  f1 <- x$lambdahat[1]*knownpdf(t)
  f2 <- x$lambdahat[2]*wkde(x$data-x$muhat, u=t-x$muhat, w=x$post[,2], bw=bw, sym=TRUE)
  f <- f1+f2
  
  if(is.null(col.dens)){
    col.dens <- hue_pal()(3)
  }
  if (length(col.dens) != 3){
    print("Please sepcify 3 colors in 'col.dens'.")
  }
  if (is.null(title)){
    title <- ""
  }
  plot <- plot_ly()%>%
    add_trace(x=t , 
              y=f , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col.dens[1]),
              name = "f", showlegend = TRUE)%>%
    add_trace(x=t , 
              y=f1 , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col.dens[2]),
              name = "f1", showlegend = TRUE)%>%
    add_trace(x=t , 
              y=f2 , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col.dens[3]),
              name = "f2", showlegend = TRUE)%>%
    plotly::layout(
      legend = list(title=list(text=legend.text,
                               font=list(size=legend.text.size)),
                    font = list(size=legend.size)),
      title = list(text = title,
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
  if (add.plot){
    plot <- plot%>%
      add_trace(x=x$data ,
                type = 'histogram', histnorm = "probability density",
                name = 'Data' , showlegend = FALSE,
                marker = list(color = col.hist,
                              line = list(color = col.hist))
      )%>%
      plotly::layout(bargap = 0.01)
  }
  print(plot)
}