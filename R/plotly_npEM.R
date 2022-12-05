plotly_spEM <- plotly_npEM <- function(
  x, blocks = NULL, hist=TRUE, addlegend=TRUE,
  scale = TRUE, title=NULL, breaks="Sturges", 
  dens.col = NULL, newplot=TRUE, ylim = NULL ,
  col.hist = "#1f77b4",
  width = 3, title.x = 0.5 , title.y = 0.95, title.size = 15,
  xlab = "X" , xlab.size = 15 , xtick.size = 15,
  ylab = "Density" , ylab.size = 15 , ytick.size = 15,
  legend.text = "Posteriors",
  legend.text.size = 15,
  legend.size = 15
  ){
  r <- NCOL(x$data)
  m <- NCOL(x$posteriors)
  if(is.null(dens.col)){
    dens.col <- hue_pal()(m)
  }
  if (length(dens.col) > m){
    print(paste("Please specify",m,"colors in 'dens.col'."))
  }
  blockid <- x$blockid
  if (is.null(blocks)) {
    if(!is.null(blockid)) {
      blocks <- 1:max(blockid)
    } else {
      blocks <- blockid <- 1:r
    }
  }
  ylim.orig <- ylim
  out <- list(x=list(), y=list())
  if (!newplot) {
    hist <- FALSE
  }
  ############################
  plot.all <- plot_ly()%>%
    plotly::layout(
      legend = list(title=list(text=legend.text,
                               font=list(size=legend.text.size)),
                    font = list(size=legend.size)),
      title = list(text = "Densities for Different Posteriors",
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
  for(i in 1:length(blocks)){
    coords <- blockid == blocks[i]
    ylim <- ylim.orig
    if (is.null(title)) {
      if (r>1) {
        tt <- paste(which(coords), collapse=",")
        tt <- paste("Coordinate", ifelse(sum(coords)>1, "s ", " "), tt, sep="")
      } else {
        tt <- "Density Curves"
      }
    } else {
      tt <- rep(title,length(blocks))[i]
    }
    dx <- dy <- NULL
    for (j in 1:m) {
      d <- density(x, component=j, block=blocks[i], scale=scale)
      dx <- cbind(dx, d$x)
      dy <- cbind(dy, d$y)
    }
    xx <- as.vector(as.matrix(x$data)[,coords])
    if (is.null(ylim)) {
      ylim=range(dy)
      if (hist) {
        ylim[2] <- max(ylim[2], hist(xx, breaks=breaks, plot=FALSE)$density)
      }
    }
    if (newplot){
      plot.new <- plot_ly()%>%
        plotly::layout(
          legend = list(title=list(text=legend.text,
                                   font=list(size=legend.text.size)),
                        font = list(size=legend.size)),
          title = list(text = tt,
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
      for (j in 1:m){
        plot.new <- plot.new %>%
          add_trace(x=dx[,j] , 
                    y=dy[,j] , type = 'scatter' , mode = 'lines',
                    line = list(width = width , color = dens.col[j]),
                    name = paste("Posterior",j), showlegend = addlegend)
      }
      if (hist){
        plot.new <- plot.new%>%
          add_trace(x=xx ,
                    type = 'histogram', histnorm = "probability density",
                    name = 'Data' , showlegend = FALSE,
                    marker = list(color = col.hist,
                                  line = list(color = col.hist))
          )%>%
          plotly::layout(bargap = 0.01)
      }
      print(plot.new)
    } else {
      if (i > 1){addlegend <- FALSE}
      for (j in 1:m){
        plot.all <- plot.all %>%
          add_trace(x=dx[,j] , 
                    y=dy[,j] , type = 'scatter' , mode = 'lines',
                    line = list(width = width , color = dens.col[j]),
                    name = paste("Posterior",j), showlegend = addlegend)
      }
    }
  }
  if (!newplot){print(plot.all)}
}