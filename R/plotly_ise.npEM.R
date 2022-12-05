plotly_ise.npEM <- function(npEMout, component=1, block=1, truepdf=dnorm, lower=-Inf,
                            upper=Inf, plots = TRUE ,
                            col = NULL , width = 3,
                            title = NULL , title.size = 15 , title.x = 0.5 , title.y = 0.95,
                            xlab = "t" , xlab.size = 15 , xtick.size = 15,
                            ylab = "" , ylab.size = 15 , ytick.size = 15,
                            legend.text = "" , legend.text.size = 15 , legend.size = 15, ...){
  coords <- npEMout$blockid == block
  bs <- sum(coords) # block size
  xx <- as.vector(npEMout$data[,coords]) # flatten data 
  wts <- rep(npEMout$post[,component],bs) # duplicate weights
  if (is.matrix(npEMout$bandwidth)){
    bw <- npEMout$bandwidth[block,component]
  } else {bw <- npEMout$bandwidth}
  integrand = function(u,...) {
    (wkde(xx,u,wts,bw) - truepdf(u,...))^2
  }
  numint <- integrate(integrand,lower,upper, ...)
  
  if (is.null(col)){
    col <- hue_pal()(2)
  }
  if (length(col) != 2){
    print("Please specify 2 colors in 'col'.")
  }
  
  if (plots){
    # plot of estimated and truepdf
    ise <- paste(round(numint$value,4))
    temp=paste(component, block, sep="")
    if (is.null(title)){
      title = substitute(expression(paste("Integrated Squared Error for ",
                                          f[temp]," = ",ise,sep="")))
    }
    
    if (!is.finite(lower)) {
      lower <- min(xx)
    }
    if (!is.finite(upper)) {
      upper <- max(xx)
    }    
    u <- seq(lower,upper, 0.01)
    fhat <- wkde(xx,u,wts,bw)
    ymax <- max(max(truepdf(u, ...)),max(fhat))
    
    plot <- plot_ly()%>%
      add_trace(x=u , 
                y=fhat , type = 'scatter' , mode = 'lines',
                line = list(width = (width/2) , color = col[2]),
                name = "Fitted", showlegend = TRUE)%>%
      add_trace(x=u , 
                y=truepdf(u, ...) , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = col[1]),
                name = "True", showlegend = TRUE)%>%
      plotly::layout(
        legend = list(title=list(text=legend.text,
                                 font=list(size=legend.text.size)),
                      font = list(size=legend.size)),
        title = list(text = eval(title),
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
  }
  print(plot)
  numint
}