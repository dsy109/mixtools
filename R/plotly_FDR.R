plotly_FDR <- function(post1, post2=NULL, lg1="FDR 1", lg2=NULL, 
                       compH0=1, alpha=0.1, complete.data =NULL, pctfdr=0.3,
                       col = NULL, width = 3 ,
                       title = NULL , title.size = 15 , title.x = 0.5 , title.y = 0.95,
                       xlab = "Index" , xlab.size = 15 , xtick.size = 15,
                       ylab = "Probability" , ylab.size = 15 , ytick.size = 15,
                       legend.text = "" , legend.text.size = 15 , legend.size = 15
                       ){
  hline <- function(y = 0, color = '#1f77b4') {
    list(
      type = "line",
      y0 = y,
      y1 = y,
      xref = "paper",
      x0 = 0,
      x1 = 1,
      line = list(color = '#1f77b4',
                  dash = "dash",
                  width = 1)
    )
  }
  
  if(is.null(col)){
    col <- hue_pal()(3)
  }
  if(length(col) != 3){
    print("Please specify 3 colors in 'col'.")
  }
  
  n <- dim(post1)[1]
  cs1 <- cumsum(post1[,compH0]) # local FDR(p_i)'s
  fdr1 <- cs1/(1:n) # FDR(p_i)'s
  if (is.null(title)) title <- paste("FDR estimate(s), n=",n)
  if (!is.null(post2)) {
    cs2 <- cumsum(post2[,compH0]) # local FDR(p_i)'s
    fdr2 <- cs2/(1:n)
    if (is.null(lg2)) {lg2 <- "FDR 2"}
  }
  i1 <- sum(fdr1<pctfdr) # plot up to pctfdr % FDR
  if (i1 == 0) {i1 <- n}   # for very bad fit, fdr[1] > pctfdr
  # cat("index",i1)
  plot <- plot_ly()%>%
    add_trace(x=seq(from = 1 , to = i1 , by = 1) , 
              y=fdr1[1:i1] , type = 'scatter' , mode = 'lines',
              line = list(width = width , color = col[1]),
              name = lg1, showlegend = TRUE)%>%
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
                   ),
      shapes = list(type = "line",
                    y0 = alpha,
                    y1 = alpha,
                    xref = "paper",
                    x0 = 0,
                    x1 = 1,
                    line = list(color = '#1f77b4',
                                dash = "dash",
                                width = 1)
      )
    )
    if (!is.null(post2)){
      plot <- plot%>%
        add_trace(x=seq(from = 1 , to = i1 , by = 1) , 
                  y=fdr2[1:i1] , type = 'scatter' , mode = 'lines',
                  line = list(width = width , color = col[2]),
                  name = lg2, showlegend = TRUE)
    }
  if (!is.null(complete.data)){
    V <- cumsum(complete.data[,1]==1) # cumulative nb of items under H0 
    trueFDR <- V/(1:n)
    plot <- plot%>%
      add_trace(x=seq(from = 1 , to = i1 , by = 1) , 
                y=trueFDR[1:i1] , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = col[3] , dash = "dash"),
                name = "True FDR", showlegend = TRUE)
  }
  print(plot)
}
