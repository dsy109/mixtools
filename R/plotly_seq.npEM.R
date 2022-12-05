plotly_seq.npEM <- function(x, col = '#1f77b4' , width = 6,
                            xlab = "Iteration" , xlab.size = 15 , xtick.size = 15,
                            ylab.size = 15 , ytick.size = 15,
                            title.size = 15 , title.x = 0.5 , title.y = 0.95){
  r <- NCOL(x$data)
  n <- NROW(x$data)
  m <- length(x$lambdahat)
  iter <- NROW(x$lambda)
  nbcol <- 1
  if (!is.null(x$symmetric) && x$symmetric){nbcol <- 2}
  # in all cases, plots the lambda's  
  for (j in 1:m){
    estim <- paste(round(x$lambdahat[j],3))
    # tt <- substitute(expression(paste("sequence of ",lambda[j],
    #                                   ", estimate ",widehat(lambda[j]),"=", estim, sep="")))
    tt <- paste("Sequence of " ,'&#955;',"<sub>",j,"</sub>",
                " (Estimated ","&#955;","<sub>",j,"</sub>","=",estim, ")",sep="")
    # ylabel <- substitute(expression(paste(lambda[j],sep="")))
    ylabel <- paste('&#955;',"<sub>",j,"</sub>",sep="")
    plot1 <- plot_ly()%>%
      add_trace(x = seq(from = 1 , to = length(x$lambda[,j]) , by = 1),
                y = x$lambda[,j] , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = col),
                showlegend = FALSE)%>%
      add_trace(x = c(0 , iter),
                y = rep(x$lambdahat[j],2) , type = 'scatter' , mode = 'lines',
                line = list(width = width , color = 'red' , dash = "dash"),
                showlegend = FALSE)%>%
      plotly::layout(
        title = list(text = tt,
                     x = title.x,
                     y = title.y,
                     font = list(size=title.size)),
        xaxis = list(title = list(text = xlab,
                                  font = list(size = xlab.size)),
                     tickfont = list(size = xtick.size)
        ),
        yaxis = list(title = list(text = eval(ylabel),
                                  font = list(size = ylab.size)),
                     tickfont = list(size = ytick.size)
        )
      )
    print(plot1)
  }
  ## for symmetric location spEM case plots mu
  if (!is.null(x$symmetric) && x$symmetric){
    for (j in 1:m){
      estim <- paste(round(x$muhat[j],3))
      # tt <- substitute(expression(paste("sequence of ",mu[j],
      #                                   ", estimate ",widehat(mu[j]),"=",estim,sep="")))
      tt <- paste("Sequence of " ,'&#956;',"<sub>",j,"</sub>",
                  " (Estimated ","&#956;","<sub>",j,"</sub>","=",estim, ")",sep="")
      # ylabel <- substitute(expression(paste(mu[j],sep="")))
      ylabel <- paste('&#956;',"<sub>",j,"</sub>",sep="")
      plot2 <- plot_ly()%>%
        add_trace(x = seq(from = 1 , to = length(x$mu[,j]) , by = 1),
                  y = x$mu[,j] , type = 'scatter' , mode = 'lines',
                  line = list(width = width , color = col),
                  showlegend = FALSE)%>%
        add_trace(x = c(0 , iter),
                  y = rep(x$muhat[j],2) , type = 'scatter' , mode = 'lines',
                  line = list(width = width , color = 'red' , dash = "dash"),
                  showlegend = FALSE)%>%
        plotly::layout(
          title = list(text = tt,
                       x = title.x,
                       y = title.y,
                       font = list(size=title.size)),
          xaxis = list(title = list(text = xlab,
                                    font = list(size = xlab.size)),
                       tickfont = list(size = xtick.size)
          ),
          yaxis = list(title = list(text = eval(ylabel),
                                    font = list(size = ylab.size)),
                       tickfont = list(size = ytick.size)
          )
        )
      print(plot2)
    }
  }
}