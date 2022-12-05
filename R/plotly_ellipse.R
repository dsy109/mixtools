plotly_ellipse <- function(mu, sigma, alpha=.05, npoints=250,
                           draw=TRUE,
                           cex = 3,
                           col = "#1f77b4",
                           lwd = 3,
                           title = "",
                           title.x = 0.5,
                           title.y = 0.95,
                           title.size = 15,
                           xlab = "X",
                           xlab.size = 15,
                           xtick.size = 15,
                           ylab = "Y",
                           ylab.size = 15,
                           ytick.size = 15
) {
  es <- eigen(sigma)
  e1 <- es$vec%*%diag(sqrt(es$val))
  r1 <- sqrt(qchisq(1-alpha,2))
  theta <- seq(0,2*pi,len=npoints)
  v1 <- cbind(r1*cos(theta),r1*sin(theta))
  pts=t(mu-(e1%*%t(v1)))
  
  updatemenus <- list(
    list(
      x = 1.2,
      y = 1,
      visible = TRUE,
      font = list(size=12),
      active = 0,
      type= 'dropdown',
      buttons = list(
        list(
          label = "Show Markers",
          method = "update",
          args = list(list(visible = c(TRUE, FALSE , FALSE)))),
        list(
          label = "Show Lines",
          method = "update",
          args = list(list(visible = c(FALSE , TRUE , FALSE)))),
        list(
          label = "Show Lines+Markers",
          method = "update",
          args = list(list(visible = c(FALSE , FALSE , TRUE)))))
    )
  )
  
  plot <- plot_ly() %>%
    add_markers(x=pts[,1] , 
                y=pts[,2] ,
                type = 'scatter' , mode = 'markers' ,
                marker = list(color = col , size = cex) , 
                name = 'Data' , showlegend = FALSE) %>%
    add_trace(x=pts[,1] , 
              y=pts[,2] , type = 'scatter' , mode = 'lines',
              line = list(color = col , width = lwd),
              name = "Data" , showlegend = FALSE) %>%
    add_trace(x=pts[,1] , 
              y=pts[,2] , type = 'scatter' , mode = 'lines+markers',
              marker = list(color = col , size = cex),
              line = list(color = col , width = lwd),
              name = "Data" , showlegend = FALSE) %>%
    plotly::layout(
      updatemenus = updatemenus,
      title = list(text = title,
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)
      ),
      xaxis = list(title = list(text = xlab,
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size)
      ),
      yaxis = list(title = list(text = ylab,
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size)
      )
    )
  print(plot)
  
  invisible(pts)
}
