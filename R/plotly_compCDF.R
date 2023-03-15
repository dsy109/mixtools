plotly_compCDF <- function(data, 
                           weights, 
                           x=seq(min(data, na.rm=TRUE), max(data, na.rm=TRUE), len=250), 
                           comp=1:NCOL(weights), 
                           makeplot=TRUE,
                           cex = 3,
                           width = 3,
                           legend.text = "Composition",
                           legend.text.size = 15,
                           legend.size = 15,
                           title = "Empirical CDF",
                           title.x = 0.5,
                           title.y = 0.95,
                           title.size = 15,
                           xlab = "Data",
                           xlab.size = 15,
                           xtick.size = 15,
                           ylab = "Probability",
                           ylab.size = 15,
                           ytick.size = 15,
                           col.comp = NULL){
  if (NROW(weights) != NROW(data)) {
    stop("data and weights arguments must have same number of rows")
  }
  if (is.null(col.comp)){
    col.comp <- hue_pal()(length(comp))
  }
  if (length(col.comp) != length(comp)){
    print(paste("Please specify",length(comp),"colors in 'col.comp'."))
  }
  # First, normalize the weights so the sum of each column is 1/NCOL(data)
  weights <- t(t(weights) / (NCOL(data) * colSums(weights)))
  # Next, give a binomial count for each row of the data and for each x
  f <- function(row, cutpt) colSums(outer(row, cutpt, "<="), na.rm = TRUE)
  bc <- apply(data, 1, f, x)
  # bc is a length(x) by n matrix; each column should be multiplied by
  # the appropriate weight(s) and then the rows summed to give the 
  # unnormalized cdf estimates. This is just a matrix product.
  cdfs <- bc %*% weights[,comp,drop=FALSE]
  
  if(makeplot) {
    plot <- plot_ly() %>%
      layout(
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
                     tickfont = list(size = ytick.size),
                     range = c(0 , 1)
        )
      )
    for (i in 1:length(comp)) {
      plot <- plot %>%
        add_trace(x=x , y=cdfs[,comp[i]] , 
                  type = 'scatter' , mode = 'lines+markers',
                  marker = list(size = cex , color = col.comp[i]),
                  line = list(width = width , color = col.comp[i]),
                  name = comp[i] , showlegend = TRUE) %>%
        layout(
          legend = list(title=list(text=legend.text,
                                   font=list(size=legend.text.size)),
                        font = list(size=legend.size))
        )
    }
    print(plot)
  }
  invisible(t(cdfs))
}
