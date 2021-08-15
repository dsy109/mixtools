library(mixtools)
library(plotly)

plotly.compCDF <- function(data, 
                           weights, 
                           x=seq(min(data, na.rm=TRUE), max(data, na.rm=TRUE), len=250), 
                           comp=1:NCOL(weights), 
                           makeplot=TRUE,
                           cex = 3,
                           width = 3,
                           legend.text = "Empirical CDF",
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
                           ytick.size = 15){
  if (NROW(weights) != NROW(data)) {
    stop("data and weights arguments must have same number of rows")
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
    plot <- plot_ly()
    for (i in 1:length(comp)) {
      plot <- add_trace(plot,
                        plot,
                        x=x , y=cdfs[,comp[i]] , type = 'scatter' , mode = 'lines+markers',
                        marker = list(size = cex),
                        line = list(width = width),
                        name = comp[i] , showlegend = TRUE) %>%
        layout(
          legend = list(title=list(text=legend.text,
                                   font=list(size=legend.text.size)),
                        font = list(size=legend.size)),
          title = list(text = title,
                       x = title.x,
                       y = title.y,
                       font = list(size=title.size)),
          xaxis = list(title = xlab,
                       tickfont = list(size = xtick.size),
                       titlefont = list(size = xlab.size)),
          yaxis = list(title = ylab,
                       tickfont = list(size = ytick.size),
                       titlefont = list(size = ylab.size),
                       range = c(0 , 1))
        )
    }
    print(plot)
  }
  t(cdfs)
}

### Example ###
## The sulfur content of the coal seams in Texas
set.seed(100)
A <- c(1.51, 1.92, 1.08, 2.04, 2.14, 1.76, 1.17)
B <- c(1.69, 0.64, .9, 1.41, 1.01, .84, 1.28, 1.59)
C <- c(1.56, 1.22, 1.32, 1.39, 1.33, 1.54, 1.04, 2.25, 1.49)
D <- c(1.3, .75, 1.26, .69, .62, .9, 1.2, .32)
E <- c(.73, .8, .9, 1.24, .82, .72, .57, 1.18, .54, 1.3)
dis.coal <- makemultdata(A, B, C, D, E,
                         cuts = median(c(A, B, C, D, E)))
temp <- multmixEM(dis.coal)
## Now plot the components' CDF via the posterior probabilities
plotly.compCDF(dis.coal$x, temp$posterior, xlab="Sulfur")
