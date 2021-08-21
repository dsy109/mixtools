library(mixtools)
library(plotly)
library(scales)

plotly.spEMN01 <- function(x, bw=x$bandwidth, knownpdf=dnorm, add.plot=FALSE,
                           width = 3 , col.dens = NULL, col.hist =  '#1f77b4',
                           title = "" , title.size = 15 , title.x = 0.5 , title.y = 0.95,
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
                   titlefont = list(size = ylab.size))
    )
  if (add.plot){
    plot <- plot%>%
      add_trace(x=x$data ,
                type = 'histogram', histnorm = "probability density",
                name = 'Data' , showlegend = FALSE,
                marker = list(color = col.hist,
                              line = list(color = col.hist))
      )%>%
      layout(bargap = 0.01)
  }
  print(plot)
}

## Probit transform of p-values
## from a Beta-Uniform mixture model
## comparion of parametric and semiparametric EM fit
## Note: in actual situations n=thousands
set.seed(50)
n=300 # nb of multiple tests
m=2 # 2 mixture components
a=c(1,0.1); b=c(1,1); lambda=c(0.6,0.4) # parameters
z=sample(1:m, n, rep=TRUE, prob = lambda)
p <- rbeta(n, shape1 = a[z], shape2 = b[z]) # p-values
o <- order(p)
cpd <- cbind(z,p)[o,] # sorted complete data, z=1 if H0, 2 if H1
p <- cpd[,2] # sorted p-values
y <- qnorm(p) # probit transform of the pvalues
# gaussian EM fit with component 1 constrained to N(0,1)
s1 <- normalmixEM(y, mu=c(0,-4),
                  mean.constr = c(0,NA), sd.constr = c(1,NA))
s2 <- spEMsymlocN01(y, mu0 = c(0,-3)) # spEM with N(0,1) fit
plotly.spEMN01(s2 , add.plot = FALSE)
