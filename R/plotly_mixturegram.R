plotly_mixturegram <- function(
  data, pmbs, method=c("pca","kpca","lda"), 
  all.n=FALSE, id.con=NULL, score=1, iter.max=50, nstart=25,
  xlab = "K" , xlab.size = 15 , xtick.size = 15,
  ylab = NULL , ylab.size = 15 , ytick.size = 15,
  cex = 12 , col.dot = "red" , width = 1 ,
  title = "Mixturegram" , title.size = 15 , title.x = 0.5 , title.y = 0.95
){
  vline <- function(x = 0, color = '#1f77b4') {
    list(
      type = "line",
      y0 = 0,
      y1 = 1,
      yref = "paper",
      x0 = x,
      x1 = x,
      line = list(color = '#1f77b4',
                  dash = "dash",
                  width = 1)
    )
  }
  
  col.blind=rep(c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#7FFF00","#7D26CD"),100)
  
  method <- match.arg(method)
  
  if(is.null(ylab)){
    if (method == "pca"){ylab <- "PC Scores"}
    else if (method == "kpca"){ylab <- "Kernel PC Scores"}
    else if (method == "lda"){ylab <- "LDC Scores"}
  }
  
  k=length(pmbs)+1
  if (is.null(id.con)) {id.con=lapply(pmbs,function(x) data%*%x/apply(x,2,sum))}
  ### Find ylim: lim ###
  ### Regular PCA ###
  if(method=="pca"){
    x.star=c(list(data), lapply(1:(k-1), function(i) cbind(data,pmbs[[i]][,order(id.con[[i]])])))
    Data=lapply(x.star,scale)
    if(score>(ncol(Data[[1]])+1)) {
      warning(paste("The largest value that can be specified for 'score' is ",
                    ncol(Data[[1]])+1,
                    ", which is what will be used for the mixturegram.",sep=""))
    }
    score<-min(score,ncol(x.star[[2]]))
    if(score==1){
      PCA=lapply(Data,function(x) x%*%princomp(x)$loadings[,score])
    } else{
      PCA <- vector("list",k)
      PCA[[1]]=Data[[1]]%*%princomp(Data[[1]])$loadings[,1]
      PCA[2:k]=lapply(1:(k-1),function(i) Data[[(i+1)]]%*%princomp(Data[[(i+1)]])$loadings[,score])
    }
    K=lapply(2:k,function(i) kmeans(PCA[[i]],i,iter.max=iter.max,nstart=nstart))
    lim=ceiling(max(abs(unlist(PCA))))
    mix.method <- PCA
  }else{
    ### Kernel PCA ###
    if(method=="kpca"){
      x.star=c(list(data), lapply(1:(k-1), function(i) cbind(data,pmbs[[i]][,order(id.con[[i]])])))
      Data=lapply(x.star,scale)
      if(score>(ncol(Data[[1]])+1)) warning(paste("The largest value that can be specified for 'score' is ",ncol(Data[[1]])+1,", which is what will be used for the mixturegram.",sep=""))
      score<-min(score,ncol(x.star[[2]]))
      if(score==1){
        kPCA=lapply(Data,function(x) cbind(pcv(kpca(x))[,1]))
      } else{
        kPCA <- vector("list",k)
        kPCA[[1]]=cbind(pcv(kpca(Data[[1]]))[,1])
        kPCA[2:k]=lapply(1:(k-1),function(i) cbind(pcv(kpca(Data[[i+1]]))[,score]))
      }
      K=lapply(2:k,function(i) kmeans(kPCA[[i]],i,iter.max=iter.max,nstart=nstart))
      lim=max(abs(unlist(kPCA)))+0.1
      mix.method <- kPCA
    } else if(method=="lda"){
      class=lapply(pmbs, function(post) apply(post,1,which.max))
      ldcdata = c(list(as.matrix(data) %*% ldc(data, rep(1,nrow(as.matrix(data))),score=score)),
                  lapply(class, function(class) as.matrix(data) %*% ldc(data, class,score=score)))
      K=lapply(2:k,function(i) kmeans(ldcdata[[i]],i,iter.max=iter.max,nstart=nstart))
      lim=ceiling(max(abs(unlist(ldcdata))))
      mix.method <- ldcdata
    }
  }
  plot <- plot_ly()%>%
    plotly::layout(
      title = list(text = title,
                   x = title.x,
                   y = title.y,
                   font = list(size=title.size)),
      xaxis = list(title = list(text = xlab,
                                font = list(size = xlab.size)),
                   tickfont = list(size = xtick.size),
                   range = c(0.25,k+0.25),
                   dtick = 1, 
                   tick0 = 1, 
                   tickmode = "linear"
      ),
      yaxis = list(title = list(text = ylab,
                                font = list(size = ylab.size)),
                   tickfont = list(size = ytick.size),
                   range = c(-lim, lim)
      ),
      shapes = lapply(1:k , function(i){vline(i)})
    )
  
  Kcol=function(x){
    temp=unlist(sapply(1:length(x$size),function(i) rep(rank(x$center)[i],x$size[i])))
    index=unlist(sapply(1:length(x$size),function(i) which(x$cluster==i)))
    K.col=replace(x$cluster,index,temp)
  }
  all.K.col=lapply(1:(k-1), function(i) Kcol(K[[i]]))
  K.col=all.K.col[[k-1]]
  n=length(K.col)
  K.centers=c(0, lapply(1:(k-1), function(i) sort(c(K[[i]]$centers))))
  
  if (all.n){
    for (i in 1:(k-1)){
      for (j in 1:length(mix.method[[i]])){
        plot <- plot%>%
          add_trace(x=c(i , i+1) , 
                    y=c(mix.method[[i]][j] , mix.method[[i+1]][j]) , type = 'scatter' , mode = 'lines',
                    line = list(width = width , color = col.blind[K.col][j]), 
                    name = paste("Method" , i),showlegend = FALSE) 
      }
    }
    for (q in 2:k){
      plot <- plot%>%
        add_trace(x=rep(q,q), 
                  y=sort(c(K[[q-1]]$centers)) , type = 'scatter' , mode = 'markers',
                  marker = list(size = cex , color = col.dot),
                  name = paste("Cluster" , q), showlegend = FALSE)
    }
    plot <- plot%>%
      add_trace(x=1, 
                y=mean(mix.method[[1]]) , type = 'scatter' , mode = 'markers',
                marker = list(size = cex , color = col.dot),
                name = "Cluster 1",showlegend = FALSE)
  } else{
    ride.medians=sapply(1:k, function(i) unlist(by(c(mix.method[[i]]),col.blind[K.col],median)))
    prop.mad=sapply(1:k, function(i) unlist(by(c(mix.method[[i]]),col.blind[K.col],mad)))
    L1=ride.medians-prop.mad
    U1=ride.medians+prop.mad
    L2=ride.medians-2*prop.mad
    U2=ride.medians+2*prop.mad
    L3=ride.medians-3*prop.mad
    U3=ride.medians+3*prop.mad
    
    srt.colors=rownames(ride.medians)
    srt.colors1=adjustcolor(srt.colors, alpha.f = 0.7)
    srt.colors2=adjustcolor(srt.colors, alpha.f = 0.4)
    srt.colors3=adjustcolor(srt.colors, alpha.f = 0.2)
    
    for (l in 1:k){
      plot <- plot%>%
        add_trace(
          x = c(1:k,k:1),
          y = c(L3[l,],rev(U3[l,])),
          type = 'scatter', mode = "lines",
          fill = 'tozeroy', ## fill type
          fillcolor = srt.colors3[l],
          hoveron = 'points+fills',
          line = list(
            color = srt.colors3[l]
          ),
          text = "Points + Fills",
          hoverinfo = 'text',
          showlegend = FALSE
        )%>%
        add_trace(
          x = c(1:k,k:1),
          y = c(L2[l,],rev(U2[l,])),
          type = 'scatter', mode = "lines",
          fill = 'tozeroy', ## fill type
          fillcolor = srt.colors2[l],
          hoveron = 'points+fills',
          line = list(
            color = srt.colors2[l]
          ),
          text = "Points + Fills",
          hoverinfo = 'text',
          showlegend = FALSE
        )%>%
        add_trace(
          x = c(1:k,k:1),
          y = c(L1[l,],rev(U1[l,])),
          type = 'scatter', mode = "lines",
          fill = 'tozeroy', ## fill type
          fillcolor = srt.colors1[l],
          hoveron = 'points+fills',
          line = list(
            color = srt.colors1[l]
          ),
          text = "Points + Fills",
          hoverinfo = 'text',
          showlegend = FALSE
        )
    }
    # for (h in 1:(k-1)){
    #   plot <- plot%>%
    #     add_trace(x=c(rep(h,k),rep((h+1),k)) , 
    #               y=c(ride.medians[,h] , ride.medians[,(h+1)]) , type = 'scatter' , mode = 'lines',
    #               line = list(width = width , color = srt.colors), 
    #               showlegend = FALSE) 
    # }
    for (j in 2:k){
      plot <- plot%>%
        add_trace(x=rep(j,j), 
                  y=sort(c(K[[j-1]]$centers)) , type = 'scatter' , mode = 'markers',
                  marker = list(size = cex , color = col.dot),
                  name = paste("Cluster" , j),showlegend = FALSE)
    }
    plot <- plot%>%
      add_trace(x=1, 
                y=mean(mix.method[[1]]) , type = 'scatter' , mode = 'markers',
                marker = list(size = cex , color = col.dot),
                name = "Cluster 1",showlegend = FALSE)
  }
  print(plot)
  props=c(1,sapply(1:length(K), function(i) K[[i]][[5]]/sum(unlist(K[[i]][5:6]))))
  print(list(stopping=props))
}
