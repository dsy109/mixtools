mixturegram <- function(data, pmbs, method=c("pca","kpca","lda"), all.n=FALSE, id.con=NULL, score=1, iter.max=50, nstart=25, ...){
  col.blind=rep(c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
                  "#7FFF00","#7D26CD"),100)
  method <- match.arg(method)
  k=length(pmbs)+1
  if (is.null(id.con)) id.con=lapply(pmbs,function(x) data%*%x/apply(x,2,sum))
  ###regular pca
  if(method=="pca"){
    x.star=c(list(data), lapply(1:(k-1), function(i) cbind(data,pmbs[[i]][,order(id.con[[i]])])))
    Data=lapply(x.star,scale)
    if(score>(ncol(Data[[1]])+1)) warning(paste("The largest value that can be specified for 'score' is ",ncol(Data[[1]])+1,", which is what will be used for the mixturegram.",sep=""))
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
    plot(1:k,rep(1e1000,k),ylim=c(-lim,lim),col='white',xlab='k',ylab='PC Scores',xaxt="n", ...)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[356])
    axis(1,at=1:k,labels=1:k)  
    abline(v=1:k,col="gray80")
    mix.method <- PCA
  }else{
    ###kernel pca
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
      plot(1:k,rep(1e1000,k),ylim=c(-lim,lim),col='white',xlab='k',ylab='Kernel PC Scores',xaxt="n", ...)
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[356])
      axis(1,at=1:k,labels=1:k)  
      abline(v=1:k,col="gray80")
      mix.method <- kPCA
    }else{
      if(method=="lda"){
        class=lapply(pmbs, function(post) apply(post,1,which.max))
        ldcdata = c(list(as.matrix(data) %*% ldc(data, rep(1,nrow(as.matrix(data))),score=score)),
                    lapply(class, function(class) as.matrix(data) %*% ldc(data, class,score=score)))
        K=lapply(2:k,function(i) kmeans(ldcdata[[i]],i,iter.max=iter.max,nstart=nstart))
        lim=ceiling(max(abs(unlist(ldcdata))))
        plot(1:k,rep(1e1000,k),ylim=c(-lim,lim),col='white',xlab='k',ylab='LDC Scores',xaxt="n", ...)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = colors()[356])
        axis(1,at=1:k,labels=1:k)  
        abline(v=1:k,col="gray80")
        mix.method <- ldcdata
      }
    }
  }
  Kcol=function(x){
    temp=unlist(sapply(1:length(x$size),function(i) rep(rank(x$center)[i],x$size[i])))
    index=unlist(sapply(1:length(x$size),function(i) which(x$cluster==i)))
    K.col=replace(x$cluster,index,temp)
  }
  all.K.col=lapply(1:(k-1), function(i) Kcol(K[[i]]))
  K.col=all.K.col[[k-1]]
  n=length(K.col)
  K.centers=c(0, lapply(1:(k-1), function(i) sort(c(K[[i]]$centers))))
  if(all.n){
    sapply(1:(k-1),function(i) segments(i,mix.method[[i]],(i+1),mix.method[[i+1]],col=col.blind[K.col]))
    points(1,mean(mix.method[[1]]),pch=21,cex=1.2,bg=colors()[553],col=1)
    sapply(2:k,function(i) points(rep(i,i),sort(c(K[[i-1]]$centers)),pch=21,cex=1.2,bg=colors()[553],col=1))
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
    
    
    invisible(sapply(1:k, function(i) polygon(c(1:k,k:1),c(L3[i,],rev(U3[i,])),col=srt.colors3[i], border=FALSE) ))
    invisible(sapply(1:k, function(i) polygon(c(1:k,k:1),c(L2[i,],rev(U2[i,])),col=srt.colors2[i], border=FALSE) ))
    invisible(sapply(1:k, function(i) polygon(c(1:k,k:1),c(L1[i,],rev(U1[i,])),col=srt.colors1[i], border=FALSE) ))
    invisible(sapply(1:(k-1),function(i) segments(rep(i,k),ride.medians[,i],rep((i+1),k),ride.medians[,(i+1)],col=srt.colors) ))
    points(1,mean(mix.method[[1]]),pch=21,cex=1.2,bg=colors()[553],col=1)
    invisible(sapply(2:k,function(i) points(rep(i,i),sort(c(K[[i-1]]$centers)),pch=21,cex=1.2,bg=colors()[553],col=1)))
  }
  props=c(1,sapply(1:length(K), function(i) K[[i]][[5]]/sum(unlist(K[[i]][5:6]))))
  print(list(stopping=props))
}
