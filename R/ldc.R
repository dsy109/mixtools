ldc <- function(data,class,score){
  data = as.matrix(data)
  K = length(unique(class)) # The number of classes (expect the classes to be labeled 1, 2, 3, ..., K-1, K 
  N = dim(data)[1] # The number of samples
  p = dim(data)[2] # The number of features
  
  # Compute the class dependent probabilities and class dependent centroids: 
  Pi = c(table(class))/N
  M = aggregate(data,list(class),mean)[,-1]
  # Compute W:
  W = cov(data) 
  # Compute M* = M W^{-1/2} using the eigen-decomposition of W :
  e = eigen(W)
  V = e$vectors
  W_Minus_One_Half = V %*% diag( 1/sqrt(e$values),nrow=p ) %*% t(V)
  MStar = as.matrix(M) %*% W_Minus_One_Half
  # Compute B* the covariance matrix of M* and its eigen-decomposition:
  if(p>1 & length(Pi)>1){
    BStar = cov(MStar)
    e = eigen(BStar)
    VStar = e$vectors
    # 1st linear discriminant coordinate
    ldc1 = W_Minus_One_Half %*% VStar[,score]
  } else if(p==1){
    ldc1 = c(W_Minus_One_Half)
  } else ldc1 = W_Minus_One_Half[,score]
  return(ldc1)
}
