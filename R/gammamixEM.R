gammamixEM <- function(x, lambda = NULL, alpha = NULL, beta = NULL, k = 2, mom.start = TRUE, fix.alpha = FALSE,  
                            epsilon = 1e-08, maxit = 1000, maxrestarts = 20, verb = FALSE){
  if(!is.null(alpha) & is.numeric(fix.alpha))
    stop(paste("Both alpha and fix.alpha cannot be numeric!", 
               "\n"))
  if(!is.null(lambda)|!is.null(beta)) k <- max(length(lambda),length(beta))
  x <- as.vector(x)
  y <- x^(1/3)
  n <- length(x)
  if(is.logical(fix.alpha)&(fix.alpha==FALSE)){
    cond <- 1 
  } else if(is.logical(fix.alpha)&(fix.alpha==TRUE)){
    cond <- 2
  } else cond <- 3 
  fn.alpha <- function(alpha,beta,z,x) log(beta)+sum(z*log(x))/sum(z)-digamma(alpha)
  fn.alpha.2 <- function(alpha,beta,z,x) (log(beta)+sum(z*log(x))/sum(z)-digamma(alpha))^2
  fn.beta <- function(z,x,alpha) sum(z)/sum(z*x)*alpha
  
  fn.alpha2 <- function(alpha,beta,z,x) sum(z%*%cbind(log(beta)))/length(x)+sum(t(z)%*%cbind(log(x)))/length(x)-digamma(alpha)
  fn.alpha2.2 <- function(alpha,beta,z,x) (sum(z%*%cbind(log(beta)))/length(x)+sum(t(z)%*%cbind(log(x)))/length(x)-digamma(alpha))^2
  
  if(mom.start){
    out <- try(normalmixEM(y,k=k,maxit=5000,epsilon=1e-5,maxrestarts=10),silent=TRUE)
    if(inherits(out, "try-error", which = TRUE)){
      tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, 
                           beta = beta, k = k)
    } else{
      z <- out$posterior
      wt.mean <- sapply(1:k,function(i) weighted.mean(x,w=z[,i]))
      wt.var <- sapply(1:k,function(i) sum(z[,i] * (x - wt.mean[i])^2)/((n-1)*sum(z[,i])/n))
      shape.mom <- wt.mean^2/wt.var
      scale.mom <- wt.var/wt.mean
      shape.mle <- try(sapply(1:k, function(i) uniroot(fn.alpha,interval=c(0.000001,1000000),beta=1/scale.mom[i],z=z[,i],x=x)$root),silent=TRUE)    
      if(inherits(shape.mle, "try-error", which = TRUE)) shape.mle <- sapply(1:k, function(i) nlminb(shape.mom[i],fn.alpha.2,lower=0,beta=1/scale.mom[i],z=z[,i],x=x)$par)    
      scale.mle <- sapply(1:k, function(i) 1/fn.beta(z=z[,i],x=x,alpha=shape.mle[i]))
      lambda.mle <- apply(z,2,mean)
      tmp <- list(lambda=lambda.mle,alpha=scale.mle,beta=shape.mle)
    }
  } else tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, 
                              beta = beta, k = k)
  lambda.mle <- tmp$lambda
  scale.mle <- tmp$beta
  if(cond==2){
    shape.mle <- rep(mean(tmp$alpha),k)
  } else if(cond==1){
    shape.mle <- tmp$alpha
  } else shape.mle <- rep(fix.alpha,k)
  
  dens <- function(x, lambda, alpha, beta) {
    k <- length(lambda)
    temp = sapply(1:k, function(j) dgamma(x, shape = alpha[j], scale = beta[j]))
    temp = t(lambda * t(temp))
    temp
  }
  iter <- 0
  mr <- 0
  diff <- epsilon + 1
  dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
  old.obs.ll <- sum(log(apply(dens1, 1, sum)))
  ll <- old.obs.ll
  if(cond==1){
    
    while (diff > epsilon && iter < maxit) {
      old.shape.mle <- shape.mle
      old.scale.mle <- scale.mle
      old.lambda.mle <- lambda.mle
      z <- dens1/apply(dens1, 1, sum)
      #M-step
      shape.mle <- try(sapply(1:k, function(i) uniroot(fn.alpha,interval=c(0.000001,10000),beta=1/old.scale.mle[i],z=z[,i],x=x)$root),silent=TRUE)
      if(inherits(shape.mle, "try-error", which = TRUE)) shape.mle <- sapply(1:k, function(i) nlminb(old.shape.mle[i],fn.alpha.2,lower=0,beta=1/scale.mle[i],z=z[,i],x=x)$par)
      scale.mle <- sapply(1:k, function(i) 1/fn.beta(z=z[,i],x=x,alpha=shape.mle[i]))
      lambda.mle <- apply(z,2,mean)
      dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
      new.obs.ll <- sum(log(apply(dens1, 1, sum)))
      diff <- new.obs.ll - old.obs.ll
      old.obs.ll <- new.obs.ll
      ll <- c(ll, old.obs.ll)
      iter = iter + 1
      if (verb) {
        cat("iteration =", iter, " log-lik diff =", diff, 
            " log-lik =", new.obs.ll, "\n")
      }
      if(is.na(diff)){
        cat("Note: Choosing new starting values.", "\n")
        if (mr == maxrestarts) 
          stop(paste("Try different number of components?", 
                     "\n"))
        mr <- mr + 1
        tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, 
                             beta = beta, k = k)
        lambda.mle <- tmp$lambda
        if(cond==2){
          shape.mle <- rep(mean(tmp$alpha),k)
        } else if(cond==1){
          shape.mle <- tmp$alpha
        } else shape.mle <- rep(fix.alpha,k)
        scale.mle <- tmp$beta
        iter <- 0
        diff <- epsilon + 1
        dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
        old.obs.ll <- sum(log(apply(dens1, 1, sum)))
        ll <- old.obs.ll
      }
    }
  } else if(cond==2){
    
    while (diff > epsilon && iter < maxit) {
      old.shape.mle <- shape.mle
      old.scale.mle <- scale.mle
      old.lambda.mle <- lambda.mle
      z <- dens1/apply(dens1, 1, sum)
      #M-step
      shape.mle <- try(rep(uniroot(fn.alpha2,interval=c(0.000001,10000),beta=1/old.scale.mle,z=z,x=x)$root,k),silent=TRUE)
      if(inherits(shape.mle, "try-error", which = TRUE)) shape.mle <- rep(nlminb(old.shape.mle[1],fn.alpha2,lower=0,beta=1/old.scale.mle,z=z,x=x)$par,k)
      scale.mle <- sapply(1:k, function(i) 1/fn.beta(z=z[,i],x=x,alpha=shape.mle[1]))
      lambda.mle <- apply(z,2,mean)
      dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
      new.obs.ll <- sum(log(apply(dens1, 1, sum)))
      diff <- new.obs.ll - old.obs.ll
      old.obs.ll <- new.obs.ll
      ll <- c(ll, old.obs.ll)
      iter = iter + 1
      if (verb) {
        cat("iteration =", iter, " log-lik diff =", diff, 
            " log-lik =", new.obs.ll, "\n")
      }
      if(is.na(diff)){
        cat("Note: Choosing new starting values.", "\n")
        if (mr == maxrestarts) 
          stop(paste("Try different number of components?", 
                     "\n"))
        mr <- mr + 1
        tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, 
                             beta = beta, k = k)
        lambda.mle <- tmp$lambda
        if(cond==2){
          shape.mle <- rep(mean(tmp$alpha),k)
        } else if(cond==1){
          shape.mle <- tmp$alpha
        } else shape.mle <- rep(fix.alpha,k)
        scale.mle <- tmp$beta
        iter <- 0
        diff <- epsilon + 1
        dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
        old.obs.ll <- sum(log(apply(dens1, 1, sum)))
        ll <- old.obs.ll
      }
    }     
    
  } else{
    
    while (diff > epsilon && iter < maxit) {
      old.scale.mle <- scale.mle
      old.lambda.mle <- lambda.mle
      z <- dens1/apply(dens1, 1, sum)
      #M-step
      scale.mle <- sapply(1:k, function(i) 1/fn.beta(z=z[,i],x=x,alpha=shape.mle[1]))
      lambda.mle <- apply(z,2,mean)
      dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
      new.obs.ll <- sum(log(apply(dens1, 1, sum)))
      diff <- new.obs.ll - old.obs.ll
      old.obs.ll <- new.obs.ll
      ll <- c(ll, old.obs.ll)
      iter = iter + 1
      if (verb) {
        cat("iteration =", iter, " log-lik diff =", diff, 
            " log-lik =", new.obs.ll, "\n")
      }
      if(is.na(diff)){
        cat("Note: Choosing new starting values.", "\n")
        if (mr == maxrestarts) 
          stop(paste("Try different number of components or different value for the fixed shape parameter?", 
                     "\n"))
        mr <- mr + 1
        tmp <- gammamix.init(x = x, lambda = lambda, alpha = alpha, 
                             beta = beta, k = k)
        lambda.mle <- tmp$lambda
        if(cond==2){
          shape.mle <- rep(mean(tmp$alpha),k)
        } else if(cond==1){
          shape.mle <- tmp$alpha
        } else shape.mle <- rep(fix.alpha,k)
        scale.mle <- tmp$beta
        iter <- 0
        diff <- epsilon + 1
        dens1 <- dens(x=x, lambda=lambda.mle, alpha=shape.mle, beta=scale.mle)
        old.obs.ll <- sum(log(apply(dens1, 1, sum)))
        ll <- old.obs.ll
      }
    }     
  }
  if (iter == maxit) {
    cat("WARNING! NOT CONVERGENT!", "\n")
  }
  cat("number of iterations=", iter, "\n")
  theta = rbind(shape.mle, scale.mle)
  rownames(theta) = c("alpha", "beta")
  colnames(theta) = c(paste("comp", ".", 1:k, sep = ""))
  a = list(x = x, lambda = lambda.mle, gamma.pars = theta, loglik = new.obs.ll, 
           posterior = z, all.loglik = ll, ft = "gammamixEM")
  class(a) = "mixEM"
  a
}  
