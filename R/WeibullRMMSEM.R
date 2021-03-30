# code for Reliability Mixture Models (RMM) with Censored data
# D. Chauveau
# ref: Bordes L. and Chauveau D. Computational Statistics (2016)

# Simulate from a Weibull mixture
# lambda = vector of component probabilities
#	shape, scale = vector of component rates
rweibullmix <- function(n,lambda=1, shape=1, scale=1) {
	m <- length(lambda) # nb of components
	z <- sample(m, n, replace=TRUE, prob=lambda) # component indicator
	rweibull(n, shape = shape[z], scale=scale[z])
	}


########################################################################
## Stochastic EM algorithm for Reliability Mixture Models (RMM) 
## with Censoring; Parametric model, for 
## univariate finite mixture of Weibull with right censoring
# x = lifetime data, censored by random c if d is not NULL, in which case
# x= min(x,c) and d = I(x <= c)
# uses parametric MLE for censored weibull data from the survival package
# caution: when fitted by a survreg object, the weibull parametrization is
#			shape=1/fit$scale and scale=exp(fit$coeff)
# averaged = TRUE if averaging of the parameters is computed at each iteration
#			 (cf Nielsen 2000)
# 	NB: averaging can be done in several ways; 
#	DEPRECATED: the current theta for E & M steps is the average over the sequence,
#	but the theta^next not averaged is stored on the sequence 
# 	CURRENT VERSION: the average is stored and used for next step
#
# About calls of survreg() from the survival package:
# 	maxit.survreg is
#	passed to survreg.control() for weibull MLE using survreg()

weibullRMM_SEM <- function (x, d=NULL, lambda = NULL, 
			shape = NULL, scale = NULL,  # dweibull() parameters
			k = 2, # default nb of components
          	maxit = 200, # maxrestarts = 20, 
          	maxit.survreg = 200, epsilon = 1e-03,
          	averaged = TRUE, verb = FALSE) {
#  warn <- options(warn=-1) # Turn off warnings
#  x <- as.vector(x)
#  require(survival)
	n <- length(x)
	if (!is.null(lambda)) k=length(lambda)
	if (!is.null(scale)) k=length(scale) # at least one
	if (!is.null(shape)) k=length(shape) # should be define !!!

	if (is.null(d))  d <- rep(1,n)  # noncensored case forced
	xx <- matrix(x, nrow=n, ncol=k) # x repeated k times, used in E-steps
	dd <- matrix(d, nrow=n, ncol=k) # idem for d
#
# note: d can be used instead of dd since automatic recycling gives
# dd*xx = d*xx, but we keep dd to avoid eventual futur change
	
# init function call to do later, in case lambda = rate = NULL  
	if (is.null(lambda)) lambda <- rep(1/k,k)
	if (is.null(scale)) scale <- rep(1,k)  # to do : init.functions(k)
	if (is.null(shape)) shape <- rep(1,k)  # 
# sequences for storing along iterations 	
	lambda_seq <- scale_seq <- shape_seq <- matrix(0, nrow = maxit, ncol = k)
 	lambda_seq[1,] <- lambda
 	scale_seq[1,] <- scale;  shape_seq[1,] <- shape
	loglik_seq <- NULL
	oldloglik <- -Inf
#    notdone <- TRUE # for handling restarts etc, ToDo later!
#    while(notdone) {
      # Initialize everything
      notdone <- FALSE
      # dll <- epsilon+1
      iter <- 1
      post <- z <- sumpost <- matrix(0, nrow = n, ncol = k)
      new_scale <- new_shape <- rep(0,k)
	  
	  ##### SEM iterations #####
      while (iter < maxit) { # SEM version
      	
		### E-step ###
		scn <- matrix(scale, n, k, byrow=T) # for vectorized post comput.
		shn <- matrix(shape, n, k, byrow=T)
		ll <- matrix(lambda, n, k, byrow=T)
		# handling censored & non-censored cases
		post <- ((ll*dweibull(xx,shn,scn))^dd)*((ll*(1-pweibull(xx,shn,scn)))^(1-dd))
		rs <- rowSums(post)
		loglik <- sum(log(rs)) # loglik without the constant term related to h pdf
		# post normalized per row
		post <- sweep(post, 1, rs, "/") # posteriors p_{ij}^t 's
		# check and solve NaN's: may cause theoretical pbs!?
		snans <- sum(is.na(post))
		if (snans > 0) {
			post[is.na(post[,1]),] <- 1/k
			cat("warning:",snans, "NaN's in post\n")
			}

		### S-step ###	
		# ~ matrix of component indicators simu checked ?
		z <- t(apply(post, 1, function(prob) rmultinom(1, 1, prob))) 
		nsets <- colSums(z) # subsets per component sizes
		
		# cat("it",iter,": sets=",nsets,"\n")
		
		### M-step ###
		# new_parameter = SEM(lambda,shape,scale) 
		new_lambda <- nsets/n # or colMeans(post) if EM version preferred
		for (j in 1:k) { # for each component; vectorize later?
			tj <- x[z[,j]==1] # subsample from component j
			dj <- d[z[,j]==1] # associated event indicator
			# passing maxit and epsilon parameters to survreg
			# and current shape & scale as init parameters
			fit=survreg(Surv(tj,dj)~1, dist = 'weibull',
						control = survreg.control(maxiter = maxit.survreg, 
											rel.tolerance=epsilon),
						init.beta=log(scale), init.scale=1/shape)
			new_scale[j] <- exp(fit$coeff)
			new_shape[j] <- 1/fit$scale
			}
		
		# Next parameter value, depending on average strategy
		if (averaged) {
			scale <- (new_scale + iter*scale_seq[iter,])/(iter+1)
			shape <- (new_shape + iter*shape_seq[iter,])/(iter+1)
			lambda <- (new_lambda + iter*lambda_seq[iter,])/(iter+1)
			
			} else {	# averaged=FALSE case, just use last update 
				scale <- new_scale
				shape <- new_shape
				lambda <- new_lambda
				} 		
		# new strategy = storing sequence of averages
		lambda_seq[iter+1, ] <- lambda
		scale_seq[iter+1, ] <- scale
		shape_seq[iter+1, ] <- shape

        dll <- loglik - oldloglik # = Inf for iter=0 1st time
        oldloglik <- loglik
        loglik_seq <- c(loglik_seq, loglik)
        if (verb) {
          cat("iteration", iter, " log-lik diff =", dll, " log-lik =", 
              loglik, "\n")
#          print(rbind(lambda, rate))		 
        }
        iter <- iter + 1
      }  # end of SEM loops over iterations
        
#    } # while notdone, if restarts control implemented
        
    # final estimates depending on average strategy
    if (averaged) {
    	final.lambda <- lambda
    	final.scale <- scale
    	final.shape <- shape
    	} else {
    		final.scale <- colMeans(scale_seq[1:iter,])
			final.shape <- colMeans(shape_seq[1:iter,])
			final.lambda <- colMeans(lambda_seq[1:iter,])
    		}
    
    cat("number of iterations =", iter, "\n")
	colnames(post) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, d=d, 
    		lambda = final.lambda, 
    		scale = final.scale, shape = final.shape,
    		loglik = loglik, 
            posterior = post, all.loglik=loglik_seq, 
            all.lambda = lambda_seq[1:iter,],
            all.scale = scale_seq[1:iter,],
            all.shape = shape_seq[1:iter,],
            ft="weibullRMM_SEM")
	class(a) = "mixEM"
  a
}



##################################################
# plot SEM sequences from weibullRMM_SEM:
# color by components, one plot per parameter type
plotweibullRMM <- function(a, title=NULL, rowstyle=TRUE, subtitle=NULL,...)
	{
	n <- length(a$x); m <- dim(a$all.lambda)[2]
	pcc <- round(100*(1-mean(a$d)),2)
	if (is.null(subtitle)) {
				subtitle <- paste("n=",n,", ", pcc, "% censored", sep="")}
	if (is.null(title)) {
		tt1 <- "Shape parameters"; tt2 <- "Scale parameters" 
		tt3 <- "Weight parameters"
		} else tt1 <- tt2 <- tt3 <- title
lwdset <- 2
if (rowstyle) par(mfrow=c(3,1)) else par(mfrow=c(1,3))

plot(a$all.shape[,1], type="l", ylim=c(0,max(a$all.shape)), 
	xlab="iterations", ylab="estimates", main=tt1, ...)
# if (truevalues) abline(shape[1],0,lty=3)
title(sub=subtitle, cex.sub = 0.75)
lgd <- expression(sh[1]); lcol <- 1
for (j in 2:m) {
	lines(a$all.shape[,j], col=j, ...)
	# if (truevalues) abline(shape[j],0,col=j,lty=3)
	lgd <- c(lgd,substitute(sh[j])); lcol <- c(lcol,j)
	}
legend("topright", lgd, col=lcol, lty=1,...)

plot(a$all.scale[,1], type="l", ylim=c(0,max(a$all.scale)), 
	xlab="iterations", ylab="estimates", main=tt2, ...)
# if (truevalues) abline(scale[1],0,lty=3)
title(sub=subtitle, cex.sub = 0.75)
lgd <- expression(sc[1]); lcol <- 1
for (j in 2:m) {
	lines(a$all.scale[,j], col=j, ...)
	# if (truevalues) abline(scale[j],0,col=j,lty=3)
	lgd <- c(lgd,substitute(sc[j])); lcol <- c(lcol,j)
	}
legend("topright", lgd, col=lcol, lty=1, ...)

plot(a$all.lambda[,1], type="l", ylim=c(0,1), 
	xlab="iterations", ylab="estimates", main=tt3, ...)
# if (truevalues) abline(lambda[1],0,lty=3)
title(sub=subtitle,cex.sub = 0.75)
lgd <- expression(lambda[1]); lcol <- 1
for (j in 2:m) {
		lines(a$all.lambda[,j], col=j, ...)
		# if (truevalues) abline(lambda[j],0,col=j,lty=3)
		lgd <- c(lgd,substitute(lambda[j]))
		lcol <- c(lcol,j)
		}
legend("topright", lgd, col=lcol, lty=1, ...)
	}
