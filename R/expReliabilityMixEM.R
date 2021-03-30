# code for Reliability Mixture Models (RMM) with Censored data
# D. Chauveau
# ref: Bordes L. and Chauveau D. Computational Statistics (2016)

# Simulate from an exponential mixture.  
# lambda = vector of component probabilities
#	rate = vector of component rates
rexpmix <- function(n,lambda=1,rate=1) {
	m <- length(lambda) # nb of components
	z <- sample(m, n, replace = TRUE, prob = lambda) # component indicator
	rexp(n, rate = rate[z])
}


# pdf of exponential mixture
dexpmixt <- function(t,lam,rate){
	m <- length(lam)
	f <- 0
	for (j in 1:m)
		f <- f + lam[j]*dexp(t,rate=rate[j])
	f
	}



##############################################
## EM algorithm for Reliability Mixture Models (RMM) with Censoring
## Parametric model, for 
## univariate finite mixture of exponentials with right censoring
# x = lifetime data, censored by random c if d is not NULL, in which case
# x= min(x,c) and d = I(x <= c)
# complete data may be (t,d,z) or (x,z), see option complete
expRMM_EM <- function (x, d=NULL, lambda = NULL, rate = NULL, k = 2, 
		  complete="tdz", epsilon = 1e-08, maxit = 1000, verb = FALSE) {
#  warn <- options(warn=-1) # Turn off warnings
#  x <- as.vector(x)
	n <- length(x)
	if (!is.null(lambda)) k=length(lambda)
	if (!is.null(rate)) k=length(rate) # either should be define !!!
	if (is.null(d))  d <- rep(1,n)  # simple EM for noncensored case
	xx <- matrix(x, nrow=n, ncol=k) # x repeated k times, used in E-steps
	dd <- matrix(d, nrow=n, ncol=k) # idem for d
# init function call to do later, in case lambda = rate = NULL  
	if (is.null(lambda)) lambda <- rep(1/k,k)
	if (is.null(rate)) rate <- rep(1,k)  # to do : init.rate(k)
	 # handle the 2 complete data models
	if (complete=="tdz") comptdz=TRUE else comptdz=FALSE
# sequences for storing along iterations 	
	lambda_seq <- rate_seq <- matrix(0, nrow = maxit, ncol = k)
 	lambda_seq[1,] <- lambda
 	rate_seq[1,] <- rate 
	loglik_seq <- NULL
	oldloglik <- -Inf
#    notdone <- TRUE
#    while(notdone) {
      # Initialize everything
      notdone <- FALSE
      dll <- epsilon+1
      iter <- 1
      post <- matrix(nrow = n, ncol = k)
      restarts <- 0     

      while (dll > epsilon && iter < maxit) { # EM iterations
		### E-step ###
		rr <- matrix(rate, n, k, byrow=T)
		ll <- matrix(lambda, n, k, byrow=T)
		# handling censored & non-censored cases
		post <- ((ll*dexp(xx,rr))^dd)*((ll*(1-pexp(xx,rr)))^(1-dd))
		rs <- rowSums(post)
		loglik <- sum(log(rs)) # loglik without the constant term
		# post normalized per row
		post <- sweep(post, 1, rs, "/") # posteriors p_{ij}^t 's

		### M-step ###
		lambda <- colMeans(post)
		lambda_seq[iter+1, ] <- lambda
		if (comptdz) rate <- colSums(post*dd)/colSums(post*xx) # complete=(t,d,z)
		# 			 rate <- colSums(post*d)/colSums(post*x) # gives same answer
		# 						cf automatic recycling
		if (!comptdz) {  # complete data = (x,z)
			mean_xi <- matrix(1/rate, n, k, byrow=T)
			rate <- colSums(post)/colSums(post*(xx + (1-dd)*mean_xi))
			# rate <- colSums(post)/colSums(post*(x + (1-d)*mean_xi)) # same answer
			}
		rate_seq[iter+1, ] <- rate

        dll <- loglik - oldloglik # = Inf for iter=0 1st time
        oldloglik <- loglik
        loglik_seq <- c(loglik_seq, loglik)
        if (verb) {
          cat("iteration", iter, " log-lik diff =", dll, " log-lik =", 
              loglik, "\n")
        }
        iter <- iter + 1
      }  # end of EM loops over iterations

      if (iter == maxit) cat("WARNING! NOT CONVERGENT!", "\n")   
    cat("number of iterations =", iter - 1, "\n")
	colnames(post) <- c(paste("comp", ".", 1:k, sep = ""))
    a=list(x=x, d=d, lambda = lambda, rate = rate, loglik = loglik, 
             posterior = post, all.loglik=loglik_seq, 
             all.lambda = lambda_seq[1:iter,],
             all.rate = rate_seq[1:iter,],
  #           restarts=restarts, 
              ft="expRMM_EM")
	class(a) = "mixEM"
  a
}





##################################################
# plot EM sequences from expRMM_EM:
# color by components, one plot per parameter type
# commented-out abline() usages are for true values for rate and lambda if available
plotexpRMM <- function(a, title=NULL, 
							rowstyle=TRUE, subtitle=NULL, ...)
	{
	n <- length(a$x); m <- dim(a$all.lambda)[2]
	pcc <- round(100*(1-mean(a$d)),2)
	sizes <- paste("n=",n,", ", pcc, "% censored", sep="")
	if (is.null(subtitle)) {
				subtitle <- paste("n=",n,", ", pcc, "% censored", sep="")}
	if (is.null(title)) {
		tt1 <- "Rate parameters" 
		tt2 <- "Weight parameters"
		} else tt1 <- tt2 <- title
lwdset <- 2
if (rowstyle) par(mfrow=c(1,2)) else par(mfrow=c(2,1))
plot(a$all.rate[,1], type="l", ylim=c(0,max(a$all.rate)), 
	xlab="iterations", ylab="estimates", main=tt1, ...)
title(sub=subtitle,cex.sub = 0.75)
lgd <- expression(xi[1]); lcol <- 1
for (j in 2:m) {
	lines(a$all.rate[,j], col=j, ...)
	# abline(rate[j],0,col=j,lty=3)
	lgd <- c(lgd,substitute(xi[j])); lcol <- c(lcol,j)
	}
legend("topright",lgd,col=lcol,lty=1,...)
plot(a$all.lambda[,1], type="l", ylim=c(0,1), 
	xlab="iterations", ylab="estimates", main=tt2, ...)
title(sub=subtitle,cex.sub = 0.75)
lgd <- expression(lambda[1]); lcol <- 1
if (m>2) {
	for (j in 2:m) {
		lines(a$all.lambda[,j], col=j, ...)
		lgd <- c(lgd,substitute(lambda[j]))
		lcol <- c(lcol,j)
		}
	}
legend("topright",lgd,col=lcol,lty=1,...)
	}

