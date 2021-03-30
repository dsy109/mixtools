###########################################################
# semi-parametric Reliability Mixture Models with Censored data
#	for 2 components                                     
# D. Chauveau
# ref: Bordes L. and Chauveau D. Computational Statistics (2016)
###########################################################

# Simulate from a lognormal scale mixture
# like model (3) in compstat paper
# lambda = vector of component probabilities
#	xi = scaling parameter
rlnormscalemix <- function(n, lambda=1, 
			meanlog=1, sdlog=1, scale=0.1) {
	m <- length(lambda) # nb of components
	z <- sample(m, n, replace = TRUE, prob = lambda) # component indicator
	x <- rlnorm(n, meanlog, sdlog)
	x[z==2] <- x[z==2]/scale
	x
}


##############################################
# Kaplan-Meier survival estimate 
# vectorized version ONLY on ** already ordered data **
# cpd = (n,2) matrix where
# cpd[,1] = t = ORDERED censored life data = min(x,c)
# cpd[,2] = d = censoring indicator, 1=observed 0=censored
# returned value = a step function object allowing to evaluate S()
# 	at any point, needed by SEM for computing S(xi*t_i)'s
KMod <- function(cpd, already.ordered=TRUE){
	# order stat of t and d ordered accordingly
	if (already.ordered) {
		n <- dim(cpd)[1]
		s <- cumprod(1 - cpd[,2]/(n - (1:n) +1))
		stepfun(cpd[,1], c(1,s))} # starts with 1 for t<=t(1)
	}


##############################################
## KM estimate integration
# s = a stepfun object as returned by KMod()
# returned value = int_0^max(t_i) S(u)du
KMintegrate <- function(s) {
	ks <- knots(s)
	n <- length(ks) # number of knots=data points
	ks2 <- c(0,ks[1:(n-1)]) # shifted right
	hs <- c(1,s(ks)[1:(n-1)]) # heights of survival
	sum((ks-ks2)*hs)
	}



################################################
###    HAZARD RATE SMOOTH KERNEL ESTIMATES   ###
################################################
# WKDE with Triangular Kernel and global bandwidth, vectorized
# returns vector of K_h(t_i-u_j), j=1,...,length(u)
# t = vector of data
# u = vector of points at which K(u) is computed
# w = weights defaults to 1/n, not forced to be normalized
# bw=bandwidth, can be a n-vector or a scalar
# nb: this is not specific for hazard rate, just a wkde 
# NB: triang_wkde is the default wkde called by HRkde() for hazard rate 
triang_wkde <- function(t, u=t, w=rep(1/length(t),length(t)), 
						bw=rep(bw.nrd0(t),length(t))) {
	n <- length(t); p <- length(u)
	xx <- outer(t, u, function(a,b) a-b) # K[i,j] = (t_i-u_j)
	xx <- abs(xx) # not necessarily for all kernels, separated from outer
	h <- matrix(bw,n,p) # bw repeated over p columns
						# works also if a scalar bw is passed!
	xx <- xx/h  # now xx is |t_i-u_j|/h_i
	K <- (xx <= 1)*(1-xx)/bw  # K((t-u)/h)/h
	mw <- matrix(w, nrow=1) # weights for matrix product
	as.vector(mw %*% K)
	}





## Hazard Rate kernel density estimate based on *ordered* censored data
# weighted kernel density estimate passed as function argument, 
# defaults to symmetric triangle kernel 
# cpd = (n,2) matrix where
# cpd[,1] = t = ORDERED censored life data = min(x,c)
# cpd[,2] = d = censoring indicator, 1=observed 0=censored
# u = vector of points at which alpha() is evaluated, 
#     defaults to t itself
# bw = bandwidth *vector* for the kernel density estimate
# kernelft = kernel definition
# return value = alpha(.) evaluated at u_j's
# NB: ordered data are necessary from weights definitions
#     (n-i+1) items at risk etc
HRkde <- function(cpd, u = cpd[,1], kernelft = triang_wkde, 
					bw = rep(bw.nrd0(as.vector(cpd[,1])), length(cpd[,1]))){
	# gaussian case (not recommended)
	# if (kernelft==gaussian) HRkdeGauss(cpd,u,bw) else {
	n <- length(cpd[,1])
	aw <- cpd[,2]/(n - (1:n) +1) # weights d(i)/(n-i+1) 
	kernelft(cpd[,1], u, w=aw, bw=bw)
# old, non vectorized version (deprecated)	
#	nu <- length(u)
#		hr <- rep(0,nu)
#		for (k in 1:nu) {
#			K <- kernelft(cpd[,1], u[k], bw)
#			hr[k] <- sum(K*aw)
#			}
#		}
#	hr
	}
	


#################################################################
#################################################################
## Stochastic EM algorithm for semiparametric Scaling
## Reliability Mixture Models (RMM) with Censoring, 2 components
# t = lifetime data, censored by random c if d is not NULL, in which case
# t = min(x,c) and d = I(x <= c)
# rate = scaling parameter
# centers = initial centers for initial call to kmeans
# averaged = TRUE if averaging performed at each iteration (cf Nielsen 2000)
# 	NB: averaging can be done in (at least) 2 ways; 
#	here the current theta for E & M steps is the average over the sequence,
#	but the theta^next not averaged is stored on the sequence 
# batchsize = number of last iterations to use for "averaging" unscaled samples 
#   for computing "average" final KM and HR estimates 
#	alpha() and S() 
#  kernelft = kernel used in HRkde for hazard rate nonparametric estimate
#
# NB: since the sequence of parameters is a stoch. process, the
# stopping criterion is usually not satisfied until maxit 

spRMM_SEM <- function (t, d = NULL, lambda = NULL, scaling = NULL, 
          centers = 2, kernelft = triang_wkde, 
          bw = rep(bw.nrd0(t),length(t)), averaged = TRUE,
          epsilon = 1e-08, maxit = 100, batchsize = 1, verb = FALSE) {
	k = 2 # fixed number of components for this model identifiability
	n <- length(t)
	if (is.null(d))  {
		d <- rep(1,n)  # but better use specific algo for noncensored data
		cat("warning: undefined censoring indicator d replaced by 1's")
		cat(" i.e. all data are assumed observed")
		cat(" better use instead a specific St-EM algorithm for uncensored data")
		}
									
	##### Initializations #####	
	# init function call to do later, in case lambda = scaling = NULL  
	if (is.null(lambda)) lambda <- rep(1/k,k) # not USED for init!?
	if (is.null(scaling)) scaling <- 1/2      # ToDo better init using centers fro kmeans!
	# sequences for storing along iterations 	
	lambda_seq <- matrix(0, nrow = maxit, ncol = k)
 	lambda_seq[1,] <- lambda
 	scaling_seq <- rep(0,maxit)
 	scaling_seq[1] <- scaling 
 	sumNaNs <- 0 # for total nb of NaN's in iterations while computing post
    # dll <- epsilon+1 # not applicable for SEM versions, always run maxit
    post <- z <- sumpost <- posthat <- matrix(0, nrow = n, ncol = k)
    
    qt=round(maxit/10); pc <- 0 # for printing of % of iteration done
	if (qt == 0) qt <- 1

    # initialization for batch storing
    if (batchsize > maxit) batchsize <- maxit
    batch_t <- NULL # for storing batchsize last unscaled samples
    batch_d	<- NULL # and associated event indicators

	# kmeans method for initial posterior and z matrix
	kmeans <- kmeans(t, centers=centers)
	for(j in 1:k) {
   		z[kmeans$cluster==j, j] <- 1
  		}
    # init values for Survival s() and hazard rate a()
	# strategy: use a M-step with z from kmeans result 
	zt <- z*t # recycling t m times => zt[i,] holds  
	zt[,2] <- scaling*zt[,2] # t if comp 1 and xi*t if comp 2
	newt <- rowSums(zt) # redressed "unscaled" sample
	newo <- order(newt)
	new.od <- cbind(newt,d)[newo,] # ordered data with associated d
	s <- KMod(new.od) # stepfun object, can be evaluated by s(t)
	# note: <=> to use KMsurvfit(newt,d) which calls survival package functions
	a1 <- HRkde(new.od, t, kernelft, bw=bw) # evaluated at the original t's
	a2  <- HRkde(new.od, scaling*t, kernelft, bw=bw) # and at scaling*t's

	####### SEM iterations #######
	iter <- 1
    # while (dll > epsilon && iter < maxit) { # EM-only version
    while (iter < maxit) {
	### E-step
		post[,1] <- ((lambda[1]*a1*s(t))^(d))*((lambda[1]*s(t))^(1-d))
		post[,2] <- (lambda[2]*scaling*a2*s(scaling*t))^d # observed
		post[,2] <- post[,2]*(lambda[2]*s(scaling*t))^(1-d) # censored
		# print(post)		
		rs <- rowSums(post) # post normalized per row
		post <- sweep(post, 1, rs, "/") # posteriors p_{ij}^t 's
		# check and solve NaN's: may cause theoretical pbs!?
		snans <- sum(is.na(post))
		if (snans > 0) {
			post[is.na(post[,1]),] <- 1/k # NaN's replaced by uniform weights
			sumNaNs <- sumNaNs + snans
			if (verb) cat(snans, "NaN's in post: ")
			}
		sumpost <- sumpost + post
				
	### S-step
	# ~ matrix of component indicators simu checked OK
		z <- t(apply(post, 1, function(prob) rmultinom(1, 1, prob))) 
		# cbind(post,z) # checking simulation
		zt <- z*t # each row of z is 0,1 hence zt holds {0,t}
		nsets <- colSums(z) # subsets per component sizes
	
	### M-step for scalar parameters		
		newlambda <- nsets/n # or colMeans(post) if EM strategy preferred!
		lambda_seq[iter+1, ] <- newlambda
		
	# update of the scaling (xi) parameter needs building of KM estimates
	# from each subsample, and integration of both
	# ToDo: use parallel on 2 cores to perform each subsample task!
		tsub1 <- t[z[,1]==1]; tsub2 <- t[z[,2]==1] # t subsamples
		dsub1 <- d[z[,1]==1]; dsub2 <- d[z[,2]==1] # d subsamples
		o1 <- order(tsub1); od1 <- cbind(tsub1,dsub1)[o1,]
		o2 <- order(tsub2); od2 <- cbind(tsub2,dsub2)[o2,]
		km1 <- KMod(od1); km2 <- KMod(od2)
		newscaling <- KMintegrate(km1)/KMintegrate(km2)
		scaling_seq[iter+1] <- newscaling # stored for final plotting
		if (averaged) {
			scaling <- mean(scaling_seq[1:(iter+1)])
			lambda <- colMeans(lambda_seq[1:(iter+1),])
			} else {	# averaged=FALSE case, just use last update
				scaling <- newscaling
				lambda <- newlambda
				} 		
		
	    ### M-step for nonparametric alpha() and s():
		# unscaling the sample and order it,
		# keeping the associated d censoring indicator
		zt[,2] <- scaling*zt[,2] # unscale t's from component 2
		newt <- rowSums(zt) # "unscaled" sample (sum 0 and t_i or scale*t_i )
		newo <- order(newt)
		new.od <- cbind(newt,d)[newo,] # with associated (unchanged) d
		s <- KMod(new.od) # stepfun object, can be evaluated by s(t)
	  	a1 <- HRkde(new.od, t, kernelft, bw=bw) # evaluated at the original t's
	  	a2  <- HRkde(new.od, scaling*t, kernelft, bw=bw) # and at scaling*t's
		
		### batch storage: collecting unscaled samples 
		if (iter >= (maxit - batchsize)) {
			batch_t <- c(batch_t, newt)
    		batch_d	<- c(batch_d, d)
			cat("-- adding unscaled sample at iteration",iter,"\n")
			}
		
		# change <- c(lambda_seq[iter,] - lambda_seq[iter-1,], 
        #	scaling_seq[iter]- scaling_seq[iter-1])
      	# dll <- max(abs(change))
      	if (verb)  {
      		cat("iteration", iter, ": sizes=", nsets," ")
      		cat("lambda_1=",lambda[1], "scaling=",scaling, "\n")
      		}
      			
      	# printing % done
		b=round(iter/qt); r=iter-b*qt
		if (r==0) {pc <- pc+10; cat(pc,"% done\n")}

      	iter <- iter + 1
      }  ###### end of SEM loops over iterations ######
    
	colnames(post) <- colnames(posthat) <- c(paste("comp", ".", 1:k, sep = ""))
	lambdahat <- colMeans(lambda_seq[1:iter,])
	scalinghat <- mean(scaling_seq[1:iter])
	posthat <- sumpost/(iter-1) # average posteriors
	hazard <- HRkde(new.od, kernelft=kernelft, bw=bw) # at final unscaled sample
	# Would it be better to return the final unscaled sample??
	
	cat(iter, "iterations: lambda=",lambdahat,", scaling=",scalinghat,"\n")
	if (sumNaNs > 0) cat("warning: ", sumNaNs,"NaN's occured in posterior computations")
	
	###### FINISHING ######
	# ToDo: average over batch set of last unscaled samples
	# compute average estimates via a S+M-step using posthat etc
	loglik <- sum(log(rs)) 	# "sort of" loglik with nonparam densities estimated
							# mimmics the parametric case, just like npEM does 
	z <- t(apply(posthat, 1, function(prob) rmultinom(1, 1, prob)))
	zt <- z*t
	zt[,2] <- scalinghat*zt[,2]
	avgt <- rowSums(zt) # unscaled sample
	avgo <- order(avgt)
	avg.od <- cbind(avgt,d)[avgo,] # with associated d
	shat <- KMod(avg.od) # stepfun object
	ahat <- HRkde(avg.od, kernelft=kernelft, bw=bw) # eval at the unscaled ordered avgt

    a=list(t=t, d=d, lambda=lambdahat, scaling=scalinghat, 
          posterior=post,     # final posterior probabilities
          all.lambda=lambda_seq[1:iter,],  # sequence of lambda's
          all.scaling=scaling_seq[1:iter], # sequence of scale
          loglik=loglik,	  # analog to parametric loglik, like npEM
          meanpost=posthat,   # posterior proba averaged over iterations
          survival=s,         # Kaplan-Meier last iter estimate (stepfun object)
          hazard=hazard,      # hazard rate final estimate evaluated at final.t 
          final.t=new.od[,1], # last unscaled sample
          s.hat=shat,         # Kaplan-Meier average estimate
          t.hat= avg.od[,1],  # ordered unscaled sample based on meanpost USEFUL??
          avg.od=avg.od,	  # t.hat with associated d (for computing ahat outside)
		  hazard.hat=ahat,    # hazard rate average estimate on t.hat
		  batch.t=batch_t,    # batch sample t (not ordered)
		  batch.d=batch_d,    # associated event indicators just rep(d,batchsize)
		  sumNaNs = sumNaNs,  # total nb of NaN's in iterations while computing post
          ft="spRMM_SEM")
	class(a) = "spRMM"
  return(a)
}




##############################################
## plot function for spRMM object : remove lognormal true pdf!
# sem = spRMM_SEM object
# other are true parameters, if available, for plotting references
# ToDo: plot true f and S only if true param availables
# 		pass the true pdf (like dlnorm) as well
plotspRMM <- function(sem, tmax = NULL){
	t <- sem$t   
	ym <- max(sem$all.scaling)
	par(mfrow=c(2,2))
	plot(sem$all.scaling, type="l", ylim=c(0, ym),
		xlab="iterations", main="scaling", ylab="")
	plot(sem$all.lambda[,1], ylim=c(0,1), type="l", xlab="iterations",
		main="weight of component 1",ylab="")
	# plots of Kaplan-Meier estimates
	# finding max time for plots
	if (is.null(tmax)){tmax <- max(sem$scaling*t) + 2}
	u <- seq(0, tmax, len=200)
	plot(sem$survival, xlim=c(0,tmax), ylim=c(0,1), pch=20, 
		verticals=F, do.points=F, xlab="time", main="Survival function estimate")
	# plot(sem$s.hat, pch=20, verticals=F, do.points=F, col=2, add = TRUE)
	# pdf estimates
	fhat <- sem$s.hat(sem$t.hat)*sem$hazard.hat
	ffinal <- sem$survival(sem$final.t)*sem$hazard
  plot(sem$final.t, ffinal, type="l", col=1, xlim=c(0,tmax), 
				xlab="time", ylab="density", main="Density estimate")
}



## S3 method of summary for class "spRMM"
summary.spRMM <- function(object, digits = 6, ...) {
  sem <- object
  if (sem$ft != "spRMM_SEM") stop("Unknown object of type ", sem$ft)
  cat("summary of", sem$ft, "object:\n")
  o <- matrix(NA, nrow = 2, ncol = 2)
  o[1,] <- sem$lambda
  o[2,2] <- sem$scaling
  colnames(o) <- paste("comp",1:2)
  rownames(o) <- c("lambda", "scaling")
  print(o, digits = digits, ...)
  cat("(pseudo) loglik at estimate: ", sem$loglik, "\n")
  pcc <- round(100*(1-mean(sem$d)), 2)
  cat(pcc, "% of the data right censored\n")
}



