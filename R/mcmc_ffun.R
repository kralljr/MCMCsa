
#####
# adaptive mcmc for F
# guess vec is vector of guesses
# iter is current iteration
# lfhist is list of length T where each item is iter X L matrix
ffun1 <- function(guessvec, lfhist, iter, bound = 80)	{
	
	lfmat <- guessvec[["lfmat"]]
	
	#set dimentions
	T1 <- nrow(lfmat)
	L <- length(guessvec[["mu"]])

	#if next iteration will be change over, set up
	if((iter + 1) == 15) {
		lfhist2 <- list()
		lfhist2[[1]] <- array(dim = c(L, L, T1))
		lfhist2[[2]] <- array(dim = c(L, T1))
		lfhist2[[3]] <- array(dim = c(L, T1))
		names(lfhist2) <- c("C1", "Xt0", "Xt1")
	}
	
	
	#values from Hackstadt 2014
	sk <- 2.4^2 / L
	eps <- 0.001
	C0 <- diag(0.39, L)
	
	
	if((iter + 1) <= 15) {
		C1 <- array(C0, dim = c(L, L, T1))
	}else{
		C1 <- lfhist[["C1"]]
		}
	
	#for each day
	for(t in 1 : T1) {
		
		
		#get mean
		mn <- lfmat[t, ]
	
		#make sure covariance is PD
		tempC1 <- make.positive.definite(C1[,, t])
		
		#sample LFs
		if(is.positive.definite(tempC1)) {
			newlf <- rmvnorm(1, mean = mn, sigma = tempC1)
			lhoodNEW <- loglhoodf(newlf, guessvec, t, bound)
		}else{
			newlf <- rmvnorm(1, mean = mn, sigma = tempC1, method = "svd")
			
			if(length(which(is.na(newlf))) > 0) {
				print("sigma not pd")
				lhoodNEW <- -Inf
			}else{
				print("svd")
				lhoodNEW <- loglhoodf(newlf, guessvec, t, bound)
				}
			
		}

		lhoodOLD <- loglhoodf(mn, guessvec, t, bound)
		
		
		
		#acceptance
		lprob <- min(0, lhoodNEW - lhoodOLD)
		lprob <- ifelse(is.na(lprob), -Inf, lprob)
		unif1 <- runif(1)
		
		#update guess
		if(log(unif1) <= lprob) {
			if(min(newlf) < -bound | max(newlf) > bound) {
				browser()
			}
			guessvec[["lfmat"]][t, ] <- newlf
		}else{
			#replace current for newlf
			newlf <- t(matrix(mn))
			
		}
				
		#bound Normal between -b, b for compact set?
		
		#record history
		if((iter + 1) < 15) {
			lfhist[[t]] <- rbind(lfhist[[t]], newlf)
			
		#do first update	
		}else if((iter + 1) == 15){
			lcurr <- rbind(lfhist[[t]], newlf)
			cov1 <- cov(lcurr)
			lfhist2[["C1"]][,, t] <- sk * cov1 + sk * eps * diag(L)
			Xt0 <- apply(lfhist[[t]], 2, mean)
			lfhist2[["Xt0"]][, t] <- Xt0
			lfhist2[["Xt1"]][, t] <- (Xt0 * (iter - 1) + newlf) / iter

		#update rest based on recursion		
		}else{
			
			newlf <- t(newlf)
			#mean of two back will be mean of one back
			Xt0 <- lfhist[["Xt1"]][, t, drop = F]
			lfhist[["Xt0"]][, t] <- Xt0
			# add in current
			Xt1 <- (Xt0 * (iter - 1) + newlf) / iter
			lfhist[["Xt1"]][, t] <- Xt1
			
			
			midd <- iter * Xt0 %*% t(Xt0) - (iter + 1) * Xt1 %*% t(Xt1)
			midd <- midd + newlf %*% t(newlf) + eps * diag(L)
			tempC1 <- (iter - 1)/iter * C1[,, t] + sk/iter * midd
			
			lfhist[["C1"]][,, t] <- tempC1
			

		}
		
			
			
			

		
	}#end loop over t
	
	
	if((iter + 1) == 15) {
		lfhist <- lfhist2
	}
	
	#print(min(guessvec[["lfmat"]]))
	out <- list(guessvec, lfhist)
	names(out) <- c("guessvec", "lfhist")
	out
}


ffun <- function(guessvec, lfhist, iter, bound = 80)	{
	
	lfmat <- guessvec[["lfmat"]]
	
	#set dimentions
	T1 <- nrow(lfmat)
	L <- length(guessvec[["mu"]])
	
	#for each day
	for(t in 1 : T1) {
		
		
		#get mean
		mn <- lfmat[t, ]

		newlf <- rnorm(L, mean = mn, sd = 0.1)
		lhoodNEW <- loglhoodf(newlf, guessvec, t, bound)
		lhoodOLD <- loglhoodf(mn, guessvec, t, bound)
		
		
		
		#acceptance
		lprob <- min(0, lhoodNEW - lhoodOLD)
		lprob <- ifelse(is.na(lprob), -Inf, lprob)
		unif1 <- runif(1)
		
		#update guess
		if(log(unif1) <= lprob) {
			if(min(newlf) < -bound | max(newlf) > bound) {
				browser()
			}
			guessvec[["lfmat"]][t, ] <- newlf
		}
				
	}

	
	#print(min(guessvec[["lfmat"]]))
	out <- list(guessvec, lfhist)
	names(out) <- c("guessvec", "lfhist")
	out
}
	
	
	
	
	

#lhood for F
# fmat is current guess for F
# guessvec is list of guesses
# t is day
loglhoodf <- function(lfmat, guessvec, t, bound) {
	
	#likelihood of data
	guessvec[["lfmat"]][t, ] <- lfmat
	ly <- logly(guessvec, t = t)
	
	
	#likelihood of f
	lfmat <- guessvec[["lfmat"]][t, ]
	xi2 <- guessvec[["xi2"]]
	mu <- guessvec[["mu"]]
	
	
	if(min(lfmat) > -bound & max(lfmat) < bound) {
		#note lfmat is only vector
		lf <- -sum(1/ (2 * xi2) * ((lfmat - mu) ^ 2))
	}else {
		lf <- -Inf
	}
	
	ly + lf
}





