
#functions for data

#guessvec is guesses
#mdls is matrix of mdls for each day/constituent
#bdls is matrix of 1/0 for whether data BDL
yfun <- function(guessvec, mdls, bdls) {


	ldat <- guessvec[["ly"]]

	#set guesses
	fmat <- exp(guessvec[["lfmat"]])
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 1, rowSums(lamstar), "/")
	sigma2 <- guessvec[["sigma2"]]
	
	
	#get dimensions
	T1 <- nrow(fmat)
	
	#mean
	mn <- log(fmat %*% lambda)
	
	
	#log mdls
	lmdls <- log(mdls)
	minmdls <- min(lmdls) - 10

	# for each day with an observation below the MDL
	for (t in 1 : T1) {
		
		smiss <- sum(bdls[t, ])
		if(smiss > 0) {
			#draw a
			wh_miss <- which(bdls[t, ] > 0)

			temp <- rtruncnorm(1, 
				lower = minmdls, 
				upper = lmdls[t, wh_miss],
				mean = mn[t, wh_miss], 
				sd = sqrt(sigma2[wh_miss]))
				
			newy <- ldat[t, ]
			newy[wh_miss] <- temp
				
			guessvec[["ly"]][t, ] <- newy	
		}
		
					
	
	}#end loop over t	
				
	guessvec			
	
}







findmean <- function(ldat, mn, sigma2, sig.misscol, wh) {

	mnzero <- t(ldat - mn[t, ])

	#set up current guesses
	cov00 <- sigma2[wh, wh]
	cov01 <- sigma2[wh, -wh]
	
	mnobs <- mn[-wh]
	mnmiss <- mn[wh]
	yobs <- ldat[-wh]

	#get conditional mean and variance
	mn <- mnmiss + cov01 %*% sig.misscol[, , wh] %*% 
		mnzero[-wh]

	sd <- sqrt(cov00 - cov01 %*% sig.misscol[, , wh] %*% cov01)

	out <- list(mn, sd)
	names(out) <- c("mn", "sd")
	out
}
	
	
	
	




# function to get likelihood of y
# guessvec is guesses
# t is which day
# p is which constituent
logly <- function(guessvec, t = NULL, p = NULL) {
	
	
	if(type == "normal") {
		llhood <- logly.normal(guessvec, t, p)
		
	}else if(type == "mixture") {
		llhood <- logly.mixture(guessvec, t, p)
		
	}
	
	llhood
	
}



logly.normal <- function(guessvec, t, p) {
	ldat <- guessvec[["ly"]]
	lfmat <- guessvec[["lfmat"]]
	fmat <- exp(lfmat)
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 1, rowSums(lamstar), "/")
	sigma2 <- guessvec[["sigma2"]]
	
	#if want all
	if(is.null(t) & is.null(p)) {
		
		#get mean
		mn <- log(fmat %*% lambda)
	
		#get in exp
		diffsq <- (ldat - mn)
		sweeps <- sweep(diffsq, 2, sigma2 * 2, "/")
		llhood <- -rowSums(sweeps)
		
	# if want one day
	}else if(is.null(p)){
		
		#get data
		ldat <- ldat[t, ]
		#get mean
		mn <- log(t(fmat[t, ]) %*% lambda)
		llhood <- -sum((ldat - mn)^2 / (2 * sigma2))
		
	# if want one constituent
	}else if(is.null(t)) {
		
		#get data
		ldat <- ldat[, p]
		
		#get mean
		mn <- log(fmat %*% lambda[, p])
		llhood <- -sum((ldat - mn)^2 / (2 * sigma2[p]))
		
	#one constituent, one day	
	} else{
		
		ldat <- ldat[t, p]
		mn <- log(sum(fmat[t, ] * lambda[, p]))
		llhood <- -(ldat - mean)^2 / (2 * sigma2[p])
	}
	
	llhood
}






logly.mixture <- function(guessvec, bdls, t, p) {
	ldat <- guessvec[["ly"]]
	lfmat <- guessvec[["lfmat"]]
	fmat <- exp(lfmat)
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 1, rowSums(lamstar), "/")
	sigma2 <- guessvec[["sigma2"]]
	psi2 <- guessvec[["psi2"]]
	
	
		#if want all
	if(is.null(t) & is.null(p)) {
		T1 <- nrow(ldat)
		
		#get mean
		mn <- log(fmat %*% lambda)
	
		#get in exp
		diffsq <- (ldat - mn)
		
		vars1 <- matrix(rep(sigma2, T1), byrow = T, nrow = T1)
		vars2 <- matrix(rep(psi2, T1), byrow = T, nrow = T1)
		vars <- vars1 + vars2 * bdls
		
		llhood <- -rowSums(diffsq / (vars * 2))
		
	# if want one day
	}else if(is.null(p)){
		
		#get data
		ldat <- ldat[t, ]
		#get mean
		mn <- log(t(fmat[t, ]) %*% lambda)
		
		vars <- sigma2 + psi2 * bdls[t, ]
		llhood <- -sum((ldat - mn)^2 / (2 * vars))
		
	# if want one constituent
	}else if(is.null(t)) {
		
		#get data
		ldat <- ldat[, p]
		
		#get mean
		mn <- log(fmat %*% lambda[, p])
		vars <- sigma2[p] + bdls[, p] * psi2[p]
		
		llhood <- -sum((ldat - mn)^2 / (2 * vars))
		
	#one constituent, one day	
	} else{
		
		ldat <- ldat[t, p]
		mn <- log(sum(fmat[t, ] * lambda[, p]))
		
		vars <- sigma2[p] + psi2[p] * bdls[t, p]
		llhood <- -(ldat - mean)^2 / (2 * vars)
	}
	
	llhood
	
	
}
	