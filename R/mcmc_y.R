
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
	
	#get inverse of gsig
	siginv <- chol2inv(chol(sigma2))
	
	#which columns have missing
	wcolmiss <- which(colSums(bdls) != 0)
	sig.misscol <- array(dim = c(nrow(sigma2) - 1, 
		ncol(sigma2) - 1, ncol(sigma2)))
		
	#for each column with missingness	
	for( i in 1 : length(wcolmiss)) {
		num <- wcolmiss[i]
		sigcut <- siginv[-num, -num] - siginv[num, -num] %*% 
			t(siginv[-num, num]) / siginv[num, num]
		#update array	
		sig.misscol[, , num] <- sigcut
	}
	
	
	# for each day with an observation below the MDL
	for (t in 1 : T1) {
		

		#if which missing greater than 0
		if(length(wh_miss[[t]]) > 0) {
		
			#for each missing value on day i
			for (j in 1 : length(wh_miss[[t]])) {
				
				datc <- ldat[t, ]
		
				#find conditional mean/var
				mnv <- findmean(dat = datc, 
					mean = mn[t, ], 
					cov1 = sigma2, 
					sig.misscol = sig.misscol, 
					wh = wh_miss[[t]][j])
	
				#propose truncated normal (0 to mdl)	 
					#log scale with cond mean and var		
				newy <- rtnorm(1, 
					lower = min(lmdls) - 10, 
					upper = lmdls[t, wh_miss[[t]][j]],
					mean = mnv[["mn"]], sd = mnv[["sd"]])
					
				if( is.na(newy)) {browser()}	
		
				#update guess
				ldat[t, wh_miss[[t]][j]] <- newy
				guessvec[["ly"]][t, wh_miss[[t]][j]] <- newy
		
			}#end loop over j
				
				
		}#end check
	
	
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
	
	#get guesses
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
	