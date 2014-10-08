
#####
# adaptive mcmc for F?
ffun <- function(dat, guessvec, iter, mnf, covf)	{
	
	T1 <- nrow(dat)
	L <- length(guessvec[["mu"]])
	
	
	for(t in 1 : T1) {
		mn <- log(guessvec[["fmat"]][t, ])
		
		sk <- 2.4^2 / L
		C0 <- diag(0.39, L)
		covf <- 
		
		C1 <- sk * covf + sk * eps * diag(L)
		cov1 <- ifelse((iter + 1) <= 15, C0, C1)
		
		newf <- rmvnorm(1, mean = mn, sigma = cov1)
		
		lhoodOLD <- lhoodf(mn, dat, guessvec, t)
		lhoodNEW <- lhoodf(newf, dat, guessvec, t)
		
		prob <- min(c(1, lhoodNEW / lhoodOLD))
		sel <- rbinom(1, 1, prob)
		
		newf <- ifelse(sel == 1, newf, mn)
		
		#may need to have normal bounded between -10000, 10000 for compact set
		
	}
	
	guessvec
}





lhoodf <- function(fmat, dat, guessvec, t) {
	
	#likelihood of data
	guessvec[["fmat"]] <- fmat
	first <- ly(dat, guessvec, t)
	
	
	#likelihood of f
	lfmat <- log(guessvec[["fmat"]][t, ])
	xi2 <- guessvec[["xi2"]]
	mu <- guessvec[["mu"]]
	second <- -sum(1/ (2 * xi2) * ((lfmat - mu) ^ 2))
	second <- exp(second)
	
	first * second
}




ly <- function(dat, guessvec, t = NULL, p = NULL) {
	ldat <- log(dat)
	
	fmat <- guessvec[["fmat"]]
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 2, colSums(lamstar), "/")
	sigma2 <- guessvec[["sigma2"]]
	
	
	if(is.null(t) & is.null(p)) {
		mn <- log(fmat %*% lambda)
	
		diffsq <- (ldat - mn)
		sweeps <- sweep(diffsq, 2, diag(sigma2) * 2, "/")
		rs <- -rowSums(sweeps)
		first <- exp(rs)
	}else if(is.null(p)){
		ldat <- ldat[t, ]
		mn <- log(t(fmat[t, ]) %*% lambda)
		first <- -sum((ldat - mn)^2 / (2 * diag(sigma2)))
		first <- exp(first)
	}else if(is.null(t)) {
		
		ldat <- ldat[, p]
		mn <- log(fmat %*% lambda[, p])
		eachday <- (ldat - mn)^2 / (2 * sigma2[p, p])
		first <- exp(-sum(eachday))
	} else{
		
		ldat <- ldat[t, p]
		mn <- log(sum(fmat[t, ] * lambda[, p]))
		first <- exp(-(ldat - mean)^2 / (2 * sigma2[p, p]))
	}
	
	first
	
}
	


