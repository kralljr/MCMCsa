#' MCMC for source apportionment
#'
#' \code{mcmcsa} performs MCMC for source apportionment
#'
#' This function aims to draw from the posterior distributions of the source
#' profiles, the source contributions/concentrations, the variance of the data
#' and the other parameters.  Eventually, this function will incorporate 
#' options for concentrations observed below the MDL
#'
#' @param data data frame of daily constituent concentrations with date as first column
#' @param L number of sources
#' @param lamcon Matrix of sources by constituents with 1 and 0 for conditions
#' @param mdls data frame of daily MDL values of the same dimensions of data
#' @param guessvec starting values (do not need to provide)
#' @param burnin number of burnin samples
#' @param N number of samples
#' @export
#' @examples
#' data(nycdat)
#' data(nycmdl)
#' #fix data totals
#' pm <- nycdat$PM25
#' whPM <- which(colnames(nycdat) == "PM25")
#' nycdat <- nycdat[, -whPM]
#' whPM <- which(colnames(nycmdl) == "PM25")
#' nycmdl <- nycmdl[, -whPM]
#' mcmcsa(nycdat, nycmdl)
mcmcsa <- function(dat, L, lamcon, mdls = NULL, 
	guessvec = NULL, burnin = 10000, N = 100000){

	if(class(dat[, 1] != "Date")) {
		stop("First column is not a date")
	}
	
	dates <- dat[, 1]
	dat <- dat[, -1]
	
	P <- ncol(dat)
	T1 <- nrow(dat)
	
	
	#get indices to replace
	cond <- which(!(is.na(lamcon)), arr.ind = T)

	#create arrays for output
	lamstar <- array(dim = c(L, P, N - burnin ))
	fmat <-  array(dim = c(T1, L, N - burnin))
	sigma2 <- array(dim = c(P, P, N - burnin))
	mu <- array(dim = c(L, N - burnin))
	xi2 <- array(dim = c(L, N - burnin))
	
	names1 <- c("lamstar", "fmat", "sigma2", "mu", "xi2")
	
	#set first guesses
	if(is.null(guessvec)) {
		
		temp <- matrix(rtruncnorm(L * P, a = 0), nrow = L)
		#fix set
		temp[cond] <- lamcon[cond]
		lamstar[,, 1] <- temp
		
		temp <- matrix(rnorm(T1 * L), nrow = T1)
		fmat[,, 1] <- exp(temp)
		sigma2[,, 1] <- diag(1, P)
		mu[, 1] <- rep(0, L)
		xi2[, 1] <- rep(1, L)
		guessvec <- list(lamstar, fmat, sigma2, mu, xi2)
		names(guessvec) <- names1
	}
	

	#for each iteration (N large)
	for (i in 1 : N) {
		
		#update all parameters
		for(j in 1 : length(guessvec)) {
			guessvec <- gibbsfun(dat = dat, guessvec = guessvec,
				type = names1[j])
		}
		
		
		#if we are past the burnin period, save
		if (i > burnin) {
			lamstar[, , i - burnin] <- guessvec[["lamstar"]]
			fmat[, , i - burnin] <- guessvec[["fmat"]]
			sigma2[, , i - burnin] <- guessvec[["sigma2"]]
			mu[, i - burnin] <- guessvec[["mu"]]
			xi2[, i - burnin] <- guessvec[["xi2"]]
		}
			
			

		
	}
	
	#create output
	listout <- list(lamstar, fmat, sigma, mu, xi)
	names(listout) <- names1
	
	
	listout
}



















##########################################
#update to new guess
##########################################
gibbsfun <- function(dat, guessvec, type) {
		
	if(type == "lamstar") {
		guessvec <- lamfun(dat, guessvec)
		
	} else if(type == "fmat") {
		
		guessvec <- ffun(dat, guessvec)
		
	} else if(type == "mu") {
		guessvec <- mufun(dat, guessvec)
		
	} else if(type == "sigma2") {
		guessvec <- sig2fun(dat, guessvec)
		
	} else if(type == "xi2") {
		guessvec <- xi2fun(dat, guessvec)
		
	}else{
		stop("Error: type not recognized")	
	}		
		
	guessvec
	}
	
	
	



	
	
#####
#check this
# normal/inv-gamma variance sampling for sigma
sig2fun <- function(dat, guessvec, dim, name)	{
	
	P <- ncol(dat)
	fmat <- guessvec[["fmat"]]
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 2, colSums(lamstar), "/")

	mean <- log(fmat %*% lambda)
	
	pra <- 0.01
	prb <- 0.01
	
	#for each constituent
	for(p in 1 : P) {
		sigmap <- invg(dat[, p], mean[, p], pra, prb)
		
		guessvec[["sigma2"]][p, p] <- sigmap
		
	}
	
	guessvec
}
	
	
	
#####
# normal/inv-gamma variance sampling for xi
xi2fun <- function(dat, guessvec)	{
	
	L <- length(guessvec[["xi2"]])
	lfmat <- log(guessvec[["fmat"]])
	mu <- guessvec[["mu"]]
	
	pra <- 0.01
	prb <- 0.01
	
	#for each constituent
	for(l in 1 : L) {
		xi2 <- invg(lfmat[, l], mu[l], pra, prb)
		
		guessvec[["xi2"]][l] <- xi2
		
	}
	
	guessvec
}
	




#### get inverse gamma
#b is scale
invg <- function(x, mu, pra, prb) {
	n <- length(x)
	
	a2 <- pra + n/2
	b2 <- prb + sum((x - mu)^2) / 2
	
	rinvgamma(1, a2, b2)
}




#####
# normal/normal mean sampling
mufun <- function(dat, guessvec)	{
	
	L <- length(guessvec[["mu"]])
	xi2 <- guessvec[["xi2"]]
	lfmat <- log(guessvec[["fmat"]])
	
	prmean <- 0
	prvar <- 100
	
	for(l in 1 : L) {
		
		vars <- 1/(1/prvar + nrow(lfmat)/xi2[l])
		
		num <- (prmean/prvar + sum(lfmat[, l])/xi2[l])
		mn <- num * vars

		mu <- rnorm(1, mean = mn, sd = sqrt(vars))

		guessvec[["mu"]] <- mu
		
	}
	
	guessvec
}
		







