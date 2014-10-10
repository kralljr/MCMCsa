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


	#first confirm first column is date
	if(class(dat[, 1]) != "Date") {
		stop("First column is not a date")
	}
	
	#separate dates and data
	dates <- dat[, 1]
	dat <- dat[, -1]
	
	
	#set dimensions
	P <- ncol(dat)
	T1 <- nrow(dat)
	
	
	#get indices correspond to id constraints
	cond <- which(!(is.na(lamcon)), arr.ind = T)

	#create arrays for output
	lamstar <- array(dim = c(L, P, N - burnin ))
	lfmat <-  array(dim = c(T1, L, N - burnin))
	# sigma2 <- array(dim = c(P, P, N - burnin))
	sigma2 <- array(dim = c(P, N - burnin))
	mu <- array(dim = c(L, N - burnin))
	xi2 <- array(dim = c(L, N - burnin))
	
	names1 <- c("lamstar", "lfmat", "sigma2", "mu", "xi2")
	
	#set first guesses
	if(is.null(guessvec)) {
		
		#truncated normal for lamstar
		lamstarG <- matrix(rtruncnorm(L * P, a = 0), nrow = L)
		lamstarG[cond] <- lamcon[cond]
		
		#lognormal for F
		lfmatG <- matrix((rnorm(T1 * L)), nrow = T1)
		# sigma2G <- diag(1, P)
		sigma2G <- rep(1, P)		
		muG <- rep(0, L)
		xi2G <- rep(1, L)
		guessvec <- list(lamstarG, lfmatG, sigma2G, muG, xi2G)
		names(guessvec) <- names1
	}
	
	
	#get lfhist
	lfhist <- list(length = T1)
	for(t in 1 : T1) {
		lfhist[[t]] <- lfmatG[t, ]
	}
	

	#for each iteration (N large)
	for (i in 1 : N) {
		print(i)
		
		#update all parameters
		for(j in 1 : length(guessvec)) {
			out <- gibbsfun(dat = dat, 
				guessvec = guessvec, lamcon = lamcon,
				type = names1[j], lfhist = lfhist, iter = i)
			lfhist <- out[["lfhist"]]
			guessvec <- out[["guessvec"]]
		}
		
		
		#if we are past the burnin period, save
		if (i > burnin) {
			lamstar[, , i - burnin] <- guessvec[["lamstar"]]
			lfmat[, , i - burnin] <- guessvec[["lfmat"]]
			# sigma2[, , i - burnin] <- guessvec[["sigma2"]]
			sigma2[, i - burnin] <- guessvec[["sigma2"]]
			mu[, i - burnin] <- guessvec[["mu"]]
			xi2[, i - burnin] <- guessvec[["xi2"]]
		}
			
			

		
	}
	
	#create output
	listout <- list(lamstar, lfmat, sigma2, mu, xi2)
	names(listout) <- names1
	
	
	listout
}



















##########################################
#update to new guess
##########################################
gibbsfun <- function(dat, guessvec, lamcon, type, lfhist, iter) {
	# print(type)
	if(type == "lamstar") {
		guessvec <- lamfun(dat, guessvec, lamcon)
		
	} else if(type == "lfmat") {
		
		out <- ffun(dat, guessvec, lfhist, iter)
		guessvec <- out[["guessvec"]]
		lfhist <- out[["lfhist"]]
		
	} else if(type == "mu") {
		guessvec <- mufun(dat, guessvec)
		
	} else if(type == "sigma2") {
		guessvec <- sig2fun(dat, guessvec)
		
	} else if(type == "xi2") {
		guessvec <- xi2fun(dat, guessvec)
		
	}else{
		stop("Error: type not recognized")	
	}		
		
	out <- list(guessvec, lfhist)
	names(out) <- c("guessvec", "lfhist")
	out
}
	
	
	



	
	
#####
# normal/inv-gamma variance sampling for sigma
# dat is data
# guessvec is list of guesses
sig2fun <- function(dat, guessvec)	{
	
	lfmat <- guessvec[["lfmat"]]
	fmat <- exp(lfmat)
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 2, colSums(lamstar), "/")
	
	#get dimensions
	P <- ncol(dat)
	T1 <- nrow(dat)

	#get mean
	mean <- log(fmat %*% lambda)
	
	#get priors
	pra <- 0.01
	prb <- 0.01
	
	#get posterior paramters
	a2 <- pra + T1 / 2
	diffs2 <- (dat - mean)^2
	b2 <- prb + colSums(diffs2) / 2
	
	#sample from inv gamma
	sigma2 <- rinvgamma(P, a2, b2)
	# sigma2 <- diag(sigma2)

	#update guess
	guessvec[["sigma2"]] <- sigma2
	
	guessvec
}
	
	
	
#####
# normal/inv-gamma variance sampling for xi
# dat is data
# guessvec is list of guesses
xi2fun <- function(dat, guessvec)	{
	
	lfmat <- guessvec[["lfmat"]]
	mu <- guessvec[["mu"]]
	
	#get dimensions
	L <- ncol(lfmat)
	T1 <- nrow(lfmat)
	
	#set prior values
	pra <- 0.01
	prb <- 0.01
	
	#get posterior paramters
	a2 <- pra + T1 / 2
	diffs2 <- (sweep(lfmat, 2, mu))^2
	b2 <- prb + colSums(diffs2) / 2
	
	#sample for all sources
	xi2 <- rinvgamma(L, a2, b2)

	#update guess
	guessvec[["xi2"]] <- xi2
	
	guessvec
}
	







#####
# normal/normal mean sampling
# dat is data
# guessvec is list of guesses
mufun <- function(dat, guessvec)	{
	

	xi2 <- guessvec[["xi2"]]
	lfmat <- guessvec[["lfmat"]]
	
	#get dimensions
	L <- ncol(lfmat)	
	T1 <- nrow(lfmat)
	
	#set prior values
	prmean <- 0
	prvar <- 100
	
	#get variance for all sources
	invars <- 1 / prvar + T1 / xi2
	vars <- 1 / invars
	
	#get mean for all sources
	num <- prmean / prvar + colSums(lfmat) / xi2
	mn <- num * vars
	
	#sample from independent normals
	mu <- rnorm(L, mean = mn, sd = sqrt(vars))
	
	#update guess
	guessvec[["mu"]] <- mu
	
	guessvec
}
		







