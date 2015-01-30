#' MCMC for source apportionment
#'
#' \code{mcmcsa} performs MCMC for source apportionment
#'
#' This function aims to draw from the posterior distributions of the source
#' profiles, the source contributions/concentrations, the variance of the data
#' and the other parameters.  Eventually, this function will incorporate 
#' options for concentrations observed below the MDL
#'
#' @param dat data frame of daily constituent concentrations with date as first column
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
#' data(lamcon)
#' #fix data totals
#' pm <- nycdat$PM25
#' whPM <- which(colnames(nycdat) == "PM25")
#' nycdat <- nycdat[, -whPM]
#' whPM <- which(colnames(nycmdl) == "PM25")
#' nycmdl <- nycmdl[, -whPM]
#' L <- 4
#' 
#' mcmcsa(nycdat, L, lamcon)
mcmcsa <- function(dat, lamcon, mdls = NULL, 
	guessvec = NULL, burnin = 10000, N = 100000, fix = NULL){


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
	L <- nrow(lamcon)
	
	
	#if mdl, then censor data
	if(!(is.null(mdls))) {
		bdls <- 1 * (dat < mdls)
		dat[which(bdls == 1, arr.ind = T)] <- NA
	}else{
		bdls <- NULL
	}
	
	
	
	#get indices correspond to id constraints
	cond <- which(!(is.na(lamcon)), arr.ind = T)

	#create arrays for output
	lamstar <- array(dim = c(L, P, N - burnin ))
	lfmat <-  array(dim = c(T1, L, N - burnin))
	# sigma2 <- array(dim = c(P, P, N - burnin))
	sigma2 <- array(dim = c(P, N - burnin))
	mu <- array(dim = c(L, N - burnin))
	xi2 <- array(dim = c(L, N - burnin))
	
	#dimnames
	cons <- colnames(dat)
	days <- seq(1, T1)
	sources <- rownames(lamcon)
	iters <- seq(1, dim(lamstar)[3])
	dimnames(lamstar) <- list(sources, cons, iters)
	dimnames(lfmat) <- list(days, sources, iters)
	dimnames(sigma2) <- list(cons, iters)
	dimnames(mu) <- list(sources, iters)
	dimnames(xi2) <- list(sources, iters)	
	
	names1 <- c("lamstar", "lfmat", "sigma2", "mu", "xi2")
	
	#set first guesses
	if(is.null(guessvec)) {
		
		#truncated normal for lamstar
		lamstarG <- matrix(rtruncnorm(L * P, a = 0), nrow = L)
		lamstarG[cond] <- lamcon[cond]
		
		#MVnormal for logF
		lfmatG <- matrix((rnorm(T1 * L)), nrow = T1)
		# sigma2G <- diag(1, P)
		sigma2G <- rep(.0001, P)		
		# muG <- rep(0, L)
		muG <- c(0.453, 1.916, 1.703, 0.582)
		#xi2G <- rep(1, L)
		xi2G <- c(0.035, 0.002, 0.003, 0.012)
		guessvec <- list(lamstarG, lfmatG, sigma2G, muG, xi2G)
		names(guessvec) <- names1
	}
	
	
	#get guess for log data
	if(!(is.null(mdls))) {
		names1 <- c(names1, "ly")
		ldat <- array(dim = c(T1, P, N - burnin))
		
		#get initial guess
		ldatG <- log(dat)
		for(t in 1 : T1) {
			for(p in 1 : P) {
				if(bdls[t, p] == 1) {
					#start with uniform between 0 and mdls
					ldatG[t, p] <- log(runif(1, 0, mdls[t, p]))
				}
			}
		}
		
		guessvec[["ly"]] <- ldatG
	}else{
		
		guessvec[["ly"]] <- log(dat)
		}

	
	#get lfhist
	lfhist <- list(length = T1)
	for(t in 1 : T1) {
		lfhist[[t]] <- lfmatG[t, ]
	}
	

	#for each iteration (N large)
	for (i in 1 : N) {
		print(i)
		
		# if(i == 53) {browser()}
		
		#update all parameters
		for(j in 1 : length(names1)) {
			out <- gibbsfun(guessvec = guessvec, 
				lamcon = lamcon,
				type = names1[j], lfhist = lfhist, iter = i,
				mdls = mdls, bdls = bdls, fix = fix)
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
			
			if(!(is.null(mdls))) {
				ldat[, , i - burnin] <- guessvec[["ly"]]
			}
		}
			
			

		
	}
	
	#create output
	listout <- list(lamstar, lfmat, sigma2, mu, xi2)
	names(listout) <- names1
	
	if(!(is.null(mdls))) {
		listout[["ly"]] <- ldat
	}
	
	
	listout
}



















##########################################
#update to new guess
##########################################
gibbsfun <- function(guessvec, lamcon, type, lfhist, 
	iter, mdls, bdls, fix = NULL) {
	# print(type)
	if(type == "lamstar" && !(type %in% fix)) {
		guessvec <- lamfun(guessvec, lamcon)
		
	} else if(type == "lfmat" && !(type %in% fix)) {
		
		out <- ffun(guessvec, lfhist, iter)
		guessvec <- out[["guessvec"]]
		lfhist <- out[["lfhist"]]
		
	} else if(type == "mu" && !(type %in% fix)) {
		guessvec <- mufun(guessvec)
		
	} else if(type == "sigma2" && !(type %in% fix)) {
		guessvec <- sig2fun(guessvec)
		
	} else if(type == "xi2" && !(type %in% fix)) {
		guessvec <- xi2fun(guessvec)
		
	} else if(type == "ly" && !(type %in% fix)) {
		guessvec <- yfun(guessvec, mdls, bdls)
		
	}else if (!(type %in% fix)){
		stop("Error: type not recognizeds")	
	}		
		
	out <- list(guessvec, lfhist)
	names(out) <- c("guessvec", "lfhist")
	out
}
	
	
	



	
	
#####
# normal/inv-gamma variance sampling for sigma
# guessvec is list of guesses
sig2fun <- function(guessvec)	{

	ldat <- guessvec[["ly"]]
	
	lfmat <- guessvec[["lfmat"]]
	fmat <- exp(lfmat)
	lamstar <- guessvec[["lamstar"]]
	lambda <- sweep(lamstar, 1, rowSums(lamstar), "/")

	
	#get dimensions
	P <- ncol(ldat)
	T1 <- nrow(ldat)

	#get mean
	mean <- log(fmat %*% lambda)
	
	#get priors
	# pra <- 0.01
	# prb <- 0.01
	# Revise prior for invgamma
	pra <- 2
	prb <- 1
	
	
	#get posterior parameters
	a2 <- pra + T1 / 2
	diffs2 <- (ldat - mean)^2
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
# guessvec is list of guesses
xi2fun <- function(guessvec)	{
		
	lfmat <- guessvec[["lfmat"]]
	mu <- guessvec[["mu"]]
	
	#get dimensions
	L <- ncol(lfmat)
	T1 <- nrow(lfmat)
	
	#set prior values
	# pra <- 0.01
	# prb <- 0.01
	# Revise prior for invgamma
	# pra <- 2
	#prb <- 1
	
	
	#pra <- 2.1
	#pra <- 3
	# prb <- c(0.04, 0.002, 0.003, 0.01)
	
	#based on check_priors
	pra <- 1.3
	prb <- 1.5

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
# guessvec is list of guesses
mufun <- function(guessvec)	{
	

	xi2 <- guessvec[["xi2"]]
	lfmat <- guessvec[["lfmat"]]
	
	#get dimensions
	L <- ncol(lfmat)	
	T1 <- nrow(lfmat)
	
	#set prior values
	prmean <- 0
	#prmean <- c(0.45, 1.9, 1.7, 0.58)
	prvar <- 10
	#prvar <- 1
	
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
		







