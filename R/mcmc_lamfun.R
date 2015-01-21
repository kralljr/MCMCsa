
#####
# random walk for lambda
#guessvec is current guesses
lamfun <- function(guessvec, lamcon)	{
	
	
	#set old guess
	lamstar <- guessvec[["lamstar"]]
	
	#set dimensions
	L <- nrow(lamstar)
	P <- ncol(lamstar)
	
	#set variance for proposal distribution
	sd1 <- 0.01
	
	#propose
	newlam <- matrix(rnorm(L * P, mean = lamstar, sd = sd1), 
		nrow = L)
	
	for(l in 1 : L) {
		for(p in 1 : P) {
			#propose
			if(newlam[l, p] >= 0 & is.na(lamcon[l, p])) {
				
				lamstarNEW <- lamstar
				lamstarNEW[l, p] <- newlam[l, p]
				
				#compute lhoods
				lhoodOLD <- loglam(lamstar, guessvec, l, p)
				lhoodNEW <- loglam(lamstarNEW, guessvec, l, p)
				
				#decide whether to accept
				lprob <- min(0, lhoodNEW - lhoodOLD)
				lprob <- ifelse(is.na(lprob), -Inf, lprob)
				unif1 <- runif(1)
				#sel <- rbinom(1, 1, prob)
	
				if(log(unif1) <= lprob) {
					#browser()
					lamstar <- lamstarNEW
					guessvec[["lamstar"]] <- lamstar
				}
			
			}

		}
	}
	
	guessvec[["lamstar"]] <- lamstar
	
	guessvec
}



#function to get lhood for lambda
#lamstar is current guess
# guess is other current guesses
loglam <- function(lamstar, guessvec, l, p) {
	
	#update guess
	guessvec[["lamstar"]] <- lamstar
	
	#likelihood of data
	ly <- logly(guessvec, p = p)
	
	# #set prior values from Nikolov 2011
	mn <- -0.5
	sd1 <- sqrt(0.588)
	lam <- dlnorm(lamstar[l, p], meanlog = mn, sdlog = sd1, log = F)
	
	
	# #get truncated normal density
	# mn <- 0.8
	# sd1 <- 0.7
	# lam <- dtruncnorm(lamstar[l, p], a = 0, mean = mn, sd = sd1)
	
	
	

	
	ly +  log(lam)
}

