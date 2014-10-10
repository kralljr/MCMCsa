
#####
# random walk for lambda
#dat is data
#guessvec is current guesses
lamfun <- function(dat, guessvec, lamcon)	{
	
	
	#set old guess
	lamstar <- guessvec[["lamstar"]]
	
	#set dimensions
	L <- nrow(lamstar)
	P <- ncol(lamstar)
	
	#set variance for proposal distribution
	sd1 <- 1
	
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
				lhoodOLD <- loglam(lamstar, dat, guessvec, p)
				lhoodNEW <- loglam(lamstarNEW, dat, guessvec, p)
				
				#decide whether to accept
				lprob <- min(0, lhoodNEW - lhoodOLD)
				lprob <- ifelse(is.na(lprob), 0, lprob)
				unif1 <- runif(1)
				#sel <- rbinom(1, 1, prob)
	
				if(log(unif1) <= lprob) {
					lamstar <- lamstarNEW
				}
			
			}

		}
	}
	
	guessvec[["lamstar"]] <- lamstar
	
	guessvec
}



#function to get lhood for lambda
#lamstar is current guess
#dat is data
# guess is other current guesses
loglam <- function(lamstar, dat, guessvec, p) {
	
	#update guess
	guessvec[["lamstar"]] <- lamstar
	
	#likelihood of data
	ly <- logly(dat, guessvec, p = p)
	
	#set prior values
	mn <- -0.5
	sd <- sqrt(0.588)
	
	#get truncated normal density
	lam <- dtruncnorm(lamstar, a = 0, mean = mn, sd = sd)
	
	ly +  log(lam)
}

