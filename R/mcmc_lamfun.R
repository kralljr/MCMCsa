
#####
# random walk for lambda?
lamfun <- function(dat, guessvec)	{
	
	lamstar <- guessvec[["lamstar"]]
	sd1 <- 0.1
	
	for(l in 1 : L) {
		for(p in 1 : P) {
			#propose
			newlam <- rnorm(1, lamstar[l, p], sd = sd1)
			lamstarNEW <- lamstar
			lamstarNEW[l, p] <- newlam
			
			#compute lhoods
			lhoodOLD <- lhoodlam(lamstar, dat, guessvec)
			lhoodNEW <- lhoodlam(lamstarNEW, dat, guessvec)
			
			#decide whether to accept
			prob <- min(1, lhoodNEW / lhoodOLD)
			sel <- rbinom(1, 1, prob)
			lamstar <- ifelse(sel == 1, lamstarNEW, lamstar)

		}
	}
	
	guessvec[["lamstar"]] <- lamstar
	
	guessvec
}




lhoodlam <- function(lamstar, dat, guessvec) {
	
	#likelihood of data
	guessvec[["lamstar"]] <- lamstar
	first <- prod(ly(dat, guessvec))
	
	mn <- -0.5
	sd <- sqrt(0.588)
	second <- dtruncnorm(lamstar, a = 0, mean = mn, sd = sd)
	
	first * second
}

