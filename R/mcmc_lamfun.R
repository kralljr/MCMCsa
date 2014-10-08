
#####
# random walk for lambda
#dat is data
#guessvec is current guesses
lamfun <- function(dat, guessvec)	{
	
	
	#set old guess
	lamstar <- guessvec[["lamstar"]]
	
	#set dimensions
	L <- nrow(lamstar)
	P <- ncol(lamstar)
	
	#set variance for proposal distribution
	sd1 <- 0.1
	
	#propose
	newlam <- rnorm(L * P, mean = lamstar, sd = sd1)
	
	for(l in 1 : L) {
		for(p in 1 : P) {
			#propose
			lamstarNEW <- lamstar
			lamstarNEW[l, p] <- newlam[l, p]
			
			#compute lhoods
			lhoodOLD <- lhoodlam(lamstar, dat, guessvec, p)
			lhoodNEW <- lhoodlam(lamstarNEW, dat, guessvec, p)
			
			#decide whether to accept
			prob <- min(1, lhoodNEW / lhoodOLD)
			sel <- rbinom(1, 1, prob)
			lamstar <- ifelse(sel == 1, lamstarNEW, lamstar)

		}
	}
	
	guessvec[["lamstar"]] <- lamstar
	
	guessvec
}



#function to get lhood for lambda
#lamstar is current guess
#dat is data
# guess is other current guesses
lhoodlam <- function(lamstar, dat, guessvec, p) {
	
	#update guess
	guessvec[["lamstar"]] <- lamstar
	
	#likelihood of data
	first <- ly(dat, guessvec, p)
	
	#set prior values
	mn <- -0.5
	sd <- sqrt(0.588)
	
	#get truncated normal density
	second <- dtruncnorm(lamstar, a = 0, mean = mn, sd = sd)
	
	first * second
}

