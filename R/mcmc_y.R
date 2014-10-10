# yfun <- function(dat, guessvec) {
	# #set old estimates
	# gsig <- guessvec[[3]]		
	# gthet <- guessvec[[2]]
	# gdat <- guessvec[[1]]
	
	
	# #set up for univariate
	# m <- 1

	
	# #get inverse of gsig
	# gsiginv <- chol2inv(chol(gsig))
	
	# #which columns have missing
	# wcolmiss <- which(colSums(nbdls) < nrow(nbdls))
	# gsig.misscol <- array(dim = c(nrow(gsig) - 1, 
		# ncol(gsig) - 1, ncol(gsig)))
		
	# #for each column with missingness	
	# for( i in 1 : length(wcolmiss)) {
		# num <- wcolmiss[i]
		# sigcut <- gsiginv[-num, -num] - gsiginv[num, -num] %*% 
			# t(gsiginv[-num, num]) / gsiginv[num, num]
		# #update array	
		# gsig.misscol[, , num] <- sigcut
	# }
	
	
	# # for each day with an observation below the MDL
	# for (k in 1 : nrow(dat)) {
		
		# #set up mean
		# mnzero <- t(gdat[k, ] - gthet)

		# #if which missing greater than 0
		# if(length(wh_miss[[k]]) > 0) {
		
			# #for each missing value on day i
			# for (j in 1 : length(wh_miss[[k]])) {
		
				# #find conditional mean/var
				# mnv <- impmisssingle(datk = gdat[k, ], 
					# gthet = gthet, 
					# gsig = gsig, gsig.misscol = gsig.misscol, 
					# wh = wh_miss[[k]][j], mnzero = mnzero)
	
				# #propose truncated normal (0 to mdl)	 
					# #log scale with cond mean and var		
				# newymiss1 <- rtnorm(1, lower = minmdls - 10, 
					# upper = mdls[k, wh_miss[[k]][j]],
					# mean = mnv[[1]], sd = sqrt(mnv[[2]]))
				# if( is.na(newymiss1)) {browser()}	
		
				# #update guess
				# guessvec[[1]][k, wh_miss[[k]][j]] <- newymiss1
		
				
	
# }