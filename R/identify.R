#' Check identifiability conditions
#' 
#' \code{identifySA} checks a identifiability conditions for source apportionment
#' 
#' This function checks conditions C1- C3-1 from Park et al. 2002 Environmetrics
#' 
#' @param lamcon Matrix of number of sources by number of constituents containing the 0s and 1s for identifiability conditions
#' @export
identifySA <- function(lamcon) {
  L <- nrow(lamcon)
  
  # Check C1
  c1 <- apply(lamcon, 1, function(x) length(which(x == 0)))
  if(min(c1) != (L - 1)) {
    stop(cat("Error: condition C1 not met for source", which.min(c1), "\n"))
  }

  # Check condition C2

  tf <- vector()
  for(i in 1 : L) {
    wh0 <- which(lamcon[i, ] == 0)  
    lam1 <- lamcon[-i, wh0]
    whNA <- which(is.na(lam1))
    N <- length(whNA)
    lam1[whNA] <- rnorm(N)
    tf[i] <- 1 * (qr(lam1)$rank == (L - 1))
  }
  diffs <- max(tf) - min(tf)
  if(diffs != 0 | min(tf) != 1) {
    stop(cat("Error: condition C2 not met for sources:", which(tf != (L - 1)), "\n"))
  }

  # Check condition C3
  c3 <- apply(lamcon, 1, function(x) length(which(x == 1)))
  diffs <- max(c3) - min(c3)
  if(diffs != 0 | min(c3) != 1) {
    stop(cat("Error: condition C1 not met for source", which.min(c1), "\n"))
  }

  cat("Conditions C1-C3 are met\n")

}

