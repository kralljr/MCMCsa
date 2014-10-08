#' Bayesian source apportionment model
#'
#' Estimates sources of particulate matter (PM) air pollution using
#' chemical constituent concentrations of PM.  The model is developed
#' using Nikolov et al. 2008 and Hackstadt and Peng 2014.  
#' 
#' To estimate sources, users need data from an ambient monitor, which is 
#' a matrix of days by number of chemical constituents.  Then, users apply
#' the function \code{\link{mcmcsa}} to the data, specifying the number
#' of sources and necessary conditions for identifiability.
#'
#' @references Margaret C. Nikolov, Brent A. Coull, Paul J. Catalano, et al.
#  (2008).  J. R. Statist. Soc. C.  
#' Statistical methods to evaluate health effects associated with 
#' major sources of air pollution: a case-study of breathing 
#' patterns during exposure to concentrated Boston air 
#' particles, 57(3) 357-378. 
#' @references Amber J. Hackstadt and Roger D. Peng (2014).  Environmetrics.
#' A Bayesian multivariate receptor model for estimating source contributions 
#' to particulate matter pollution using national databases, 25(7) 513-527.
#' @name MCMCsa
#' @docType package
#' @import MCMCpack
#' @import truncnorm
#' @import mvtnorm
NULL
