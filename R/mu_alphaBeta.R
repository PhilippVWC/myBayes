#' @title Transform shape parameters of Beta distribution
#' @description This function transforms the shape paramters alpha and beta
#' of a Beta distribution to a mean value.
#' @param alpha double - shape parameter >0
#' @param beta double - shape paramter >0
#' @return the resulting mean value
#' @author Philipp van Wickevoort Crommelin
#' @export
mu_alphaBeta = function(alpha,beta) 1/(1+(beta/alpha))

