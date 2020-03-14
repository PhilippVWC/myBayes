#' @title Transform shape parameters of Beta distribution
#' @description This function transforms the shape paramters alpha and beta
#' of a Beta distribution to a "relative standard deviation"
#' @param alpha double - shape parameter >0
#' @param beta double - shape paramter >0
#' @return the resulting "relative standard deviation"
#' @author Philipp van Wickevoort Crommelin
#' @export
sig_alphaBeta = function(alpha,beta) 1/sqrt(alpha+beta+1)
