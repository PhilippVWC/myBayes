#' @title Parametrize shape parameter beta for beta distribution
#' @description This simple function computes beta as a function of a relative standard deviation sdrel
#' and a mean value mu
#' @param mu double - mean value within (0,1)
#' @param sdrel double - relative standard deviation
#' @return double - The shape parameter beta
#' @references Kruschke, John K. (2011). Doing Bayesian data analysis: A tutuorial with R and BUGS
#' @export
beta_muSDrel = function(mu,sdrel){
  return(max(1e-10,(1-mu)*myBayes::nu(sdrel)))
}
