#' @title Parametrize shape parameter alpha for beta distribution
#' @description This simple function computes alpha as a function of a relative standard deviation sdrel
#' and a mean value mu
#' @param mu double - mean value within (0,1)
#' @param sdrel double - relative standard deviation
#' @return double - The shape parameter alpha
#' @references Kruschke, John K. (2011). Doing Bayesian data analysis: A tutuorial with R and BUGS
#' @export
alpha_muSDrel = function(mu,sdrel){
  return(max(1e-10,mu*myBayes::nu(sdrel)))
}
