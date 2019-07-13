#' @title The Beta distribution
#' @description The beta distribution is a continuous propability distribution defined in the interval [0,1].
#' @param x Double - A value within the intervall [0,1].
#' @param alpha Double - A value greater than zero.
#' @param beta Double - A value greater than zero.
#' @return Double - The corresponding probability.
#' @author J.C. Lemm, P.v.W. Crommelin
#' @export
betaDist = function(x,alpha,beta){
  if (x<=1 & x>=0 & alpha > 0 & beta > 0) {
    B = gamma(alpha+beta)/(gamma(alpha)*gamma(beta))
    return(B * x^(alpha-1)*(1-x)^(beta-1))
  }
}
