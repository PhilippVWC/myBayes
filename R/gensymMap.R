# general symmetric map (S. Sprott: Chaos and Time-Series Analysis, S. 39)
# (0 <= x0 <= 1, 0 <= r <= 1, logisticMap: alpha = 2, lorenzMap: alpha = 1/2 etc. )

#' @title General symmetric map
#' @description a generalization of a one dimensional map, which incorporates (among others) the lorenz-, logistic and cubeMap
#' @param x input value for one iteration
#' @author J.C. Lemm, P.v.W. Crommelin
#' @references S. Sprott, Chaos and Time-series analysis
#' @param r control parameter
#' @param alpha exponent
#' @return The result of one single iteration
#' @export
gensymMap = function(x,r,alpha){
  return(r*(1-(2*abs(x-0.5))^alpha))
}
