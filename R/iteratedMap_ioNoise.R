#' @title Time series with input and measurement noise
#' @details Iterated map with input and output noise
#' @param x0 double - starting value of timeseries
#' @param r double - bifurcation parameter
#' @param n integer - length of time series
#' @param funct function of two scalar values - map
#' @param irfunct function of scalar value - function for input noise
#' @param orfunct function of scalar value - output noise function
#' @return the noisy time series
#' @examples
#' x0 = 0.7
#' r = 0.8
#' N = 100
#' y = iteratedMap_ioNoise(x0 = x0,
#'                         r = r,
#'                         n = N,
#'                         funct=function(x,r) myBayes::gsm_cpp(x,r,alpha = 1.7),
#'                         orfunct = function(u) min(1,max(0,rnorm(n = 1,
#'                                                                 mean = u,
#'                                                                 sd = 0.1))))
#' plot(x = 1:N,
#'      y = y,
#'      type = "l",
#'      ylim = c(0,1))
#'
#' @export
iteratedMap_ioNoise = function(x0
                               ,r
                               ,n
                               ,funct = function(x,r) myBayes::gsm_cpp(x = x0,r = r,alpha = 2)
                               ,irfunct = identity # input noise function(u) rnorm(n=1,mean=u,sd=SD0)
                               ,orfunct = identity # output noise function(u) rnorm(n=1,mean=u,sd=SD0)
){
  y = myBayes::iteratedMap(x = irfunct(x0),
                           r = r,
                           n = n,
                           funct = funct)
  return(sapply(X = y,
                FUN = orfunct))
}
