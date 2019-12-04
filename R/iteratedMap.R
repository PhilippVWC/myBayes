#' @title Time series creation
#' @param x0 double - starting value
#' @param r double - bifurcation parameter
#' @param n integer - length of resulting time series
#' @param funct function of starting value x0 - the corresponding map
#' @details This routine takes a map an applies it iteratively
#' @return vector of type double and length n, i.e. the resulting time series.
#' @examples
#' r0 = 0.918999
#' x0 = 0.5
#' N = 100
#' y = iteratedMap(x = x0,
#'                 r = r0,
#'                 n = N,
#'                 funct = function(x,r) myBayes::gsm_cpp(x,
#'                                                        r,
#'                                                        alpha = 2.0))
#' plot(x = 1:N,
#'      y = y,
#'      ylim = c(0,1),
#'      type = "l")
#' @export
iteratedMap = function(x,
                       r,
                       n,
                       funct){
  S = x
  if(n>1){
    for (i in 2:n){
      S = c(S,funct(S[i-1],r))
    }
  }
  return(S)
}
