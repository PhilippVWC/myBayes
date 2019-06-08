#' @export
CeBePa = function(mean,sigmaRel){
  N = 1/sigmaRel^2 - 1.0
  alpha = max(1, mean*N)
  beta = max(1, (1-mean)*N)
  return(c(alpha,beta))
}
