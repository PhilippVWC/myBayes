# Gaussian white chaotic map (S. Sprott: Chaos and Time-Series Analysis, S. 420)
#' @export
gaussianWhiteChaoticMap = function(x,r){
  return(r*erfinv(1-2*erf(x/r)))
}
