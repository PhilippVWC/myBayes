#' @export
LLik = function(X,Y,sigma){
  # strict result = -sum( 0.5*((Y-X)/sigma)^2 ) + length(X)*log(1/(sqrt(2*pi)*sigma))
  return(-sum( 0.5*((Y-X)/sigma)^2 ))
}
