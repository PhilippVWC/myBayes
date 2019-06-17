#' @export
LikM = function(mat,x0,Y,sigma=0.1,skipFirst = TRUE){
  d = length(Y)
  X = discreteMap_iter(N = d,
                       x0 = X0,
                       mat = mat,
                       skipFirst = skipFirst)
  return(exp(-sum( 0.5*((X-Y)/sigma)^2 ))/(sqrt(2*pi)*sigma)^d )
}
