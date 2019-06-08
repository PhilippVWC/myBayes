#not a likelihood function
#' @export
Lik = function(X,Y,sigma){
  return(-LLik(X = X,
               Y = Y,
               sigma = sigma))
}
