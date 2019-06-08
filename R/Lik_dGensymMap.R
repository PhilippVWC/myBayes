#' @export
Lik_dGensymMap = function(alpha,r,x0,Y,sigma,N_discr){
  n = length(Y)
  X = dGensymMap_iter(N = n,
                      x0 = x0,
                      r = r,
                      alpha = alpha,
                      N_discr = N_discr,
                      skipFirst = TRUE)
  L = exp(LLik(X = X,
               Y = Y,
               sigma = sigma)/(sqrt(2*pi)*sigma))
  return(L)
}
