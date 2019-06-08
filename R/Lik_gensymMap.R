#' @export
Lik_gensymMap = function(alpha,r,x0,Y,sigma,N_discr = 0){
  # N_discr only for compatibility reasons
  n = length(Y)
  X = gensymMap_iter(N = n,
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
