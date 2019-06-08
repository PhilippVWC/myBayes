#' @export
blurrVector = function(vec,lambda,sigmaRel = 0.1){
  n = length(vec) - 1
  mu = which(vec == max(vec)) - 1
  N_samp = 1/sigmaRel^2 - 1
  y_s = min(1/(N_samp*mu/n),1/(N_samp*(1-mu/n))) # ensures that (alpha,beta) == 1 <==> lambda == 0
  lambda_t = 1-(1-y_s)*lambda #transform lambda
  alpha = max(1,lambda_t * mu/n * N_samp)
  beta = max(1,lambda_t * (1-mu/n) * N_samp)
  x = 0:n
  result = sapply(X = x,
                  FUN = function(x) (1-lambda)*kroneckerDelta(x = x,
                                                              x0 = mu) + lambda*betaBin(x = x,
                                                                                        alpha = alpha,
                                                                                        beta = beta,
                                                                                        n = n))
  return(result)
}
