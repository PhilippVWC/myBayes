#' @title  blurrVector statistically disperses a deterministic vector
#' @description This function is a convex combination of the Kronecker Delta function and the Betabinomial distribution.
#' @param vec Vector of type integer - A deterministic vector containing only zeros, except for one exposed component containing a one.
#' @param lambda A control parameter within [0,1].
#' @param sigmaRel Double - This value controlls the width of the maximum and eventually the dispersion of the resulting vector.
#' @return A vector of type double - A probability vector, that is statistically dispersed around the exposed index position and normed to one.
#' @author P.v.W. Crommelin
#' @examples
#' N = 150
#' deterministicVec = rep(0,N)
#' deterministicVec[N%/%4] = 1
#' probabilisticVec = blurrVector(vec = deterministicVec,
#'                                lambda = 0.1,
#'                                sigmaRel = 0.1)
#' plot(x = NULL,
#'      y = NULL,
#'      xlim = c(1,N),
#'      ylim = c(0,1),
#'      type = "p",
#'      cex = 3,
#'      col = "black",
#'      pch = 16)
#' M = 5
#' lambda = seq(from = 0,
#'              to = 1,
#'              length.out = M)
#' col = rainbow(M)
#' sapply(X = 1:M,
#'        FUN = function(i){
#'          points(x = 1:N,
#'                y = blurrVector(vec = deterministicVec,
#'                                lambda = lambda[i],
#'                                sigmaRel = 0.1),
#'                col = col[i])
#'        })
#' #' @export
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
