#' @export
betaBin = function(x,alpha,beta,n){  # Beta binomial distribution
  choose(n,x) * beta(alpha+x,beta+n-x) / beta(alpha,beta)
}
