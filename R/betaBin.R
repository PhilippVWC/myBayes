#' @title Beta binomial distribution
#' @description The Beta Binomial distribution is a discrete probability distribution.
#' @param x Double - A value within the interval [0,1].
#' @param alpha Double - A value greater than zero.
#' @param beta Double - A value greater than zero.
#' @param n Integer - Number of Bernoulli trials.
#' @return Double - The corresponding probability.
#' @author J.C. Lemm, P.v.W. Crommelin
#' @export
betaBin = function(x,alpha,beta,n){  # Beta binomial distribution
  choose(n,x) * beta(alpha+x,beta+n-x) / beta(alpha,beta)
}
