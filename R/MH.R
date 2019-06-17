#' @title Metropolis Hastings Algorithm
#' @param T double - temperature
#' @param x0 double - starting value
#' @param trKern function of two values - Transition kernel as a function from two double values
#' @param sampler sampling distribution
#' @param N integer - size of time series
#' @param bounds array of size two - bounds for the search in 1D
#' @return The resulting Markov Chain as an array of size N
#' @author P. van Wickevoort Crommelin
#' @description Metropolis Hastings Algoritm. Create a random variable, which is distributed according to a custom probability distribution.
#' @export
MH = function(T = 1,x0,trKern=dnorm,sampler,N=1000,bounds){
  a = bounds[1] # left bound of domain to sample from
  b = bounds[2] # right bound
  chain = rep(0,N) # Markov chain
  chain[1] = x0
  for (i in 2:N){
    x = chain[i-1]
    y = x + T*runif(n = 1,
                    min = -(x-a),
                    max = b-x) # sample next candidate within bounds
    A = min( 1 , (sampler(y)*trKern(y,x) / (sampler(x)*trKern(x,y)) ) )
    if (A == 1) { chain[i] = y }
    else {
      if ( A > runif(n = 1,
                     min = 0,
                     max = 1)){
        chain[i] = y
      }
      else {
        chain[i] = x
      }
    }
  }
  return(chain)
}
