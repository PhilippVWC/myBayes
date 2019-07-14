#' @title Simulated annealing in one dimension
#' @description Simulated annealing is a Monte Carlo method for (statistical) optimization purposes.
#' The underlying Metropolis Hastings Algoritm is a Markov Chain Monte Carlo method
#' for creation of a series of random variables.
#' Every random variable is distributed according to a supplied probability distribution (the sampler) and the
#' temperature T controls the jumping width within each iteration.
#' @param T Double - Temperature.
#' @param x0 Double - Starting value.
#' @param trKern Function of two values - Transition kernel as a function of two values of type double. Defaults
#' to the normal distribution (Equal probabilities for jumps in both directions.)
#' @param sampler The sampling distribution from which to draw the sample as function of one value of type double.
#' @param N Integer - Size of the resulting time Markov Chain.
#' @param bounds Two-component vector of type double - Bounds for the search in one dimension.
#' @return The resulting Markov Chain as an array of size N.
#' @author Philipp van Wickevoort Crommelin
#' @examples
#' sampler = function(x){
#'   # Sampler as a sum of five Gaussian distributions.
#'   x0 = c(2,6,-1,-7,4)
#'   sigma = c(1,0.3,0.7,0.2,0.1)
#'   result = 0
#'   l = length(x0)
#'   for (i in 1:length(x0)){
#'     result = result + dnorm(x = x,
#'                             mean = x0[i],
#'                             sd = sigma[i])
#'   }
#'   return(result/l)
#' }
#'
#' trKern = function(x,y){
#'   SIGMA = 1
#'   return(dnorm(x = x,
#'                mean = y,
#'                sd = SIGMA))
#' }
#' N_x = 1000
#' x_lower = -10
#' x_upper = 10
#' dx = (x_upper-x_lower)/(N_x-1)
#' x = seq(from = x_lower,
#'         to = x_upper,
#'         by = dx)
#' y = sampler(x)
#' maxY = max(y)
#' scale = maxY/N_x
#' plot(x = NULL,
#'      y = NULL,
#'      xlim = c(min(x),max(x)),
#'      ylim = c(0,maxY),
#'      main = "Simmulated Annealing",
#'      xlab = "x",
#'      ylab = "sampler")
#' lines(x = x,
#'       y = y,
#'       col = "blue")
#' N_steps = 1000
#' T = .5
#' x0 = 2
#' MCMC = MH(T = T,
#'           x0 = x0,
#'           trKern = trKern,
#'           sampler = sampler,
#'           N = N_steps,
#'           bounds = c(x_lower,x_upper))
#' lines(x = rep(x0,2),
#'       y = c(0,maxY),
#'       col = "grey",
#'       lty = 2,
#'       lwd = 2)
#' lines(x = MCMC,
#'       y = seq(from = 0,
#'               to = maxY,
#'               by = maxY/(N_steps-1)),
#'       col = "green",
#'       lwd = 0.5)
#' FinalVal = MCMC[N_steps]
#' lines(x = rep(FinalVal,2),
#'       y = c(0,maxY),
#'       col = "red",
#'       lty = 1,
#'       lwd = 3)
#'
#' legend("topleft",
#'        leg = c("Sampling distribution","Markov chain","Starting value","Final value"),
#'        col = c("blue","green","grey","red"),
#'        lty = c(1,1,2,1),
#'        lwd = c(1,0.5,2,3))
#' @export
MH = function(T = 1,x0,trKern=dnorm,sampler,N=1000,bounds){
  a = bounds[1] # Left bound of domain to sample from.
  b = bounds[2] # Right bound.
  chain = rep(0,N) # Markov chain.
  chain[1] = x0
  dt = T/N
  t = seq(from = T,
             to = dt,
             by = -dt)
  for (i in 2:N){
    x = chain[i-1]
    y = x + t[i]*runif(n = 1,
                    min = -(x-a),
                    max = b-x) # Sample next candidate within bounds.
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
