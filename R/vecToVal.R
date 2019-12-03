#' @export
#' @title Projection of a vector onto a given axis
#' @param vec vector of type double
#' @param x_lower double - lower bound
#' @param x_upper double - upper bound
#' @return the inner product of the given vector vec and the given axis.
#' @author Philipp van Wickevoort Crommelin
vecToVal = function(vec,x_lower=0,x_upper=1){
  N = length(vec)
  dx = (x_upper-x_lower)/(N-1)
  ref = seq(from = x_lower,
            to = x_upper,
            by = dx)
  return(ref%*%vec)
}
