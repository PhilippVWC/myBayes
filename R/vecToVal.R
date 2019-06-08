#' @export
vecToVal = function(vec,x_lower=0,x_upper=1){
  N = length(vec)
  dx = (x_upper-x_lower)/(N-1)
  ref = seq(from = x_lower,
            to = x_upper,
            by = dx)
  #s = sum(vec)
  #if (s!=1) print(paste0("Warning, vector is not normed to one. Sum = ",s))
  return(ref%*%vec)
}
