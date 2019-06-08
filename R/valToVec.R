#' @export
valToVec = function(val,N=10,x_lower=0,x_upper=1){
  if (val > x_upper | val < x_lower ) print("warning: val is not within bounds of x_lower and x_upper")
  if (N>1){
    dx = (x_upper-x_lower)/(N-1)
    ref = seq(from = x_lower,
              to = x_upper,
              by = dx)
    diff = abs(ref - rep(val,N))
    minInd = which(diff == min(diff))
    result = rep(0,N)
    result[minInd[1]] = 1 #if value lies exactly between to grid points -> condense to left one
    return(result)
  }else{
    print(paste0("N needs to be greater than 1. N is currently set to ",N))
    return(NA)
  }
}
