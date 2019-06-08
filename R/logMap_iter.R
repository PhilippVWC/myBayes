# logistic map
#' @export
logMap_iter = function(N,r,x0,skipFirst=TRUE){
  if (N >= 1) {
    X=rep(0.0,N)
    if(skipFirst) X[1] = 4*r*x0*(1-x0)
    else X[1] = x0
    for (n in 2:N){
      X[n] = 4*r*X[n-1]*(1-X[n-1])
    }
    return(X)
  } else {
    cat("invalid number N =",N,"\n")
    return(NULL)
  }
}
