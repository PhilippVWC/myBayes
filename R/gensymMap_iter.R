# general symmetric map iteration and time series generation

#' @export
gensymMap_iter = function(N,x0,r,alpha,N_discr = 0,skipFirst=TRUE){
  X = rep(x0,N)
  if(skipFirst){
    X[1] = if(N_discr != 0){
      round(gensymMap(x = x0,
                      r = r,
                      alpha = alpha)*N_discr)/N_discr
    }else{
      gensymMap(x = x0,
                r = r,
                alpha = alpha)
    }
  }
  if(N>1){
    if(N_discr != 0){
      for (i in 2:N) X[i] = round(gensymMap(x = X[i-1],
                                            r = r,
                                            alpha = alpha)*N_discr)/N_discr
    }else{
      for (i in 2:N) X[i] = gensymMap(x = X[i-1],
                                      r = r,
                                      alpha = alpha)
    }
  }
  return(X)
}
