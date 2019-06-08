#' @export
discretizeMap = function(N,map,xlim){
  #--- parameters ---#
  #N      number of discrete points (equal to range of matrix)
  #map    map to be discretized. Should be reduced to function of one variable
  #xlim   vector with two entries containing the bounds of the problem
  x_lower = xlim[1]
  x_upper = xlim[2]
  dx = (x_upper-x_lower)/(N-1)
  domain = seq(from = x_lower,
               to = x_upper,
               by = dx)
  codomain = sapply(X = domain,
                    FUN = function(x) map(x))

  helperMatrix = matrix(codomain,
                        ncol = 1) %*% rep(1,N)

  # "condense" function values to grid points
  indices = apply(X = helperMatrix,
                  MARGIN = 1,
                  FUN = function(x) {
                    dist = abs(x - domain) # distance from grid points
                    index = which(dist == min(dist)) # index of
                    return(index)
                  })
  A = matrix(rep(0,N*N),
             ncol = N,
             nrow = N)
  for (i in 1:N){
    A[indices[i],i] = 1
  }
  return(A)
}
