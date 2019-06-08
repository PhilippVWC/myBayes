#' @export
discreteMap_iter = function(N,x0,mat,skipFirst = TRUE){
  x0 = valToVec(val = x0,
                N = dim(mat)[2])
  if (skipFirst) x0 = mat %*% x0
  Series = matrix(data = rep(x0,N),
                  nrow = length(x0),
                  ncol = N,
                  byrow = FALSE)
  for (i in 2:N){
    Series[,i] = mat %*% Series[,i-1]
  }
  result = apply(X = Series,
                 MARGIN = 2,
                 FUN = function(colVec){
                   vecToVal(vec = colVec)
                 })
  return(result)
}
