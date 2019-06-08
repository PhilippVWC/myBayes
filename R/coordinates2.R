#' @export
coordinates2 = function(X,Y,surface){
  N_X = length(X)
  N_Y = length(Y)
  X = matrix(rep(X,N_Y),nrow =  N_X)
  Y = c(t(matrix(Y) %*% rep(1,N_X)))
  return(matrix(c(X,Y,Z),ncol=3))
}
