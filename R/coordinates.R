
# Takes X and Y coordinates and creates a grid in the form of a vector
# X and Y are both 1D, surface is 2D
# output is a matrix with 3 columns containing the coordinates X,Y,Z as column vectors
#' @export
coordinates = function(X,Y,surface){
  N_X = length(X)
  N_Y = length(Y)
  Z = c(surface)
  X = rep(X,N_Y)
  Y = c(t(matrix(Y) %*% rep(1,N_X)))
  return(matrix(c(X,Y,Z),ncol=3))
}
