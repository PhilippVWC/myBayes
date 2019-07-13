#' @title Generate coordinates
#' @description This function takes two vectors X and Y that define a domain and a surface Z,
#' and creates a matrix containing the coordinates.
#' @param X Vector of type integer - x domain
#' @param Y Vector of type integer - y domain
#' @description  X and Y are both 1D, surface is 2D
#' @return output is a matrix with 3 columns containing the coordinates X,Y,Z as column vectors
#' @export
coordinates = function(X,Y,surface){
  N_X = length(X)
  N_Y = length(Y)
  Z = c(surface)
  X = rep(X,N_Y)
  Y = c(t(matrix(Y) %*% rep(1,N_X)))
  return(matrix(c(X,Y,Z),ncol=3))
}
