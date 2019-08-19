#' @title Definiteness of a 2D real-valued matrix
#' @description This routine checks the definiteness of a two dimensional, real-valued matrix.
#' @param x matrix of type double - A real-valued, 2D matrix.
#' @return character string - A character string indicating the definiteness of the given matrix.
#' @details This routine takes the symmetric part of the given matrix and
#' decomposes its spectrum to find an answer.
#' @author Philipp van Wickevoort Crommelin
#' @examples
#' # positive definite
#' A = matrix(data = c(1,0,0,0,7,0,0,0,2),
#'            ncol = 3)
#' definiteness(A)
#'
#' # negative definite
#' definiteness(-A)
#'
#' # positive semi-definite
#' A = matrix(data = c(1,0,0,0,0,0,0,0,2),
#'            ncol = 3)
#' definiteness(A)
#'
#' # negative semi-definite
#' definiteness(-A)
#'
#' # Example, that demonstrates the importance of symmetry of the given matrix:
#' B = matrix(data = c(4,1,9,4),
#'            ncol = 2) #Eigenvalues all positive but not symmetric
#' definiteness(B)
#' eigen(B)$values
#' # Scalar product maps to negative value
#' v = c(-1,1)
#' t(v)%*%B%*%v
#'
#' @export
definiteness = function(x){
  A = 0.5*(x+t(x)) #Get symmetric part of x
  if (min(eigen(x = A)$values)>0) return("positive definite")
  if (min(eigen(x = A)$values)>=0) return("positive semi-definite")
  if (max(eigen(x = A)$values)<0) return("negative definite")
  if (max(eigen(x = A)$values)<=0) return("negative semi-definite")
  return("indefinite")
}
