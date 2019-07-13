#' @title Remove equal and adjacent components of a vector
#' @description This functions takes a vector and removes duplicate, neighbouring entries. Indices are preserved under
#' this operation. If you encounter a sequence of equal values within a vector, the first one is preserved and its right
#' neighbours are removed.
#' @param vec Vector of type double - The input vector.
#' @param order Boolean - Value indicating wether resulting dataframe should be sorted in decreasing order (defaults to TRUE).
#' @return A data frame with two columns. 1st column $values: The resulting vector. 2nd column $indices:
#' The corresponding indices matching the input vector.
#' @author P.v.W. Crommelin
#' @examples
#' vec = c(1,1,3,4,5,6,7,7,7,7,7,8,9,10,9,8,7,6,5,4,4,4,4,3,2,1)
#' N = length(vec)
#' plot(x = 1:N,
#'      y = vec,
#'      type = "l",
#'      main = "demonstration")
#' points(x = 1:N,
#'        y = vec)
#' collapsedVec = collapse(vec = vec,
#'                         order = FALSE)
#' lines(x = collapsedVec$indices,
#'       y = collapsedVec$values,
#'       col = "red")
#' points(x = collapsedVec$indices,
#'        y = collapsedVec$values,
#'        col = "red")
#' legend("topright",
#'        leg = c("original vector","collapsed vector"),
#'        col = c("black","red"),
#'        lty = 1)
#' @export
collapse = function(vec,order = TRUE){
  deltaVec = vec - c(vec[1]+1,vec[1:(length(vec)-1)])
  indices = which(deltaVec != 0)
  df = data.frame(values = vec[indices], indices = indices)
  if(order){
    df = df[order(df$values,
                  decreasing = TRUE),]
  }
  return(df)
}
