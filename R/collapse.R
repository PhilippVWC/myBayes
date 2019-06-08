#' @export
collapse = function(vec){
  deltaVec = vec - c(vec[1]+1,vec[1:(length(vec)-1)])
  indices = which(deltaVec != 0)
  df = data.frame(values = vec[indices], indices = indices)
  return(df[order(df$values,
                  decreasing = TRUE),])
}
