#shifts every element of a vector >>offset<< times to the right
#' @export
shift = function(vec,offset){
  return(c(vec[(1+offset):length(vec)],vec[1:offset]))
}
