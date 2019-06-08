#' @export
Heaviside = function(x,x0){
  if (x<x0) return(0)
  else return(1)
}
