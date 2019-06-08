#' @export
Heaviside = function(x,x0){
  if (x<x0) return(0)
  else return(1)
}
#' @export
rectFun = function(x,a,b,z){
  result = z * (Heaviside(x = x,
                          x0 = a) - Heaviside(x = x,
                                              x0 = b))
  return(result)
}
