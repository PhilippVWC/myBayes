#' @export
rectFun = function(x,a,b,z){
  result = z * (Heaviside(x = x,
                          x0 = a) - Heaviside(x = x,
                                              x0 = b))
  return(result)
}
