#### vor Iteration gekappte Variante
#' @export
logisticMapPrecut = function(x,r){
  x=min(1,max(0,x))
  return(min(1,max(0,4*r*x*(1-x))))
}
