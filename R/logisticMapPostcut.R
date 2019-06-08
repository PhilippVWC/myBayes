#### nach Iteration gekappte Variante
#' @export
logisticMapPostcut = function(x,r){
  return(min(1,max(0,4*r*x*(1-x))))
}
