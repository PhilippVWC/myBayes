#' @export
discreteLogisticMap =function(x,r,m=100){
  #  0 < x <= 1, 0 <= r <= 1
  return(round(logisticMap(round(x*(m-1))/(m-1),r)*(m-1))/(m-1))
}
