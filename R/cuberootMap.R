#### Kubische Wurzel-Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
cuberootMap = function(x,r){
  return(r*(1-(2*abs(x-0.5))^(1/3)))
}
