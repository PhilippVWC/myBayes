#### kubische Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
cubeMap = function(x,r){
  return(r*(1-8*abs(x-0.5)^3)) # = r*(1-(2*abs(x-1/2))^3)
}
