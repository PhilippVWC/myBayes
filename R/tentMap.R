#### tent map (Zelt-Abbildung)
# (0 <= x0 <= 1, 0 <= r <= 1)
#' @export
tentMap = function(x,r){
  return(r*(1-2*abs(x-0.5)))
}
