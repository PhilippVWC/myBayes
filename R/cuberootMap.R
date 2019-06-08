#' @export
cuberootMap = function(x,r){
	# (0 <= x0 <= 1, 0 <= r <= 1)
  return(r*(1-(2*abs(x-0.5))^(1/3)))
}
