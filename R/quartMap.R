# quartische Abbildung
# (0 <= x0 <= 1, 0 <= r <= 1)

#' @export
quartMap = function(x,r){
  return(r*(1-16*(x-0.5)^4))
  # r*(1-(2*(x-1/2))^4) =  r*(1-16*(x-0.5)^4)
  # = r*(1-16*(x^2 -x +1/4)^2) = r*(1 - 16*(x^4 +x^2 +1/16 -2*x^3 +x^2/2 -x/2))
  # = r*(1 -1 -2^3*x -3*2^3*x^2 -2^5*x^3 -2^4*x^4)
  # = 2^3*r*(x -3*x^2 -2^2*x^3 -2*x^4) = 8*r*(x -3*x^2 +4*x^3 -2*x^4)
  # = 8*r*(x -x^2 -2*x^2 +2*x^3 +2*x^3 -2*x^4)
  # = 8*r*x*(1-x)*(1-2*x+2*x^2) (hier ist die Nullstelle bei x=0 und x=1 direkt abzulesen)
  # = 8*r*(x-x^2)*(1-2*x+2*x^2)
  # = logisticMap(x,r)*2*(1-2*x+2*x^2)
}
