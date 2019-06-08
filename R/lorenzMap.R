#### Lorenz-Abbildung (Wurzel)
# (0 <= x0 <= 1, 0 <= r <= 1)

#' @export
lorenzMap = function(x,r){
  return(r*(1-sqrt(abs(x-0.5))*sqrt(2))) # = r*(1-(2*abs(x-1/2))^(1/2))
}
