#' @export
discr = function(val,N_discr,lower = 0.0,upper = 1.0){
  dx = (upper-lower)/(N_discr-1)
  ref = seq(from = lower,
            to = upper,
            by = dx)
  diff = abs(ref - rep(val,N_discr))
  factor = which( diff == min(diff) )[1] # if value is exactly between two grid points --> condense to left one
  return((factor-1)*dx)
}
