#' @export
betaDist = function(x,alpha,beta){
  if (x<=1 & x>=0 & alpha > 0 & beta > 0) {
    B = gamma(alpha+beta)/(gamma(alpha)*gamma(beta))
    return(B * x^(alpha-1)*(1-x)^(beta-1))
  }
}
