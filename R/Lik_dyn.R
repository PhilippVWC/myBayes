#Likelihood function for dynamic noise (beta distributed)
#returns a modell Vector X of dimension Dim(Y)-1
#with Y being the data vector, respectively
#' @export
Lik_dyn = function(a,x0,Y,sigmaRel){
  ly = length(Y)
  n = ly-1
  X = a*Y[2:ly]*(1-Y[2:ly])
  Lik = exp(-sum( ( 0.5*(X-Y[2:ly])/(sigmaRel) )^2  ))/(sigmaRel*sqrt(2*pi))^n
  return(Lik)
}
