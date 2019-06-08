# dynamisches Rauschen (intrinsisches Rauschen)
# logistische Abbildung mit dynamischem (additiven, betaverteiltem) Rauschen
# (0<x0<1, 0<a<4, 0 < sig)

#' @export
logMap_dyn = function(N,a,x0,sigmaRel=0.1){
  X = c(x0)
  if(N>=1){
    for (i in 1:N){
      if( (X[i] > 0) & (X[i] < 1) & (sigmaRel > 0) ) {
        #Var = max(0,min(sigma^2,x*(1-x))) # Varianz aus absolutem Sigma, muss bei maximaler Varianz gekappt werden
        #Var = x*(1-x)*(sigmaRel^2) # Varianz aus relativem Sigma
        #nu = x*(1-x)/Var - 1.0 # "sample size" aus Mittelwert (mu=x) und Varianz = sigma^2
        nu = 1/sigmaRel^2 - 1.0 # "sample size"
        alpha = max(1e-10, X[i] * nu) # alpha aus sample size (nu) und Mittelwert (mu=x)
        beta = max(1e-10, (1-X[i]) * nu) # beta aus sample size (nu) und Mittelwert (mu=x)
        X[i] = rbeta(1,alpha,beta) # neues x mit Rauschen um mu=x
        #cat("i=",i,"sigma=",sigma,"nu=",nu,"alpha=",alpha, "beta=",beta," x=",x,"\n")
      }
      X[i+1] = a*X[i]*(1-X[i])
    }
  }
  return(X[-1])
}
