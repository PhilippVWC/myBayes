#### allgemeines dynamisches (Zustands- und Parameter-)Rauschen
# (ersetzt noisyStateMap und noisyParmMap)
#' @export
noisyMap = function(x,r
                    ,funct  = logisticMap  # definiert deterministische univariate Abbildung
                    ,rfunct = function(z) rnorm(n=1,mean=z,sd=SD0)  # definiert dynamisches Zustandsrauschen nach Iterationsschritt
                    ,prfunct = function(z) rnorm(n=1,mean=z,sd=SD0) # definiert dynamisches Parameterrauschen
                    ,xrfunct = identity # definiert dynamisches Zustandsrauschen vor Iterationsschritt
){
  #  return(rfunct(funct(x,prfunct(r))))
  return(rfunct(funct(xrfunct(x),prfunct(r)))) # Variante, bei auch auch x vor der Iteration durch xrfunct() verrauscht wird
}
# Examples:
#noisyMap(x0,r0)
#noisyMap(x0,r=r0,funct=logisticMap)
#noisyMap(x0,r=r0,funct=tentMap)
