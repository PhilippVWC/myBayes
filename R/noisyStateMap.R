# probabilistische Einschritt-Abbildungen
# allgemeines dynamisches Zustandsrauschen

#' @export
noisyStateMap = function(x,r
                         ,funct  = logisticMap
                         ,rfunct = function(z) rnorm(n=1,mean=z,sd=SD0) ){
  return(rfunct(funct(x,r)))
}
# Examples:
#noisyStateMap(x0,r0)
#noisyStateMap(x0,r0,rfunct = function(z) rnorm(n=1,mean=z,sd=1)) # Bsp. Normalverteilung
#noisyStateMap(x0,r0,rfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=1)))) # Bsp. Normalverteilung mit Postcut
#noisyStateMap(x0,r0,rfunct = function(z) rbeta_MuSDrel(n=1, x=x0, sdRel=SDrel0)) # Bsp. Betaverteilung
#noisyStateMap(x0,r0,rfunct = rbeta_MuSDrel0) # Bsp. Betaverteilung mit festem sdRel=SDrel0
