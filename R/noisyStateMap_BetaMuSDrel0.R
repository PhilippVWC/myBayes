#### Spezialfall: dynamisches Zustandsrauschen mit Betaverteilung
#' @export
noisyStateMap_BetaMuSDrel0 = function(x,r,funct=logisticMap) noisyStateMap(x, r, funct, rfunct=rbeta_MuSDrel0)
# Examples:
#noisyStateMap_BetaMuSDrel0(x0,r0)
#noisyStateMap_BetaMuSDrel0(x0,r=1)
