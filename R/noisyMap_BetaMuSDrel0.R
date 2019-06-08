#### Spezialfall mit Betaverteilung
# (mit vorgegebener SD f?r dynamisches Zustands- und Parameterrauschen)
#' @export
noisyMap_BetaMuSDrel0 = function(x,r,funct=logisticMap)
  noisyMap(x,r,funct
           , rfunct = rbeta_MuSDrel0
           ,prfunct = rbeta_MuSDrel0
  )
# Examples:
#noisyMap_BetaMuSDrel0(x0,r0)
#noisyMap_BetaMuSDrel0(x0,r0,funct=tentMap)
