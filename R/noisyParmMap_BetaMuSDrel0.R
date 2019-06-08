# Spezialfall: dynamisches Parameterrauschen mit Betaverteilung

#' @export
noisyParmMap_BetaMuSDrel0 = function(x,r,funct=logisticMap) noisyParmMap(x, r, funct, prfunct=rbeta_MuSDrel0)
# Examples:
#noisyParmMap_BetaMuSDrel0(x0,r0)
#noisyParmMap_BetaMuSDrel0(x0,r=1)
