#### Spezialfall mit gekappter Normalverteilung
# (mit vorgegebener SD f?r dynamisches Zustands- und Parameterrauschen)
#' @export
noisyMap_NormalCut0 = function(x,r,funct=logisticMap)
  noisyMap(x,r,funct
           , rfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=SD0)))
           ,prfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=SD0)))
  )
# Examples:
#noisyMap_NormalCut0(x0,r0)
#noisyMap_NormalCut0(x0,r0,funct=tentMap)
