#### Iteration mit Inputrauschen und abschlie?endem Beobachtungsrauschen (Ausgangsrauschen)
#' @export
iteratedMap_ioNoise = function( x
                                ,r
                                ,n
                                ,funct  = logisticMap
                                ,irfunct = identity # input noise function(u) rnorm(n=1,mean=u,sd=SD0)
                                ,orfunct = identity # output noise function(u) rnorm(n=1,mean=u,sd=SD0)
){
  #  return(rfunct(funct(x,r)))
  return(orfunct(iteratedMap(irfunct(x),r,n,funct)))
}
# Examples:
#iteratedMap_ioNoise(x0,r0,1)
#iteratedMap_ioNoise(x0,r0,1,irfunct=function(u) min(1,max(0,rnorm(n=1,mean=u,sd=0.1))))
#iteratedMap_ioNoise(x0,r0,1,irfunct=function(u) rbeta_MuSDrel(n=1,mean=u,sdRel=0.1))
#iteratedMap_ioNoise(x0,r0,1,orfunct=function(u) rnorm(n=1,mean=u,sd=0.1))
#iteratedMap_ioNoise(x0,r0,1,funct=function(y,r) iteratedMap(y,r,n=10,funct=logisticMap))
#iteratedMap_ioNoise(x0,r0,10,funct=logisticMap)
#iteratedMap_ioNoise(x0,r0,10,funct=tentMap)
