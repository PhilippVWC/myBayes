# Zeitreihe
# (eindimensionaler, homogener, diskreter hidden Markovprozess)

#' @export
oneDimhHMM = function(x,r,n
                      ,funct=logisticMap # function(x,r)
                      ,irfunct = identity # function(x)
                      ,orfunct = identity# function(x)
){
  y=rep(NA,n)
  x = irfunct(x)
  if(n>=1){
    for (i in 1:n){
      x = funct(x,r)
      y[i] = orfunct(x)
    }
  }
  return(y)
}
# Examples:
#oneDimhHMM(x0,r0,10)
#plot(oneDimhHMM(x0,0.999,100),type="l",col="blue",lwd=1)
#lines(oneDimhHMM(x0,0.999,100,funct=discreteLogisticMap),type="l",col="orange",lty=3)
#lines(oneDimhHMM(x0,0.999,100,funct=tentMap),type="l",col="red",lwd=1,lty=2)
#lines(oneDimhHMM(x0,0.999,100,funct=lorenzMap),type="l",col="green",lwd=2)
