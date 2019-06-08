#### allgemeines dynamisches Parameterrauschen
#' @export
noisyParmMap = function(x,r
                        ,funct  = logisticMap
                        ,prfunct = function(z) rnorm(n=1,mean=z,sd=SD0) ){
  return(funct(x,prfunct(r)))
}
# Examples:
#noisyParmMap(x0,r0)
#noisyParmMap(x0,r0,prfunct = function(z) rnorm(n=1,mean=z,sd=1))
#noisyParmMap(x0,r0,prfunct = function(z) min(1,max(0,rnorm(n=1,mean=z,sd=1))))
