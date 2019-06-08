# einfache Iteration
# (auch geeignet fuer dynamische Rauschvarianten)

#' @export
iteratedMap = function(x,r,n
                       ,funct=logisticMap ){
  if(n>=1){
    for (i in 1:n){
      x = funct(x,r)
    }
  }
  return(x)
}

# Examples:
#iteratedMap(x0, r0, 2)
#iteratedMap(x0, r0, 2, funct = logisticMap)
#iteratedMap(x0, r0, 2, funct = function(x,r) noisyMap(x, r, funct=logisticMap))
#iteratedMap(x0, r0, 2, funct = noisyMap_BetaMuSDrel0)
