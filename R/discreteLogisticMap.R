#### discrete logistic map
# Tabelle diskrete logistische Abbildung (bei festem r so schneller)
# Dlogvals = round(logisticMap((0:(nLevels-1))/(nLevels-1),r=r0)*(nLevels-1))
# --- Diskrete Funktion
#' @export
discreteLogisticMap =function(x,r,m=nLevels){
  #  0 < x <= 1, 0 <= r <= 1
  #return(Dlogvals[round(x*nLevels)]/nLevels)
  return(round(logisticMap(round(x*(m-1))/(m-1),r)*(m-1))/(m-1))
}
#discreteLogisticMap(x=seq(0,1,length.out = nLevels),r0)
