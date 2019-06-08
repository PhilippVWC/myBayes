nLevels = 100
#' @export
discreteRandomMap = function(x,r){
	#Drvals = pmax(0,ceiling(runif(nLevels)*nLevels)-1)/(nLevels-1) # Levels von 0 bis nLevels-1 jew. zwischen 0 und 1
	Drvals = pmax(0,ceiling(runif(nLevels)*nLevels)-1) # Levels von 0 bis nLevels-1
	Drvals = Drvals - min(Drvals) # Starte immer bei 0
  #  0 < x <= 1, 0 <= r <= 1
  return(ceiling(r*Drvals[pmax(1,ceiling(x*nLevels))])/(nLevels-1))
}
