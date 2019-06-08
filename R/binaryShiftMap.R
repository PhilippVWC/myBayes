#### Binary shift map
#' @export
binaryShiftMap = function(x,r){
  return((2*r*x)%%1)
}

#plotBifurcation (funct=binaryShiftMap, xStart=c(0.2,0.4), nrLevels=100, nPlot = nLevels, nStart=nLevels)
#plotIteratedMaps(mapNames=list("binaryShiftMap"),r=r0,nList=1:5)
#plotTimeSeries(mapName=list("binaryShiftMap"),x=0.2,r=1,n=200)
#plotTimeSeries(mapName=list("binaryShiftMap"),x=0.2,r=1,n=20,cobweb=TRUE)
