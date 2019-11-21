#' @export
resetPar = function(){
  png()
  freshPar = par(no.readonly = TRUE)
  dev.off()
  return(freshPar)
}
