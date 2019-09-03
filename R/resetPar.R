#' @export
resetPar = function(){
  pdf()
  freshPar = par(no.readonly = TRUE)
  dev.off()
  return(freshPar)
}
