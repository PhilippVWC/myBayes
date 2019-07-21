#' @export
toPDF = function(file = "plot.pdf",
                 openPDF = FALSE){
  dev.copy2pdf(file = file,
               out.type = "pdf")
  if(openPDF) browseURL(file)
}
