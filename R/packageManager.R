#' @export
packageManager = function(necessaryPackages){
  installedPackages = installed.packages()[,"Package"]
  missingPackages = necessaryPackages[!(necessaryPackages %in% installedPackages)]
  if (length(missingPackages) > 0){
    cat("Installation of the following packages:", missingPackages ,"\n")
    cat("This may take a while\n")
    install.packages(missingPackages)
    installedPackages = installed.packages()[,"Package"]
  }
  successfullyInstalled = missingPackages[missingPackages %in% installedPackages]
  notInstalled = missingPackages[!(missingPackages %in% installedPackages)]
  if (length(successfullyInstalled) > 0){
    cat("The following packages were installed:",successfullyInstalled,"\n")
  }
  if (length(notInstalled) > 0) {
    cat("The following packages were not installed:",notInstalled,"\n")
  }
  lapply(necessaryPackages,require,character.only=TRUE) # load library after installation
  return("packageManager finished")
}
