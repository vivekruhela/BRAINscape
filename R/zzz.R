.onLoad <- function(libname, pkgname) {
  # List of required packages
  required_packages <- c("shiny", "ggplot2", "dplyr", "CMplot")
  
  # Install missing packages
  new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
  
  if(length(new_packages)) {
    install.packages(new_packages, dependencies = TRUE)
  }
  
  # Install Bioconductor dependencies
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  bioc_packages <- c("DESeq2", "limma")
  
  for (pkg in bioc_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
  }
}
